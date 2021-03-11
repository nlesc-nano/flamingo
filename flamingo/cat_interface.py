"""Interface with `CAT <https://cat.readthedocs.io/en/latest/>`_.

API
---

.. autofunction:: compute_bulkiness
"""
import logging
import shutil
import tempfile
import uuid
from contextlib import redirect_stderr
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, NamedTuple, Tuple, Union

import h5py
import numpy as np
import pandas as pd
import yaml
from CAT.base import prep
from dataCAT import prop_to_dataframe
from more_itertools import chunked
from nanoCAT.recipes import run_fast_sigma
from scm.plams import Settings

from .utils import Options, normalize_smiles

__all__ = ["compute_bulkiness"]

# Long types
BatchResult = Union[np.ndarray, pd.DataFrame]
Callback = Callable[[pd.Series, Mapping[str, Any], pd.Index], BatchResult]
Reducer = Callable[[List[BatchResult]], BatchResult]

# Starting logger
# logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger('CAT')
logger.propagate = False
handler = logging.FileHandler("cat_output.log")
logger.addHandler(handler)


class PropertyMetadata(NamedTuple):
    """Name and metadata of the property computed by CAT."""

    name: str
    dset: str  # Dset in the HDF5


def call_cat(
        smiles: pd.Series, opts: Mapping[str, Any],
        cat_properties: Dict[str, Any]) -> Tuple[Path, Path]:
    """Call cat with a given `config` and returns a dataframe with the results.

    Parameters
    ----------
    molecules
        Pandas Series with the smiles to compute
    opts
        Options for the computation
    cat_properties
        Dictionary with the name of the properties to compute
    chunk
        Name of the chunk (frame) being computed

    Returns
    -------
    Path to the HDF5 file with the results

    Raises
    ------
    RuntimeError
        If the Cat calculation fails
    """
    # create workdir for cat
    chunk_name = uuid.uuid1().hex
    path_workdir_cat = Path(opts["workdir"]) / "cat_workdir" / chunk_name
    path_workdir_cat.mkdir(parents=True, exist_ok=True)

    path_smiles = (path_workdir_cat / "smiles.txt").absolute().as_posix()

    # Save smiles of the candidates
    smiles.to_csv(path_smiles, index=False, header=False)

    input_cat = yaml.load(f"""
path: {path_workdir_cat.absolute().as_posix()}

input_cores:
    - {Path(opts['core']).absolute().as_posix()}:
        guess_bonds: False

input_ligands:
    - {path_smiles}

optional:
    qd:
       {generate_bulkiness_section(cat_properties)}
    ligand:
       functional_groups:
          ['{opts["anchor"]}']
    database:
        thread_safe: True
""", Loader=yaml.FullLoader)

    with open(path_workdir_cat / "cat_input.yml", 'w') as handler:
        yaml.dump(input_cat, handler)

    inp = Settings(input_cat)
    with open("cat_output.log", 'a') as f:
        with redirect_stderr(f):
            prep(inp)

    path_hdf5 = path_workdir_cat / "database" / "structures.hdf5"

    if not path_hdf5.exists():
        raise RuntimeError(f"There is not hdf5 file at:{path_hdf5}")
    else:
        return path_hdf5, path_workdir_cat


def generate_bulkiness_section(cat_properties: Dict[str, Any]) -> str:
    """Generate the CAT bulkiness input section."""
    def replace_None(x):
        return "NULL" if x is None else x
    if "bulkiness" not in cat_properties:
        return "bulkiness: False"
    bulkiness = cat_properties['bulkiness']
    string = "bulkiness:\n"
    for key in {"h_lim", "d"}:
        string += f"{' ':>10}{key}: {replace_None(bulkiness[key])}\n"
    return string


def compute_property_using_cat(
        smiles: pd.Series, opts: Mapping[str, Any], metadata: PropertyMetadata) -> pd.Series:
    """Compute the bulkiness for the candidates."""
    # Properties to compute using cat
    cat_properties = opts['filters']

    # run cat
    path_hdf5, path_workdir_cat = call_cat(smiles, opts, cat_properties)

    df = extract_dataframe_from_hdf5(path_hdf5, metadata)
    if path_workdir_cat.exists():
        shutil.rmtree(path_workdir_cat)

    return df


def extract_dataframe_from_hdf5(path_hdf5: Path, metadata: PropertyMetadata) -> pd.DataFrame:
    """Get a Dataframe from the CAT HDF5 results."""
    with h5py.File(path_hdf5, 'r') as f:
        dset = f[metadata.dset]
        df = prop_to_dataframe(dset)

    # flat the dataframe and remove duplicates
    df.reset_index(inplace=True)

    # make anchor atom neutral to compare with the original
    # TODO make it more general
    df.ligand = df.ligand.str.replace("[O-]", "O", regex=False)

    # remove duplicates
    df.drop_duplicates(subset=['ligand'], keep='first', inplace=True)

    return df


def compute_batch_bulkiness(
        smiles: pd.Series, opts: Mapping[str, Any], indices: pd.Index) -> pd.Series:
    """Compute bulkiness using CAT."""
    chunk = smiles[indices]

    # Transform the smiles to standard representation
    chunk = chunk.apply(normalize_smiles)

    # compute and extract the bulkiness
    metadata = PropertyMetadata("bulkiness", 'qd/properties/V_bulk')
    df = compute_property_using_cat(chunk, opts, metadata)

    bulkiness = pd.merge(chunk, df, left_on="smiles", right_on="ligand")["V_bulk"]
    if len(indices) != len(bulkiness):
        msg = "There is an incongruence in the bulkiness computed by CAT!"
        logger.error(f"There was an error processing chunk:\n{msg}")
        values = np.repeat(np.nan, len(indices))
    else:
        values = bulkiness.to_numpy()
    return values


def map_reduce(smiles: pd.Series, opts: Options,
               callback: Callback, reduce: Reducer) -> BatchResult:
    """Distribute the properties computation in batches."""
    worker = partial(callback, smiles, opts.to_dict())

    with Pool() as p:
        results = list(p.imap_unordered(worker, chunked(smiles.index, 10), 1))

    return reduce(results)


def compute_bulkiness(smiles: pd.Series, opts: Options) -> np.ndarray:
    """Compute a ligand/quantum dot bulkiness using CAT.

    It creates several instances of CAT using multiprocessing.

    Parameters
    ----------
    smiles
        `pandas.Series` with the smiles to compute
    opts
        Options to call CAT

    Returns
    -------
    numpy.ndarray
        Array with the computed properties
    """
    results = map_reduce(smiles, opts, compute_batch_bulkiness, np.concatenate)

    if len(smiles.index) != results.size:
        msg = "There is an incongruence in the bulkiness computed by CAT!"
        raise RuntimeError(msg)

    return results


def compute_cosmo_rs(
        molecules: pd.DataFrame, solvents: Dict[str, str], workdir: str) -> pd.DataFrame:
    """Compute Cosmo Rs properties using CAT.

    Parameters
    ----------
    molecules
        Pandas.DataFrame with the smiles to compute
    solvents
        Dictionary with the Paths to the solvents data
    workdir
        Directory to write the temporal results

    Returns
    -------
    pandas.DataFrame 
        Values

    example:
    >>> smiles = pd.Series(['CO'])
    >>> solvents = {"hexane": "$AMSRESOURCES/ADFCRS/Hexane.coskf",
                    "toluene": "$AMSRESOURCES/ADFCRS/Toluene.coskf"}
    >>> compute_cosmo_rs(smiles, solvents, Path("."))

    """
    try:
        with tempfile.TemporaryDirectory(prefix="cosmo_rs_", dir=workdir) as output_dir:
            rs = run_fast_sigma(molecules.smiles, solvents, output_dir=output_dir, return_df=True)
            molecules = pd.merge(molecules, rs, left_on="smiles", right_index=True)
    except RuntimeError:
        pass

    return molecules