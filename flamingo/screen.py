"""
Module to screen smile by functional group and other properties.

API
---
.. autofunction:: split_filter_in_batches
.. autofunction:: apply_filters

"""

import argparse
import logging
import sys
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Callable, Iterator, List, Mapping, Optional, Set, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Fragments

from .cat_interface import compute_bulkiness
from .features.featurizer import generate_fingerprints
from .log_config import configure_logger
from .models.scscore import SCScorer
from .schemas import validate_input
from .utils import (Options, get_number_of_smiles, read_molecules_in_batches,
                    take)

logger = logging.getLogger(__name__)


def split_filter_in_batches(opts: Options) -> None:
    """Split the computations into smaller batches that fit into memory.

    Parameters
    ----------
    opts
        :class:`swan.utils.Options` options to run the filtering

    """
    # Create folder to store the output
    result_path = Path(opts.output_path)
    result_path.mkdir(exist_ok=True, parents=True)

    # Read smiles using an iterator
    smiles_file = Path(opts.smiles_file)
    smiles = read_molecules_in_batches(smiles_file, opts.batch_size)

    # Check precomputed batches
    computed_batches = search_for_computed_batches(opts.output_path)
    if computed_batches > 0:
        logger.info(f"There are already {computed_batches} batches computed!")
        # removed the computed batches
        take(smiles, computed_batches - 1)


    tasks = enumerate(smiles, start=computed_batches)
    if not opts.parallel:
        # Run batches sequential in a single CPU
        for k, batch in tasks:
            compute_batch(opts.to_dict(), result_path, (k, batch))
    else:
        # Run in Multiple CPUs
        worker = partial(compute_batch,opts.to_dict(), result_path)
        with Pool() as p:
            list(p.imap_unordered(worker, tasks, 1))


def compute_batch(opts: Mapping[str, Any], result_path: Path, data: Tuple[int, List[str]]) -> None:
    """Compute a single filtering batch."""
    opts = Options(opts)
    number, smiles = data
    batch = pd.DataFrame({"smiles": map(lambda x: x.rstrip(), smiles)})

    logger.info(f"computing batch: {number}")
    output_file = create_ouput_file(result_path, number)
    try:
        apply_filters(batch, opts, output_file)
    except RuntimeError as err:
        print("Error applying filter:\n", err)
        raise
    except Exception as ex:
        error, msg, _ = sys.exc_info()
        logger.error(f"Error processing batch: {number}\n{error} {msg}", exc_info=ex)


def apply_filters(molecules: pd.DataFrame, opts: Options, output_file: Path) -> None:
    """Apply a set of filters to the given smiles.

    Parameters
    ----------
    molecules
        :class:`pandas.Dataframe` with the molecular data.
    opts
        :class:`flamingo.utils.Options` options to run the filtering
    output_file
        :class:`pathlib.Path`
    """
    logger.debug("converting smiles to rdkit molecules")
    # Create rdkit representations
    converter = np.vectorize(read_smile_and_sanitize)
    molecules["rdkit_molecules"] = converter(molecules.smiles)

    # Remove invalid molecules
    molecules = molecules[molecules.rdkit_molecules.notnull()]

    # Convert smiles to the standard representation
    back_converter = np.vectorize(Chem.MolToSmiles)
    molecules.smiles = back_converter(molecules.rdkit_molecules)

    # Apply all the filters
    available_filters = {
        "include_functional_groups": include_functional_groups,
        "exclude_functional_groups": exclude_functional_groups,
        "bulkiness": filter_by_bulkiness,
        "scscore": filter_by_scscore
        }

    for key in opts.filters.keys():
        if key in available_filters:
            molecules = available_filters[key](molecules, opts)
            if molecules.empty:
                print("There no more molecules to perform the filter in the batch!")

    columns = [x for x in ("smiles", 'scscore', 'bulkiness') if x in molecules.columns]
    molecules.to_csv(output_file, columns=columns)
    logger.debug(f"The filtered candidates has been written to the {output_file} file!")


def include_functional_groups(molecules: pd.DataFrame, opts: Options) -> pd.DataFrame:
    """Check that the molecules contain some functional groups."""
    groups = opts["filters"]["include_functional_groups"]["groups"]
    maximum = opts["filters"]["include_functional_groups"]["maximum"]
    logger.debug(f"including molecules that contains the groups: {groups}")
    if maximum == 1:
        return filter_single_group(molecules, groups, False)
    else:
        return filter_by_functional_group(molecules, groups, False)


def exclude_functional_groups(molecules: pd.DataFrame, opts: Options) -> pd.DataFrame:
    """Check that the molecules do not contain some functional groups."""
    groups = opts["filters"]["exclude_functional_groups"]["groups"]
    logger.debug(f"exclude molecules that contains the groups: {groups}")
    return filter_by_functional_group(molecules, groups, True)


def filter_single_group(
    molecules: pd.DataFrame, groups: List[str], exclude: bool) -> pd.DataFrame:
    """Check that the molecule has a single functional group."""
    # First filter the molecules that have the target functional groups
    molecules = filter_by_functional_group(molecules, groups, False)

    # Now remove the molecules that have more than 2 different functional groups
    patterns = {Chem.MolFromSmarts(f) for f in groups}
    pattern_check = np.vectorize(partial(has_single_substructure, patterns))
    molecules = molecules[pattern_check(molecules["rdkit_molecules"])]

    # Finally exclude molecules that have more that 2 identical functional groups
    return exclude_same_functional_groups(molecules, patterns)

def filter_by_functional_group(
    molecules: pd.DataFrame, groups: List[str], exclude: bool) -> pd.DataFrame:
    """Search for a set of functional_groups."""
    # Transform functional_groups to rkdit molecules
    patterns = {Chem.MolFromSmarts(f) for f in groups}

    # Function to apply predicate
    pattern_check = np.vectorize(partial(has_substructure, patterns))

    # Check if the functional_groups are in the molecules
    has_pattern = pattern_check(molecules["rdkit_molecules"])
    if exclude:
        has_pattern = ~has_pattern

    return molecules[has_pattern]


def exclude_same_functional_groups(molecules: pd.DataFrame, patterns: Set[Chem.rdchem.Mol]) -> pd.DataFrame:
    """Exclude molecules that has the same functional group more than one time."""
    # Iterate over all functional groups
    for counter_function in generate_fragment_counters(patterns):
        function_fragments = np.vectorize(counter_function)
        # Count the number of a specific functional group in a molecule
        number_of_groups = function_fragments(molecules["rdkit_molecules"])
        # Keep only molecules with a single functional group
        molecules = molecules[number_of_groups == 1]
    return molecules


def generate_fragment_counters(patterns: Set[Chem.rdchem.Mol]) -> Iterator[Callable]:
    """Search for a functional group in pattern and return a function that this group in a molecule."""
    functional_groups = (
        ("[OH]C=O", "fr_COO"),          #  Carboxylic acid
        ("CO", "fr_Al_OH"),             #  Aliphatic Hydroxyl
        ("Oc1ccccc1", "fr_Ar_OH"),      #  Aromatic Hydroxyl
        ("C1=CC=C(C=C1)N", "fr_Ar_NH"), #  Aromatic amines
        ("CN", "fr_NH2"),               #  Primary amines
        ("CNC","fr_NH1"),               #  Secondary amines
        ("CP(=O)(O)O", "fr_phos_acid"), #  Phosphonic acid
        ("SC", "fr_SH"),                #  Thiols
        ("CS(=O)(=O)O", "fr_sulfone"),  #  Sulfonic acid
    )
    # Select the functional groups in patterns
    pairs = [(has_substructure(patterns, Chem.MolFromSmiles(smile)), getattr(Fragments, fr))
            for smile, fr in functional_groups]

    # Check that the functional group is present and select its fragment counter
    return map(lambda t: t[1], filter(lambda t: t[0], pairs))


def has_substructure(patterns: Set[Chem.rdchem.Mol], mol: Chem.Mol) -> bool:
    """Check if there is any element of `pattern` in `mol`."""
    return False if mol is None else any(mol.HasSubstructMatch(p) for p in patterns)


def has_single_substructure(patterns: Set, mol: Chem.Mol) -> bool:
    """Check if there a single functional pattern in mol."""
    acc = 0
    for pat in patterns:
        has_pattern = mol.HasSubstructMatch(pat)
        if has_pattern and acc == 0:
            acc += 1
        # There is more than one pattern
        elif has_pattern and acc > 0:
            return False
    return True


def filter_by_bulkiness(molecules: pd.DataFrame, opts: Options) -> pd.DataFrame:
    """Filter the ligands that have a given bulkiness.

    The bulkiness is computed using the CAT library: https://github.com/nlesc-nano/CAT
    The user must specify whether the bulkiness should be lower_than, greater_than
    or equal than a given value.
    """
    logger.info("Filtering by bulkiness")
    if opts.core is None:
        raise RuntimeError("A core molecular geometry is needed to compute bulkiness")

    opts.bulkiness = True
    molecules["bulkiness"] = compute_bulkiness(molecules.smiles, opts)
    logger.debug("CAT has been called!")
    return apply_predicate(molecules, "bulkiness", opts)


def filter_by_scscore(molecules: pd.DataFrame, opts: Options) -> pd.DataFrame:
    """Compute the `SCScore` for `molecules` and filter those that fulfill the predicate."""
    logging.info("filtering by scscore")

    scorer = SCScorer("1024bool")
    fingerprints = generate_fingerprints(molecules.rdkit_molecules, "morgan", 1024, use_chirality=True)
    molecules["scscore"] = scorer.compute_score(fingerprints)

    return apply_predicate(molecules, "scscore", opts)


def apply_predicate(molecules: pd.DataFrame, feature: str, opts: Options) -> pd.Series:
    """Apply `predicate_type` on `column_name`."""

    # Check if the molecules fulfill the bulkiness predicate
    property_info = opts.filters[feature]
    keywords = property_info.keys()

    predicate = "lower_than" if "lower_than" in keywords else "greater_than"
    limit = property_info[predicate]
    if predicate == "lower_than":
        has_pattern = molecules[feature] <= limit
    else:
        has_pattern = molecules[feature] >= limit

    logger.info(f"Keep molecules that have {feature} {predicate} {limit}")

    return molecules[has_pattern]


def create_ouput_file(result_path: Path, k: int) -> Path:
    """Create path to print the resulting candidates."""
    parent = result_path / f"batch_{k}"
    parent.mkdir(exist_ok=True)
    return parent / "candidates.csv"


def read_smile_and_sanitize(smile: str) -> Optional[Chem.rdchem.Mol]:
    """Try to read and sanitize a given smile"""
    sanitize = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
    try:
        mol = Chem.MolFromSmiles(smile)
        Chem.rdmolops.SanitizeMol(mol, sanitizeOps=sanitize)
    except:
        mol = None
    return mol


def merge_result():
    """Merge all the results file into a single one in the CWD."""
    files = [pd.read_csv(p,index_col=0) for p in Path("results").glob("batch_*/candidates.csv")]
    if files:
        results = pd.concat(files)
        results.to_csv("FinalResults.csv", index=False)


def search_for_computed_batches(output_path: Path) -> int:
    """Check for batches that have been already computed."""
    path = Path(output_path)
    if path.exists():
        computed = len(set(path.glob("batch_*")))
        return computed - 1 if computed > 0 else 0
    else:
        return 0


def main():
    """Parse the command line arguments to screen smiles."""
    parser = argparse.ArgumentParser("smiles_screener")
    # configure logger
    parser.add_argument('-i', required=True,
                        help="Input file with options")
    args = parser.parse_args()

    configure_logger(Path("."), "flamingo")

    # parse command line options and run workflow
    options = validate_input(args.i, action="screen")
    split_filter_in_batches(options)
    merge_result()
