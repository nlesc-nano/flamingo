"""Utility functions."""
from pathlib import Path
from subprocess import getoutput
from typing import Any, Dict, Iterator, List, Optional

import h5py
import numpy as np
import pandas as pd
from more_itertools import chunked
from rdkit import Chem

__all__ = ["Options", "convert_to_standard_representation", "normalize_smiles",
           "read_molecules", "read_smile_and_sanitize"]


class Options(dict):
    """Extend the base class dictionary with a '.' notation.

    example:
    .. code-block:: python
       d = Options({'a': 1})
       d['a'] # 1
       d.a    # 1
    """

    def __init__(self, *args, **kwargs):
        """ Create a recursive Options object"""
        super().__init__(*args, **kwargs)
        for k, v in self.items():
            if isinstance(v, dict):
                self[k] = Options(v)

    def __getattr__(self, attr):
        """ Allow `obj.key` notation"""
        return self.get(attr)

    def __setattr__(self, key, value):
        """ Allow `obj.key = new_value` notation"""
        self.__setitem__(key, value)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to a normal dictionary."""
        def converter(var):
            return var.to_dict() if isinstance(var, Options) else var

        return {k: converter(v) for k, v in self.items()}


def retrieve_hdf5_data(path_hdf5: Path, paths_to_prop: str) -> np.ndarray:
    """Read Numerical properties from ``paths_hdf5``.

    Parameters
    ----------
    path_hdf5
        path to the HDF5
    path_to_prop
        str or list of str to data

    Returns
    -------
    np.ndarray
        array or list of array

    Raises
    ------
    RuntimeError
        The property has not been found

    """
    try:
        with h5py.File(path_hdf5, 'r') as f5:
            return f5[paths_to_prop][()]
    except KeyError:
        msg = f"There is not {paths_to_prop} stored in the HDF5\n"
        raise KeyError(msg)
    except OSError:
        msg = f"there is no {path_hdf5} file!"
        raise OSError(msg)


def normalize_smiles(smile: str) -> str:
    """Write a smile in its normal form."""
    mol = Chem.MolFromSmiles(smile)
    if mol is not None:
        return Chem.MolToSmiles(mol)
    else:
        return smile


def read_molecules(input_file: Path) -> pd.DataFrame:
    """Read data (e.g. smiles) from a csv-like file."""
    df = pd.read_csv(input_file).reset_index(drop=True)
    # remove unnamed columns
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    return df


def read_molecules_in_batches(input_file: Path, size: int) -> Iterator[Any]:
    """Read a file into chunks."""
    f = open(input_file, 'r')
    # Skip first line with the header
    f.readline()
    return chunked(f.readlines(), size)


def take(it: Iterator[Any], n: int) -> List[Any]:
    """Take n elements of the iterator."""
    return [next(it) for _ in range(n)]


def read_smile_and_sanitize(smile: str) -> Optional[Chem.rdchem.Mol]:
    """Try to read and sanitize a given smile"""
    sanitize = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
    try:
        mol = Chem.MolFromSmiles(smile)
        Chem.rdmolops.SanitizeMol(mol, sanitizeOps=sanitize)
    except:
        mol = None
    return mol


def convert_to_standard_representation(molecules: pd.DataFrame, sanitize_smiles: bool) -> pd.DataFrame:
    """Convert Smiles to standard representation."""
    molecules.dropna(inplace=True)

    # Convert smiles to the standard representation
    if sanitize_smiles:
        back_converter = np.vectorize(mol2smile)
        molecules.smiles = back_converter(molecules.rdkit_molecules)
        molecules.dropna(inplace=True)
    return molecules


def mol2smile(mol: Chem.rdchem.Mol) -> Optional[str]:
    """Return a smile representation if possible."""
    try:
        smile = Chem.MolToSmiles(mol)
    except RuntimeError:
        smile = None
    return smile
