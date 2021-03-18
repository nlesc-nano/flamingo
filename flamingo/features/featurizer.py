"""Compute the fingerprints of an array of smiles."""

from itertools import chain

import numpy as np
import pandas as pd
from rdkit.Chem import AllChem

__all__ = ["generate_fingerprints"]


#: Floating point used to stored the features
DTYPE = np.float32

dictionary_functions = {
    "morgan": AllChem.GetMorganFingerprintAsBitVect,
    "atompair": AllChem.GetHashedAtomPairFingerprintAsBitVect,
    "torsion": AllChem.GetHashedTopologicalTorsionFingerprintAsBitVect
}


def generate_fingerprints(molecules: pd.Series, fingerprint: str, bits: int,
                          use_chirality: bool = False) -> np.ndarray:
    """Generate the Extended-Connectivity Fingerprints (ECFP).

    Available fingerprints:
    * morgan https://doi.org/10.1021/ci100050t
    * atompair
    * torsion
    """
    size = len(molecules)

    it = (compute_fingerprint(molecules[i], fingerprint, bits, use_chirality) for i in molecules.index)
    result = np.fromiter(
        chain.from_iterable(it),
        DTYPE,
        size * bits
    )

    return result.reshape(size, bits)


def compute_fingerprint(molecule, fingerprint: str, nbits: int, use_chirality: bool) -> np.ndarray:
    """Calculate a single fingerprint."""
    # Select the fingerprint calculator
    fingerprint_calculator = dictionary_functions[fingerprint]
    if fingerprint == "morgan":
        bit_vector = fingerprint_calculator(molecule, 2, nBits=nbits, useChirality=use_chirality)
    else:
        bit_vector = fingerprint_calculator(molecule, nBits=nbits, includeChirality=use_chirality)
    return np.fromiter(bit_vector.ToBitString(), DTYPE, nbits)
