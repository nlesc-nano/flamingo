"""Test the features generation functionality."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from flamingo.features.featurizer import generate_molecular_features


def test_molecular_features():
    """Test that the atomic and bond features are properly created."""
    # Generate mol and add conformers
    mol = Chem.MolFromSmiles("CC(=O)O")
    AllChem.EmbedMolecule(mol)

    atomic, bond = generate_molecular_features(mol)

    # There are four heavy atoms with 17 atomic features each
    assert atomic.shape == (4, 17)

    # There are 3 x 2 bidirectional edges (Bonds) with 7 features each
    assert bond.shape == (6, 7)

    # All entries in the matrix are different of Nan
    assert not np.all(np.isnan(atomic))
    assert not np.all(np.isnan(bond))
