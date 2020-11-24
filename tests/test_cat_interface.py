"""Check that the CAT library is called correctly."""

from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd
from pytest_mock import MockFixture

from flamingo.cat_interface import compute_batch_bulkiness
from flamingo.utils import read_molecules

from .utils_test import PATH_TEST


def create_options(tmp_path: Path) -> Dict[str, Any]:
    """Create options dictionary."""
    return {
        "workdir": tmp_path,
        "core": PATH_TEST / "Cd68Se55.xyz",
        "anchor": "O(C=O)[H]",
        "filters": { "bulkiness": {"h_lim": 10, "d": "auto"}}
        }


def test_compute_bulkiness(tmp_path: Path):
    """Check that computation of bulkiness using CAT."""
    df = read_molecules(PATH_TEST / "smiles_carboxylic.csv")
    opts = create_options(tmp_path)

    data = compute_batch_bulkiness(df.smiles, opts, pd.Index(list(range(4))))
    assert not np.all(np.isnan(data))


def test_wrong_len(mocker: MockFixture, tmp_path: Path):
    """Check that an error is raise if an incomplete result is return."""
    df = read_molecules(PATH_TEST / "smiles_carboxylic.csv")
    opts = create_options(tmp_path)

    mocked_df = pd.DataFrame.from_dict({'ligand': ['CO', 'CC'], 'V_bulk': [1, 2]})
    mocker.patch(
        "flamingo.cat_interface.compute_property_using_cat", return_value=mocked_df)

    data = compute_batch_bulkiness(df.smiles, opts, pd.Index(list(range(4))))

    assert np.all(np.isnan(data))
