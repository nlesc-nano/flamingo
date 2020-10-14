"""Check that the CAT library is called correctly."""

from pathlib import Path

import numpy as np
import pandas as pd

from flamingo.cat_interface import compute_batch_bulkiness
from flamingo.utils import read_molecules

from .utils_test import PATH_TEST


def test_compute_bulkiness(tmp_path: Path):
    """Check that computation of bulkiness using CAT."""
    df = read_molecules(PATH_TEST / "smiles_carboxylic.csv")
    opts = {
        "workdir": tmp_path,
        "core": PATH_TEST / "Cd68Se55.xyz",
        "anchor": "O(C=O)[H]"}

    data = compute_batch_bulkiness(df.smiles, opts, pd.Index(list(range(4))))
    assert not np.all(np.isnan(data))
