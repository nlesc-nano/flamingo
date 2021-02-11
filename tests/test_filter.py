"""Test the screening functionality."""

import argparse
import shutil
from pathlib import Path
from typing import Any, Iterable, Mapping, Set

import numpy as np
import pandas as pd
import pytest
import schema
import yaml
from pytest_mock import MockFixture

from flamingo.screen import main, split_filter_in_batches
from flamingo.utils import Options, read_molecules

from .utils_test import PATH_TEST

PATH_INPUT_TEST_FILTER = PATH_TEST / "input_test_filter.yml"


def run_workflow(opts: Options) -> pd.DataFrame:
    """Apply the filters and read the output."""
    split_filter_in_batches(opts)
    path = Path("results/batch_0/candidates.csv")
    filter_mols = read_molecules(path)
    print("Result path: ", path)
    return filter_mols


def create_options(filters: Mapping[str, Any], smiles_file: str, tmp_path: Path) -> Options:
    """Create Options object to filter."""
    opts = Options()
    opts.smiles_file = (PATH_TEST / smiles_file).absolute().as_posix()
    opts.filters = filters
    opts.output_path = "results"
    opts.workdir = tmp_path
    opts.batch_size = 100
    opts.parallel = np.random.choice([True, False])

    return opts


def remove_output(output_path: str) -> None:
    """Remove the output file if exists."""
    path = Path(output_path)
    if path.exists():
        shutil.rmtree(path)


def check_expected(opts: Options, expected: Set[str]) -> None:
    """Run a filter workflow using `opts` and check the results."""
    try:
        computed = run_workflow(opts)
        print("expected:\n", expected)
        print(f"candidates:\n{computed.smiles.values}")
        assert all(mol in computed.smiles.values for mol in expected)
        assert len(computed.smiles) == len(expected)

    finally:
        remove_output(opts.output_path)


def test_filter_cli(mocker: MockFixture) -> None:
    """Test that the CLI works correctly."""
    mocker.patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(
        i=PATH_INPUT_TEST_FILTER))

    mocker.patch("flamingo.screen.split_filter_in_batches", return_value=None)
    main()


def test_invalid_input(mocker: MockFixture, tmp_path: Path):
    """Check that an error is raised if an invalid input is provided."""
    invalid_input = {"smiles_file": "non-existing", "prop1": "invalid"}
    path_input = tmp_path / "invalid.yml"

    with open(path_input, 'w') as handler:
        yaml.dump(invalid_input, handler)

    mocker.patch("argparse.ArgumentParser.parse_args", return_value=argparse.Namespace(
        i=path_input.absolute().as_posix()))

    with pytest.raises(schema.SchemaMissingKeyError) as info:
        main()

    error = info.value.args[0]
    assert "Missing key" in error


def test_contain_functional_groups(tmp_path: Path) -> None:
    """Test that the functional group filter is applied properly."""
    smiles_file = "smiles_functional_groups.csv"
    filters = {"include_functional_groups": {"groups": ["[CX3](=O)[OX2H1]"],"maximum": 1}}
    opts = create_options(filters, smiles_file, tmp_path)
    expected = {"O=C(O)C1CNC2C3CC4C2N4C13", "C#CC12CC(CO1)NCC2C(=O)O",
                "CCCCCCCCC=CCCCCCCCC(=O)O", "CC(=O)O",
                "O=C(O)Cc1ccccc1", "CC(O)C(=O)O"}
    check_expected(opts, expected)


def test_exclude_functional_groups(tmp_path: Path) -> None:
    """Test that some functional groups are excluded correctly."""

    smiles_file = "smiles_functional_groups.csv"
    filters = {"exclude_functional_groups": {"groups": [
        "[#7][#6](=[OX1])", "C#C", "[#6](=[OX1])[OX2][#6]", "[NX3]"], "maximum": 1}}
    opts = create_options(filters, smiles_file, tmp_path)
    expected = {"c1ccccc1", "CCO", "CCCCCCCCC=CCCCCCCCC(=O)O",
                "CC(=O)O", "O=C(O)Cc1ccccc1", "CC(O)C(=O)O"}
    check_expected(opts, expected)


def test_filter_bulkiness(tmp_path: Path) -> None:
    """Test that the bulkiness filter is applied properly."""
    smiles_file = "smiles_carboxylic.csv"
    filters = {"bulkiness": {"h_lim": None, "d": "auto", "lower_than": 20}}
    opts = create_options(filters, smiles_file, tmp_path)
    opts.core = PATH_TEST / "Cd68Se55.xyz"
    opts.anchor = "O(C=O)[H]"

    expected = {"CC(=O)O", "CC(O)C(=O)O"}
    check_expected(opts, expected)


def test_filter_bulkiness_no_core(tmp_path: Path) -> None:
    """Test that the bulkiness filter is applied properly."""
    smiles_file = "smiles_carboxylic.csv"
    filters = {"bulkiness": {"h_lim": None, "lower_than": 20}}
    opts = create_options(filters, smiles_file, tmp_path)
    opts.anchor = "O(C=O)[H]"

    expected = set()  # type: Set[str]
    with pytest.raises(RuntimeError) as err:
        check_expected(opts, expected)

    error = err.value.args[0]
    print(err)
    assert all(x in error for x in ("bulkiness","core"))


def test_filter_scscore_lower(tmp_path: Path) -> None:
    """Test that the scscore filter is applied properly."""
    smiles_file = "smiles_carboxylic.csv"
    filters = {"scscore": {"lower_than": 1.3}}
    opts = create_options(filters, smiles_file, tmp_path)

    expected = {"CC(=O)O"}
    check_expected(opts, expected)


def test_filter_scscore_greater(tmp_path: Path) -> None:
    """Test that the scscore filter is applied properly."""
    smiles_file = "smiles_functional_groups.csv"
    filters = {"scscore": {"greater_than": 3.0}}
    opts = create_options(filters, smiles_file, tmp_path)

    expected = {"O=C(O)C1CNC2C3CC4C2N4C13"}
    check_expected(opts, expected)


def test_single_carboxylic(tmp_path: Path) -> None:
    """Check that only molecules with a single Carboxylic acids are included."""
    smiles_file = "smiles_carboxylic.csv"
    filters = {"include_functional_groups": {"groups": ["[CX3](=O)[OX2H1]"], "maximum": 1}}
    opts = create_options(filters, smiles_file, tmp_path)
    opts.anchor = "O(C=O)[H]"

    expected = {"CCCCCCCCC=CCCCCCCCC(=O)O", "CC(=O)O", "O=C(O)Cc1ccccc1", "CC(O)C(=O)O"}

    check_expected(opts, expected)


def test_single_functional_group(tmp_path: Path) -> None:
    """Check that molecules with a single functional group are filtered."""
    smiles_file = "smiles_multiple_groups.csv"
    filters = {"include_functional_groups": {
        "groups": ["[CX3](=O)[OX2H1]", "[#16X2H]", "[NX3;H2,H1;!$(NC=O)]"], "maximum": 1}}
    opts = create_options(filters, smiles_file, tmp_path)
    opts.anchor = "O(C=O)[H]"

    expected = {"NCCc1ccncc1", "O=C(O)C1CCC1(F)F"}
    check_expected(opts, expected)


def test_multiple_anchor(tmp_path: Path) -> None:
    """Check that molecules with multiple Carboxylic acids are included."""
    smiles_file = "smiles_carboxylic.csv"
    filters = {"include_functional_groups": {"groups": ["[CX3](=O)[OX2H1]"], "maximum": 2}}
    opts = create_options(filters, smiles_file, tmp_path)
    opts.anchor = "O(C=O)[H]"

    expected = {"CCCCCCCCC=CCCCCCCCC(=O)O", "CC(=O)O", "O=C(O)Cc1ccccc1", "CC(O)C(=O)O", "O=C(O)c1cccc(C(=O)O)c1"}

    check_expected(opts, expected)

