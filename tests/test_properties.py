"""Test the properties calculations with CAT."""

import json
import os
import shutil
from pathlib import Path

from pytest_mock import MockFixture
from scm.plams import Molecule

from flamingo.properties import compute_properties_with_cat

from .utils_test import PATH_TEST


def test_properties_calculation(mocker: MockFixture, tmp_path: Path):
    """Check that the properties are read properly from the HDF5."""
    # Mock call to CAT
    mocker.patch("flamingo.properties.prep", return_value=None)

    smile = "O=C(O)Cc1ccccc1"
    input_file = PATH_TEST / "input_compute_properties_CAT.yml"

    # copyt the hdf5 to the results
    database = tmp_path / "database"
    os.mkdir(database)
    shutil.copy(PATH_TEST / "properties/database/structures.hdf5", database / "structures.hdf5")
    #  Read from the hdf5
    compute_properties_with_cat(smile, input_file, tmp_path.as_posix())

    results = tmp_path / "results.json"
    inputs = tmp_path / "inputs.json"
    geometry = tmp_path / "geometry.xyz"

    # Check that all the files exists
    assert all(p.exists() for p in {results, inputs, geometry})

    # Check that the files can be read
    for file in (results, inputs):
        with open(file, 'r') as handler:
            json.load(handler)

    # read molecule
    Molecule(geometry)
