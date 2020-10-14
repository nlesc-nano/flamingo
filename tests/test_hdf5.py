"""Check HDF5 interface."""

import numpy as np
import pytest

from flamingo.utils import retrieve_hdf5_data

from .utils_test import PATH_TEST

path_cat_hdf5 = PATH_TEST / "output_cat.hdf5"
PATH_BULKINESS = "qd/properties/V_bulk"


def test_hdf5_interface():
    """Check the hdf5 reading interface."""
    values = retrieve_hdf5_data(path_cat_hdf5, PATH_BULKINESS)

    assert not np.all(np.isnan(values))


def test_nonexisting_dataset():
    """Check that an error is raised if no dataset."""
    with pytest.raises(KeyError):
        retrieve_hdf5_data(path_cat_hdf5, "non/existing/dataset")


def test_nonexisting_hdf5():
    """Check that an error is raised if not hdf5 file."""
    with pytest.raises(OSError):
        retrieve_hdf5_data("nonexistingfile.hdf5", PATH_BULKINESS)
