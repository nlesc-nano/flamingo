"""
Module to compute molecular properties.

API
---
.. autofunction:: compute_properties_with_cat

"""
import argparse
import json
import logging
from contextlib import redirect_stderr
from pathlib import Path
from typing import Any, Dict, List, Union

import h5py
import numpy as np
import pandas as pd
import yaml
from CAT.base import prep
from scm.plams import Settings

from .log_config import configure_logger

logger = logging.getLogger(__name__)


def get_property_from_hdf5(handler: h5py.File, prop: str, name: str) -> Dict[str, Any]:
    """Extract a given property from the hdf5."""
    root = "ligand/properties"
    dset = f"{root}/{prop}"
    if dset not in handler:
        print(f"There is no {dset} dataset in hdf5!")
        return {}

    values = handler[dset][()]
    keys = handler[f"{root}/{name}"][()]
    return {k.decode(): v for k, v in zip(keys, values.flatten())}


def extract_ligand_properties(path_hdf5: Path) -> Dict[str, Any]:
    """Extract the properties as a dictionary of dictionaries."""
    properties = ['E_solv', 'pKa', 'gamma', 'LogP', 'cdft']
    names = [f"{prop}_names" for prop in properties]

    results = {}
    with h5py.File(path_hdf5, 'r') as handler:
        for prop, name in zip(properties, names):
            results[prop] = get_property_from_hdf5(handler, prop, name)

    return results


def read_input_data(data: np.ndarray) -> str:
    """Read a given input ``name``."""
    if data.size == 0:
        return ""

    # Concatenate the data and covert it to str
    data = data.flatten()
    return "\n".join(map(lambda x: x.decode(), data))


def extract_input_files(path_hdf5: Path) -> Dict[str, str]:
    """Get the input used to run the calculations."""
    inputs = {}
    with h5py.File(path_hdf5, 'r') as handler:
        names = [name for name in handler.keys() if name.startswith("job_settings_")]
        for name in names:
            inputs[name] = read_input_data(handler[name][()])

    return inputs


def extract_optimized_geometry(path_hdf5: Path) -> str:
    """Get the optimized geometry from the HDF5."""
    with h5py.File(path_hdf5, 'r') as handler:
        atoms = handler['ligand/atoms'][()]

    data = atoms.flatten()
    data = data[['symbol','x', 'y', 'z']]

    geometry = f"{len(data)}\n\n"
    for symbol, x, y, z in data:
        geometry += f"{symbol.decode()} {x:.8f} {y:.8f} {z:.8f}\n"

    return geometry

def compute_properties_with_cat(
        smile: str, input_file: Union[str, Path], workdir: str) -> None:
    """Compute properties for the given smile and write then down in JSON format.

    Parameters
    ----------
    smile
        String representing the molecular for which the properties are going to be computed
    input_file
        YAML file with the settings to run the simulation
    workdir
        Folder where the simulation is going to write the output

    """
    # Add the smile to the CAT input
    inp = generate_input(smile, input_file)

    path = Path(workdir)

    with open(path / "cat_output.log", 'a') as f:
        with redirect_stderr(f):
            prep(Settings(inp))

    path_hdf5 = path / "database" / "structures.hdf5"

    results = extract_ligand_properties(path_hdf5)
    inputs = extract_input_files(path_hdf5)
    geometry = extract_optimized_geometry(path_hdf5)

    # Store input and results as json
    for data, name in [(results, "results"), (input, "inputs")]:
        with open(path / f"{name}.json", 'w') as handler:
            json.dump(results, handler, indent=4)

    # Store geometry in xyz format
    with open(path / "geometry.xyz", "w") as handler:
        handler.write(geometry)


def generate_input(smile: str, input_file: Union[str, Path]) -> Dict[str, Any]:
    """Insert the smile string into the CAT input."""
    with open(input_file, 'r') as handler:
        inp = handler.read()

    new_input = inp.replace("{smile}", smile)
    return yaml.load(new_input, Loader=yaml.FullLoader)


def main():
    """Parse the command line arguments to screen smiles."""
    parser = argparse.ArgumentParser("compute_properties")
    parser.add_argument('-s', '--smile', required=True, help="smile to compute properties")
    parser.add_argument('-i', '--input', help="YAML settings for CAT")
    parser.add_argument('-w', '--workdir', help="Work directory", default=".")
    args = parser.parse_args()

    # configure logger
    configure_logger(Path("."), "flamingo")

    compute_properties_with_cat(args.smile, args.input, args.workdir)
