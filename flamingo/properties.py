"""
Module to compute molecular properties.

API
---


"""
import argparse
import json
import logging
from contextlib import redirect_stderr
from pathlib import Path
from typing import Any, Dict

import h5py
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


def compute_properties_with_cat(smile: str, input_file: str, workdir: str) -> None:
    """Compute properties for the given smile and write then down in JSON format.

    Parameters
    ----------
    smile
        String representing the molecular for which the properties are going to be computed

    """
    # Add the smile to the CAT input
    inp = generate_input(smile, input_file)

    with open("cat_output.log", 'a') as f:
        with redirect_stderr(f):
            prep(Settings(inp))

    path_hdf5 = Path(workdir) / "database" / "structures.hdf5"

    results = extract_ligand_properties(path_hdf5)

    with open("results.json", 'w') as handler:
        json.dump(results, handler, indent=4)


def generate_input(smile: str, input_file: str) -> Dict[str, Any]:
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
