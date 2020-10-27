"""
Module to compute molecular properties.

API
---


"""
import argparse
import logging
from contextlib import redirect_stderr
from pathlib import Path
from typing import Any, Dict

import pandas as pd
import yaml
from CAT.base import prep
from scm.plams import Settings

from .log_config import configure_logger
from .cat_interface import PropertyMetadata, extract_dataframe_from_hdf5

logger = logging.getLogger(__name__)


def compute_properties_with_cat(smile: str, input_file: str, workdir: str) -> None:
    """Compute properties for the given smile and write then down in csv format.

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

    metadata  =  PropertyMetadata("bulkiness", 'qd/properties/V_bulk')
    df = extract_dataframe_from_hdf5(path_hdf5, metadata)

    df.to_csv("results.csv", columns=["V_bulk"])


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
