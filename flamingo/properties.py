"""
Module to compute molecular properties.

API
---


"""
import argparse
import logging
from pathlib import Path

from rdkit import Chem

import pandas as pd

from .features.featurizer import generate_fingerprints
from .log_config import configure_logger
from .models.scscore import SCScorer

logger = logging.getLogger(__name__)


def compute_properties(smile: str, workdir: str) -> None:
    """Compute properties for the given smile.

    Parameters
    ----------
    smile
        String representing the molecular for which the properties are going to be computed

    """
    mol = Chem.MolFromSmiles(smile)
    fingerprint = generate_fingerprints(mol, "morgan", 1024, use_chirality=True).reshape(1, 1024)
    scorer = SCScorer("1024bool")
    prediction = scorer.compute_score(fingerprint)

    # Store results
    path_results = Path(workdir) / "result.csv"
    table = pd.DataFrame.from_dict({"scscore": prediction})
    table.to_csv(path_results)


def main():
    """Parse the command line arguments to screen smiles."""
    parser = argparse.ArgumentParser("compute_properties")
    parser.add_argument('-s', '--smile', required=True, help="smile to compute properties")
    parser.add_argument('-w', '--workdir', required=True, help="Work directory")
    args = parser.parse_args()

    # configure logger
    configure_logger(Path("."), "flamingo")

    compute_properties(args.s, args.w)
