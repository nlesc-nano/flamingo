"""Functions use for testing."""

import os
from pathlib import Path

import flamingo

__all__ = ["PATH_FLAMINGO", "PATH_TEST"]

# Environment data
PATH_FLAMINGO = Path(os.path.dirname(flamingo.__file__))
ROOT = PATH_FLAMINGO.parent

PATH_TEST = ROOT / "tests" / "files"
