"""Functions use for testing."""

from pathlib import Path

import pkg_resources as pkg

__all__ = ["PATH_FLAMINGO", "PATH_TEST"]

# Environment data
PATH_FLAMINGO = Path(pkg.resource_filename('flamingo', ''))
ROOT = PATH_FLAMINGO.parent

PATH_TEST = ROOT / "tests" / "files"
