"""Logger configuration."""

import logging
import sys
import importlib
from pathlib import Path

__all__ = ["configure_logger"]

logger = logging.getLogger(__name__)


def configure_logger(workdir: Path, package_name: str) -> None:
    """Set the logging infrasctucture."""
    pkg = sys.modules.get(package_name)
    if pkg is None:
        pkg = importlib.import_module(package_name)

    file_log = workdir / f'{package_name}_output.log'
    logging.basicConfig(filename=file_log, level=logging.INFO,
                        format='%(asctime)s  %(message)s',
                        datefmt='[%I:%M:%S]')
    handler = logging.StreamHandler()
    handler.terminator = ""

    version = getattr(pkg, "__version__", "UNKNOWN")
    path = Path(pkg.__file__).parent

    logger.info(f"Using {package_name} version: {version}\n")
    logger.info(f"{package_name} path is: {path}\n")
    logger.info(f"Working directory is: {workdir}")


class LoggerWriter:
    """Modify the default behaviour of the logger."""

    def __init__(self, level):
        self.level = level

    def write(self, message):
        # if statement reduces the amount of newlines that are
        # printed to the logger
        if message != '\n':
            self.level(message)

    def flush(self):
        # create a flush method so things can be flushed when
        self.level(sys.stderr)
