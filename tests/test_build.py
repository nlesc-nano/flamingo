import os
import shutil
import subprocess
from pathlib import Path
from typing import Container

BUILD = Path('build')
FLAMINGO = Path('flamingo')
FLAMINGO_BUILD = BUILD / 'lib' / 'flamingo'


def _recursive_compare(
    folder1: Path, folder2: Path, ignore: Container[str], allowed_ext: Container[str]
) -> None:
    folder_list1 = os.listdir(folder1)
    folder_list2 = os.listdir(folder2)
    for f1 in folder_list1:
        if f1 in ignore:
            continue
        elif os.path.isdir(folder1 / f1):
            abs_folder1 = folder1 / f1
            abs_folder2 = folder2 / f1
            assert os.path.isdir(abs_folder2)
            _recursive_compare(
                abs_folder1, abs_folder2, ignore=ignore, allowed_ext=allowed_ext
            )

        _, ext = os.path.splitext(f1)
        if ext not in allowed_ext:
            continue
        assert f1 in folder_list2


def test_isfile() -> None:
    """Check that all files ``.py`` and ``.gz`` files in the root directory are also present when the package is build."""  # noqa: E501
    try:
        out = subprocess.run(f'python setup.py build', shell=True, check=True)
        out.check_returncode()

        _recursive_compare(
            FLAMINGO, FLAMINGO_BUILD, ignore={'__pycache__'}, allowed_ext={'.py', '.gz'}
        )
    finally:
        if os.path.isdir(BUILD):
            shutil.rmtree(BUILD)
