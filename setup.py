#!/usr/bin/env python

import importlib
import os

from setuptools import find_packages, setup

HERE = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(HERE, 'flamingo', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.rst') as readme_file:
    README = readme_file.read()

try:
    importlib.import_module("rdkit")
    importlib.import_module("h5py")
except ModuleNotFoundError:
    exc = ModuleNotFoundError(
        """'flamingo' requires the 'rdkit' package: https://anaconda.org/conda-forge/rdkit
and the h5py package: https://anaconda.org/conda-forge/h5py"""
    )
    exc.__cause__ = None
    raise exc


setup(
    name='flamingo',
    version=version['__version__'],
    description="Compute and filter molecular properties",
    long_description=README + '\n\n',
    author="Felipe Zapata",
    author_email='f.zapata@esciencecenter.nl',
    url='https://github.com/https://github.com/nlesc-nano/flamingo',
    packages=find_packages(),
    include_package_data=True,
    license="Apache Software License 2.0",
    zip_safe=False,
    keywords='flamingo',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    install_requires=[
        'CAT@git+https://github.com/nlesc-nano/CAT',
        'nano-CAT@git+https://github.com/nlesc-nano/nano-CAT@master',
        'data-CAT@git+https://github.com/nlesc-nano/data-CAT',
        'mendeleev', 'more_itertools', 'numpy', 'pandas',
        'pyyaml>=5.1.1', 'schema', 'typing_extensions'],
    entry_points={
        'console_scripts': [
            'smiles_screener=flamingo.screen:main',
            'compute_properties=flamingo.properties:main'
        ]
    },
    package_data={
        'flamingo': ['data/scscore/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz',
        'data/scscore/full_reaxys_model_2048bool/model.ckpt-10654.as_numpy.json.gz']
    },
    data_files=[('citation/flamingo', ['CITATION.cff'])],
    extras_require={
        'test': ['coverage', 'mypy', 'pycodestyle', 'pytest>=3.9', 'pytest-cov',
                 'pytest-mock'],
        'doc': ['sphinx', 'sphinx-autodoc-typehints', 'sphinx_rtd_theme',
                'nbsphinx']
    }
)
