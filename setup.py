#!/usr/bin/env python

import os

from setuptools import setup

HERE = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(HERE, 'flamingo', '__version__.py')) as f:
    exec(f.read(), version)

with open('README.rst') as readme_file:
    README = readme_file.read()

setup(
    name='flamingo',
    version=version['__version__'],
    description="Compute and filter molecular properties",
    long_description=README + '\n\n',
    author="Felipe Zapata",
    author_email='f.zapata@esciencecenter.nl',
    url='https://github.com/https://github.com/nlesc-nano/flamingo',
    package_dir={'flamingo': 'flamingo'},
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
        'CAT@git+https://github.com/nlesc-nano/CAT@master',
        'nano-CAT@git+https://github.com/nlesc-nano/nano-CAT@master',
        'data-CAT@git+https://github.com/nlesc-nano/data-CAT@master',
        'mendeleev', 'more_itertools', 'numpy', 'pandas',
        'pyyaml>=5.1.1', 'retry', 'schema'],
    entry_points={
        'console_scripts': [
            'smiles_screener=flamingo.screen:main'
        ]
    },
    data_files=[('citation/flamingo', ['CITATION.cff'])],
    extras_require={
        'test': ['coverage', 'mypy', 'pycodestyle', 'pytest>=3.9', 'pytest-cov',
                 'pytest-mock'],
        'doc': ['sphinx', 'sphinx-autodoc-typehints', 'sphinx_rtd_theme',
                'nbsphinx']
    }
)
