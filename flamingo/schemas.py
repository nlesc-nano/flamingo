"""Module with the schemas to validate the user input.

.. autodata:: SCHEMA_SCREEN
    :annotation: : schema.Schema
.. autodata:: SCHEMA_FILTERS
    :annotation: : schema.Schema
.. autodata:: SCHEMA_ORDERING
    :annotation: : schema.Schema

"""
import tempfile
from numbers import Integral, Real

import yaml
from schema import Optional, Or, Schema, SchemaError

from .utils import Options

#: Schema to validate the ordering keywords
SCHEMA_ORDERING = Or(
    Schema({"greater_than": Real}),
    Schema({"lower_than": Real}))

#: Schema to validate the bulkiness parameters Check cat documentation
#: https://cat.readthedocs.io/en/latest/4_optional.html#optional.qd.bulkiness
SCHEMA_BULKINESS = Schema({
    "lower_than": Real,
    Optional("d", default="auto"): Or(str, Real, None),
    Optional("h_lim", default=10): Or(Real, None)})


#: Schema to include/exclude multiple functional groups
SCHEMA_FUNCTIONAL_GROUPS = Schema({
    #Functional Groups to include/exclude
    "groups": Schema([str]),
    Optional("maximum", default=1): Integral
})

#: Schema to compute the properties described at:
#: https://www.rdkit.org/docs/source/rdkit.Chem.QED.html
SCHEMA_DRUG_LIKENESS = Schema({
    Optional("MW", default=None): Or(SCHEMA_ORDERING, None),
    Optional("ALOGP", default=None): Or(SCHEMA_ORDERING, None),
    Optional("HBA", default=None): Or(SCHEMA_ORDERING, None),
    Optional("HBD", default=None): Or(SCHEMA_ORDERING, None),
    Optional("PSA", default=None): Or(SCHEMA_ORDERING, None),
    Optional("ROTB", default=None): Or(SCHEMA_ORDERING, None),
    Optional("AROM", default=None): Or(SCHEMA_ORDERING, None),
    Optional("ALERTS", default=None): Or(SCHEMA_ORDERING, None)
})

SCHEMA_COSMO_RS = Schema({
    "solvents": dict,
})

#: Schema to validate the filters to apply for screening
SCHEMA_FILTERS = Schema({
    # Include or exclude one or more functional group using smiles
    Optional("include_functional_groups"): SCHEMA_FUNCTIONAL_GROUPS,
    Optional("exclude_functional_groups"): SCHEMA_FUNCTIONAL_GROUPS,

    # Select smiles >, < or = to some value
    Optional("bulkiness"): SCHEMA_BULKINESS,

    Optional("scscore"): SCHEMA_ORDERING,

    Optional("drug_likeness", default=None): Or(SCHEMA_DRUG_LIKENESS, None),

    Optional("cosmo_rs"): SCHEMA_COSMO_RS
})

#: Schema to validate the input for screening
SCHEMA_SCREEN = Schema({
    # Load the dataset from a file
    "smiles_file": str,

    # Constrains to filter
    "filters": SCHEMA_FILTERS,

    # Functional group used as anchor
    Optional("anchor", default="O(C=O)[H]"): str,

    # path to the molecular coordinates of the Core to attach the ligands
    Optional("core"): str,

    # path to the workdir
    Optional("workdir", default=tempfile.mkdtemp(prefix="flamingo_workdir_")): str,

    # File to print the final candidates
    Optional("output_path", default="results"): str,

    # Number of molecules per batch
    Optional("batch_size", default=1000): int,

    # Run filters in Parallel
    Optional("parallel", default=True): bool,

    # Clean up the smiles
    Optional("sanitize_smiles", default=False): bool

})

DICT_ACTIONS = {"screen": SCHEMA_SCREEN}


def validate_input(file_input: str, action: str) -> Options:
    """Check the input validation against an schema."""
    action_schema = DICT_ACTIONS[action]
    with open(file_input, 'r') as handler:
        dict_input = yaml.load(handler.read(), Loader=yaml.FullLoader)
    try:
        data = action_schema.validate(dict_input)
        return Options(data)
    except SchemaError as err:
        msg = f"There was an error in the input yaml provided:\n{err}"
        print(msg)
        raise
