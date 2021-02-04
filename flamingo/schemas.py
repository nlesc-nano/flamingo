"""Module with the schemas to validate the user input.

.. autodata:: SCHEMA_SCREEN
    :annotation: : schema.Schema
.. autodata:: SCHEMA_FILTERS
    :annotation: : schema.Schema
.. autodata:: SCHEMA_ORDERING
    :annotation: : schema.Schema

"""
import tempfile
from numbers import Real

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

#: Schema to validate the filters to apply for screening
SCHEMA_FILTERS = Schema({
    # Filter out molecules with more a given anchor
    Optional("single_anchor", default=True): bool,

    # Include or exclude one or more functional group using smiles
    Optional("include_functional_groups"): Schema([str]),
    Optional("exclude_functional_groups"): Schema([str]),

    # Select smiles >, < or = to some value
    Optional("bulkiness"): SCHEMA_BULKINESS,

    Optional("scscore"): SCHEMA_ORDERING

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
