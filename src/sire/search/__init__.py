__all__ = [
    "approx_equal",
    "approx_greater",
    "approx_greater_equal",
    "approx_less",
    "approx_less_equal",
    "approx_not_equal",
    "delete_all_tokens",
    "delete_token",
    "get_approx_epsilon",
    "get_min_protein_residues",
    "get_protein_residue_names",
    "get_token",
    "has_token",
    "set_approx_epsilon",
    "set_min_protein_residues",
    "set_protein_residue_names",
    "set_token",
]

# Imported so that it is pythonized
from ..legacy import Search as _Search  # noqa: F401

from .. import use_new_api as _use_new_api

from ..legacy.Search import (
    approx_equal,
    approx_greater,
    approx_greater_equal,
    approx_less,
    approx_less_equal,
    approx_not_equal,
    set_approx_epsilon,
    get_approx_epsilon,
    get_min_protein_residues,
    set_min_protein_residues,
    get_protein_residue_names,
    set_protein_residue_names,
    set_token,
    has_token,
    get_token,
    delete_token,
    delete_all_tokens,
)

_use_new_api()
