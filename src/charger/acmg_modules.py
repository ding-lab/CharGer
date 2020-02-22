# To avoid circular import, we only import the classes required by type
# checking when TYPE_CHECKING is TRUE. And use PEP 563 to avoid forward
# reference (__future__ imports on python 3.7+)
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .classifier import CharGerResult


def run_pvs1_module(result: CharGerResult) -> None:
    pass
