from csv import DictReader
from pathlib import Path
from typing import Any

import pytest

from charger.io import read_csv


def test_read_csv(test_root: Path):
    pth = test_root / "test_files" / "example.csv"
    reader: Any = read_csv(pth)
    assert next(reader) == ["Col A", "Col_B"]
    assert next(reader) == ["A", "1.35"]
    assert next(reader) == ["1.2", "Long test"]
    assert reader.gi_yieldfrom.line_num == 3
    with pytest.raises(StopIteration):
        assert next(reader)


def test_read_csv_asdict(test_root: Path):
    pth = test_root / "test_files" / "example.csv"
    # Disable the type check because we are exploiting the "yield from" internals
    reader: Any = read_csv(pth, as_dict=True)
    assert next(reader) == {"Col A": "A", "Col_B": "1.35"}
    assert reader.gi_yieldfrom.fieldnames == ["Col A", "Col_B"]
    assert isinstance(reader.gi_yieldfrom, DictReader)
    assert next(reader) == {"Col A": "1.2", "Col_B": "Long test"}
    assert reader.gi_yieldfrom.line_num == 3
    with pytest.raises(StopIteration):
        assert isinstance(reader.gi_yieldfrom, DictReader)
        assert next(reader)
