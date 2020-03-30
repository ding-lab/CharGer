from csv import DictReader
from pathlib import Path

import pytest

from charger.io import read_csv, read_lines


def test_read_csv(test_root: Path):
    pth = test_root / "test_files" / "example.csv.gz"
    reader = read_csv(pth)
    assert next(reader) == ["Col A", "Col_B"]
    assert next(reader) == ["A", "1.35"]
    assert next(reader) == ["1.2", "Long test"]
    assert reader.gi_yieldfrom.line_num == 3  # type: ignore
    with pytest.raises(StopIteration):
        assert next(reader)


def test_read_csv_asdict(test_root: Path):
    pth = test_root / "test_files" / "example.csv.gz"
    # Disable the type check because we are exploiting the "yield from" internals
    reader = read_csv(pth, as_dict=True)
    assert next(reader) == {"Col A": "A", "Col_B": "1.35"}
    assert reader.gi_yieldfrom.fieldnames == ["Col A", "Col_B"]  # type: ignore
    assert isinstance(reader.gi_yieldfrom, DictReader)
    assert next(reader) == {"Col A": "1.2", "Col_B": "Long test"}
    assert reader.gi_yieldfrom.line_num == 3
    with pytest.raises(StopIteration):
        assert isinstance(reader.gi_yieldfrom, DictReader)
        assert next(reader)


def test_read_lines(test_root: Path):
    pth = test_root / "test_files" / "example.csv.gz"
    lines = read_lines(pth)
    assert len(lines) == 3
    assert lines[0] == '"Col A",Col_B'
    assert lines[1] == "A,1.35"
    assert lines[2] == '1.2,"Long test"'
