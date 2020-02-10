import argparse

import pytest

from charger.argtype import ModuleScoreOverrideType, PathType


@pytest.fixture
def example_folder(tmp_path):
    """A example folder fixture to test file existence.

    It creates one existing folder and file under the ``tmp_path`` directory. It returns a dict of existent and
    non-existent paths.
    """
    folder = tmp_path / "existing_folder"
    folder.mkdir()
    file = tmp_path / "existing_file.txt"
    with file.open("w") as f:
        f.write("This file exists.\n")

    return {
        "exist": {"folder": folder, "file": file},
        "not_exist": {"file": tmp_path / "new_path"},
    }


def test_pathtype_must_exist(example_folder):
    MustExistPathType = PathType(exists=True)
    for pthtype, pth in example_folder["exist"].items():
        assert MustExistPathType(str(pth)) == pth

    for pthtype, pth in example_folder["not_exist"].items():
        with pytest.raises(argparse.ArgumentTypeError) as excinfo:
            MustExistPathType(str(pth))
        assert "Path does not exist" in excinfo.value.args[0]
        assert str(pth) in excinfo.value.args[0]


def test_pathtype_must_not_exist(example_folder):
    MustNotExistPathType = PathType(exists=False)
    for pthtype, pth in example_folder["not_exist"].items():
        assert MustNotExistPathType(str(pth)) == pth

    for pthtype, pth in example_folder["exist"].items():
        with pytest.raises(argparse.ArgumentTypeError) as excinfo:
            MustNotExistPathType(str(pth))
        assert excinfo.value.args[0].startswith("Path already exists")
        assert str(pth) in excinfo.value.args[0]


def test_pathtype_type_requirement(example_folder):
    FileOnlyPathType = PathType(type="file", exists=True)
    DirOnlyPathType = PathType(type="dir", exists=True)
    exist_paths = example_folder["exist"]

    assert FileOnlyPathType(str(exist_paths["file"])) == exist_paths["file"]

    with pytest.raises(argparse.ArgumentTypeError) as excinfo:
        FileOnlyPathType(str(exist_paths["folder"]))
    assert excinfo.value.args[0].startswith("Path is not a file")

    assert DirOnlyPathType(str(exist_paths["folder"])) == exist_paths["folder"]

    with pytest.raises(argparse.ArgumentTypeError) as excinfo:
        DirOnlyPathType(str(exist_paths["file"]))
    assert excinfo.value.args[0].startswith("Path is not a dir")


def test_pathtype_invalid_type(example_folder):
    PipeOnlyPathType = PathType(type="pipe")
    with pytest.raises(ValueError):
        PipeOnlyPathType(str(example_folder["exist"]["folder"]))
    with pytest.raises(ValueError):
        PipeOnlyPathType(str(example_folder["exist"]["file"]))


def test_pathtype_custom_type(example_folder):
    TxtPathType = PathType(type=lambda p: p.suffix == ".txt")
    exist_txt_pth = example_folder["exist"]["file"]
    notexist_txt_pth = example_folder["exist"]["folder"] / "new.txt"

    assert TxtPathType(str(exist_txt_pth)) == exist_txt_pth
    assert TxtPathType(str(notexist_txt_pth)) == notexist_txt_pth

    with pytest.raises(argparse.ArgumentTypeError) as excinfo:
        TxtPathType(example_folder["not_exist"]["file"])
    assert excinfo.value.args[0].startswith("Path is invalid")


def test_pathtype_custom_type_must_exist(example_folder):
    TxtPathType = PathType(exists=True, type=lambda p: p.suffix == ".txt")
    exist_txt_pth = example_folder["exist"]["file"]
    notexist_txt_pth = example_folder["exist"]["folder"] / "new.txt"

    assert TxtPathType(str(exist_txt_pth)) == exist_txt_pth
    with pytest.raises(argparse.ArgumentTypeError):
        TxtPathType(str(notexist_txt_pth))


@pytest.fixture
def override_score():
    """A fixture of default module scores."""
    return ModuleScoreOverrideType(defaults={"A": 0, "B": 10})


def test_modulescoreoverridetype(override_score):
    scores = override_score("A=3 B=10")
    assert scores["A"] == 3
    assert scores["B"] == 10


def test_modulescoreoverridetype_invalid_module(override_score):
    with pytest.raises(argparse.ArgumentTypeError, match="Module does not exist: PPAP"):
        override_score("A=10 PPAP=999")


def test_modulescoreoverridetype_invalid_score(override_score):
    with pytest.raises(
        argparse.ArgumentTypeError, match="New score of module A is not an integer:"
    ):
        override_score("A=0.99")
