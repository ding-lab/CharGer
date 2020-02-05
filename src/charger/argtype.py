from argparse import ArgumentTypeError
from pathlib import Path
from typing import Union, Callable, Optional


class PathType:
    """Parse file path arguments as a pathlib.Path object.

    Type factory for ArgumentParser.add_argument() that validates the file path arguments.
    Adapted from https://stackoverflow.com/a/33181083

    Args:
        exists: The path must exist when True, and the path must not exist when False.
            Skip the check when None.
        type: ``file`` or ``dir`` will be checked by Path.is_file() and Path.is_dir(), respectively.
            It also accepts a function that takes Path object and returns True for valid paths.
            Skip the check when None.

    Examples:
        >>> parser = argparse.ArgumentParser()

        Add an argument that must be an existing file, but can also be specified as a dash ('-') in the command,
        >>> parser.add_argument('existing_file', type=PathType(exists=True, type='file'))

        Add an argument for a new path that must NOT exist, but the parent folder exists,
        >>> from pathlib import Path
        >>> CreatablePathType = PathType(exists=False, type=lambda p: p.parent.is_dir())
        >>> parser.add_argument('non_existing_folder', type=CreatablePathType)
    """

    def __init__(
        self,
        exists: Optional[bool] = None,
        type: Union[None, str, Callable[[Path], bool]] = None,
    ):
        self._exists = exists
        self._type = type

    def __call__(self, pth: str) -> Path:
        p = Path(pth)
        if self._exists is not None:
            # Check if the path must or must not exist
            pth_exists = p.exists()
            if self._exists:
                if not pth_exists:
                    raise ArgumentTypeError(f"Path does not exist: {pth}")
            elif pth_exists:
                raise ArgumentTypeError(f"Path already exists: {pth}")

        if self._type is None:
            pass
        elif isinstance(self._type, str):
            # Check if the path is the asserted type
            if self._type == "file":
                if not p.is_file():
                    raise ArgumentTypeError(f"Path is not a file: {pth}")
            elif self._type == "dir":
                if not p.is_dir():
                    raise ArgumentTypeError(f"Path is not a directory: {pth}")
            else:
                raise ValueError(f"Unknown type {self._type!r}")
        elif not self._type(p):
            # Use the given function to check the path type
            raise ArgumentTypeError(f"Path is invalid: {pth}")

        return p
