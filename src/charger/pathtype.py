from argparse import ArgumentTypeError
from pathlib import Path
from typing import Union, Callable, Optional


class PathType:
    """Path type of argparse argument.

    Adapted from https://stackoverflow.com/a/33181083

    Examples:
        >>> parser = argparse.ArgumentParser()

        Add an argument that must be an existing file, but can also be specified as a dash ('-') in the command,
        >>> parser.add_argument('existing_file', type=PathType(type='file', accept_dash=True, exists=True))

        Add an argument for a new folder that must NOT exist, but the parent folder exists,
        >>> from pathlib import Path
        >>> NonexistingFolderType = PathType(exists=False, type=lambda p: p.parent.is_dir())
        >>> parser.add_argument('non_existant_folder', type=NonexistingFolderType)
    """

    exists: Optional[bool]
    check_type: Union[None, str, Callable[[Path], bool]]

    def __init__(self, exists = False, type = None):
        """
        exists:
            True: a path that does exist
            False: a path that does not exist, in a valid parent directory
            None: don't care
        type: 'file', 'dir', None, or a function returning True for valid paths
            None: don't care
        """

        assert exists in (True, False, None)
        assert type in ("file", "dir", None) or hasattr(type, "__call__")

        self.check_exist = exists
        self.check_type = type

    def __call__(self, pth: str) -> Path:
        p = Path(pth)
        if self.check_exist is not None:
            # Check if the path must or must not exist
            pth_exists = p.exists()
            if self.check_exist and not pth_exists:
                raise ArgumentTypeError(f"path does not exist: {pth}")
            elif pth_exists:
                raise ArgumentTypeError(f"path already exists: {pth}")

        if self.check_type is None:
            pass
        elif isinstance(self.check_type, str):
            # Check if the path is the asserted type
            if self.check_type == "file" and not p.is_file():
                raise ArgumentTypeError(f"path is not a file: {pth}")
            elif self.check_type == "dir" and not p.is_dir():
                raise ArgumentTypeError(f"path is not a directory: {pth}")
        elif not self.check_type(p):
            # Use the given function to check the path type
            raise ArgumentTypeError(f"path is invalid: {pth}")

        return p
