from argparse import ArgumentTypeError as err
from pathlib import Path
import os.path
from typing import Union, Callable


class PathType:
    """Path type of argparse argument.

    Adapted from https://stackoverflow.com/a/33181083
    And argparse.FileType https://github.com/python/cpython/blob/cb1c0746f277052e45a60d6c436a765e34722821/Lib/argparse.py#L1227

    Examples:
        >>> parser = argparse.ArgumentParser()

        Add an argument that must be an existing file, but can also be specified as a dash ('-') in the command.
        >>> parser.add_argument('existing_file', type=PathType(type='file', dash_ok=True, exists=True))

        Add an argument for a folder that must NOT exist.
        >>> parser.add_argument('non_existant_folder', type=PathType(type='dir', exists=False))

        Add an argument for EITHER a folder or a file, but can be dash_ok, and don't check if it exists.
        >>> existing_dir_or_file_type = PathType(type=('dir','file').__contains__, dash_ok=False, exists=None)
        >>> parser.add_argument('exiting_file_or_folder', type=existing_dir_or_file_type)
    """

    _type: Union[None, str, Callable[[str], bool]]

    def __init__(
        self,
        exists: bool = False,
        type: Union[None, str, Callable[[str], bool]] = "file",
        dash_ok: bool = False,
    ):
        """
        exists:
            True: a path that does exist
            False: a path that does not exist, in a valid parent directory
            None: don't care
        type: 'file', 'dir', 'symlink', None, or a function returning True for valid paths
            None: don't care
        dash_ok: whether to allow "-" as stdin/stdout
        """

        assert exists in (True, False, None)
        assert type in ("file", "dir", "symlink", None) or hasattr(type, "__call__")

        self._exists = exists
        self._type = type
        self._dash_ok = dash_ok

    def __call__(self, string: str) -> str:
        if string == "-":
            # the special argument "-" means sys.std{in,out}
            if self._type == "dir":
                raise err("standard input/output (-) not allowed as directory path")
            elif self._type == "symlink":
                raise err("standard input/output (-) not allowed as symlink path")
            elif not self._dash_ok:
                raise err("standard input/output (-) not allowed")
        else:
            e = os.path.exists(string)
            if self._exists:
                if not e:
                    raise err("path does not exist: '%s'" % string)

                if self._type is None:
                    pass
                elif self._type == "file":
                    if not os.path.isfile(string):
                        raise err("path is not a file: '%s'" % string)
                elif self._type == "symlink":
                    if not os.path.islink(string):
                        raise err("path is not a symlink: '%s'" % string)
                elif self._type == "dir":
                    if not os.path.isdir(string):
                        raise err("path is not a directory: '%s'" % string)
                elif not self._type(string):
                    raise err("path is not valid: '%s'" % string)
            else:
                if not self._exists and e:
                    raise err("path exists: '%s'" % string)

                p = os.path.dirname(os.path.normpath(string)) or "."
                if not os.path.isdir(p):
                    raise err("parent path is not a directory: '%s'" % p)
                elif not os.path.exists(p):
                    raise err("parent directory does not exist: '%s'" % p)

        return string
