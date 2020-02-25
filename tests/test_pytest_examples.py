from textwrap import dedent

from loguru import logger


def test_charger_version():
    # make sure the version is the same on PyPI and in __init__.py
    import pkg_resources
    from charger import __version__

    assert pkg_resources.get_distribution("charger").version == __version__


def test_log_capture(caplog):
    logger.info("Log something")
    assert len(caplog.records) == 1
    assert "Log something" in caplog.text


# @pytest.mark.xfail(reason="Demo of long string comparison")
def test_long_string():
    long_string_a = dedent(
        """
        Basic research is performed without thought of practical ends.
        It results in general knowledge and an understanding of nature and its laws.

        A worker in basic scientific research is motivated by a driving curiosity about the unknown.
        """
    )
    long_string_b = dedent(
        """
        Basic research is performed without ends of practical thoughts.
        It results in general knowledge and an understanding of nature and its laws.

        A worker in basic scientific research is motivated by a driving unknown about curiosity.
        """
    )
    assert long_string_a != long_string_b
