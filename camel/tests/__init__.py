import logging
import os
import unittest
from collections.abc import Callable

from camel.app.config import config


# noinspection PyPep8Naming
def longRunningTest() -> Callable:
    """
    Custom decorator for long-running tests, used to limit the duration when running all tests.
    :return: Decorator
    """
    is_skipped = (os.environ.get('CAMEL_SKIP_LONG_TESTS')) == '1'
    if is_skipped:
        return unittest.skip("Skipping long test")
    return lambda func: func


# noinspection PyPep8Naming
def resourceIntensiveTest(reason: str | None = None) -> Callable:
    """
    Custom decorator for resource intensive tests, used to limit resource allocation for RAM / CPU intensive tests.
    :param reason: (Optional) clarification for resource allocation
    :return: Decorator
    """
    is_skipped = (os.environ.get('CAMEL_SKIP_RESOURCE_INTENSIVE_TESTS')) == '1'
    if is_skipped:
        return unittest.skip(f"Resource intensive test ({reason})")
    return lambda func: func


def get_ubuntu_release_codename() -> str | None:
    """
    Retrieves the codename of the Ubuntu release version.
    :return: Ubuntu release version (if available)
    """
    try:
        with open('/etc/os-release') as handle:
            for line in handle.readlines():
                if line.startswith("VERSION_CODENAME="):
                    return line.strip().split('=')[1].strip('"')
    except FileNotFoundError:
        logging.warning(
            "Cannot determine release code name, '/etc/os-release' not found"
        )
    return None


# noinspection PyPep8Naming
def minOSVersion(codename_min: str) -> Callable:
    """
    Custom decorator to skip tests on older OS versions.
    Notes:
    - this method currently only supports the Ubuntu OS
    - versions are compared using string comparisons (considering the alphabetical order of version codenames)
    :param codename_min: Release version code name (e.g. 'focal' or 'jammy')
    :return: Decorator
    """
    codename_os = get_ubuntu_release_codename()
    if codename_os < codename_min:
        return unittest.skip(
            f"Skipping test because of OS version (min={codename_min}, current={codename_os})"
        )
    return lambda func: func


def requires_dependency_service(*services: str):
    """
    Decorator to skip a test when the configured dependency service is not in the given services.
    :param services: One or more accepted dependency services ('lmod', 'pixi')
    :return: Test decorator
    """
    return unittest.skipUnless(
        config.dependency_service in services,
        f"Requires dependency service: {', '.join(services)} (current: {config.dependency_service})",
    )
