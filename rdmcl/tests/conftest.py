#!/usr/bin/env python3
# coding=utf-8
""" Fixtures for py.test  """
import pytest

from .__init__ import RESOURCE_PATH
from . import __init__ as init


# #################################  -  Helper functions  -  ################################## #
@pytest.fixture(scope="session")
def hf():
    """
    Collection of helper methods
    """
    return init.HelperMethods()
