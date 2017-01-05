# coding=utf-8
import os
import sys
from hashlib import md5

SEP = os.sep

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, "%s%s.." % (DIRECTORY_SCRIPT, SEP))
RESOURCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'unit_test_resources') + SEP


# #################################  -  Helper class  -  ################################## #
class HelperMethods(object):
    def __init__(self):
        self.resource_path = RESOURCE_PATH
        self.sep = SEP

    @staticmethod
    def string2hash(_input):
        return md5(_input.encode("utf-8")).hexdigest()
