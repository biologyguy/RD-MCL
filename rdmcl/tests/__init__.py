# coding=utf-8
import os
import sys
import pandas as pd
from hashlib import md5
from copy import deepcopy
from buddysuite import SeqBuddy as Sb

SEP = os.sep

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, "%s%s.." % (DIRECTORY_SCRIPT, SEP))
RESOURCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'unit_test_resources') + SEP

cteno_panxs = Sb.SeqBuddy("%s%sCteno_pannexins.fa" % (RESOURCE_PATH, SEP))
ids = [rec.id for rec in cteno_panxs.records]
sim_scores = pd.read_csv("%sCteno_pannexins_sim.scores" % RESOURCE_PATH, "\t", index_col=False, header=None)
sim_scores.columns = ["seq1", "seq2", "score"]


# #################################  -  Helper class  -  ################################## #
class HelperMethods(object):
    def __init__(self):
        self.sep = SEP
        self.resource_path = RESOURCE_PATH
        self._cteno_panxs = cteno_panxs
        self._cteno_ids = ids
        self._cteno_sim_scores = sim_scores

    @staticmethod
    def string2hash(_input):
        return md5(_input.encode("utf-8")).hexdigest()

    def base_cluster_args(self):
        return list(self._cteno_ids), deepcopy(self._cteno_sim_scores)

    def get_data(self, data):
        if data == "cteno_panxs":
            return deepcopy(self._cteno_panxs)
        elif data == "cteno_ids":
            return deepcopy(self._cteno_ids)
        elif data == "cteno_sim_scores":
            return deepcopy(self._cteno_sim_scores)
        else:
            raise AttributeError("Unknown data type: %s" % data)
