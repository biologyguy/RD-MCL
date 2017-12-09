# coding=utf-8
import os
import sys
import pandas as pd
from hashlib import md5
from copy import deepcopy
from buddysuite import SeqBuddy as Sb
from buddysuite import AlignBuddy as Alb
from io import StringIO

SEP = os.sep

DIRECTORY_SCRIPT = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, "%s%s.." % (DIRECTORY_SCRIPT, SEP))
RESOURCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'unit_test_resources') + SEP

'''
cteno_pannexins.fa contains 134 sequences, 16 of which are removed as paralog cliques, leaving 118

taxon   total_paralogs  minus_cliques
BOL     8
Bab     5
Bch     5
Bfo     10
Bfr     4
Cfu     6
Dgl     9
Edu     9
Hca     8
Hru     5
Hvu     14
Lcr     12
Lla     3
Mle     12
Oma     4
Pba     7
Tin     6
Vpa     7
'''
cteno_panxs = Sb.SeqBuddy("%s%sCteno_pannexins.fa" % (RESOURCE_PATH, SEP))
cteno_panxs_aln = Alb.AlignBuddy("%s%sCteno_pannexins_aln.fa" % (RESOURCE_PATH, SEP))
ids = sorted([rec.id for rec in cteno_panxs.records])
sim_scores = pd.read_csv("%sCteno_pannexins_sim.scores" % RESOURCE_PATH, index_col=False, header=None)
sim_scores.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]


# #################################  -  Helper class  -  ################################## #
class HelperMethods(object):
    def __init__(self):
        self.sep = SEP
        self.resource_path = RESOURCE_PATH
        self._cteno_panxs = cteno_panxs
        self._cteno_panxs_aln = cteno_panxs_aln
        self._cteno_ids = ids
        self._cteno_sim_scores = sim_scores

    @staticmethod
    def string2hash(_input):
        return md5(_input.encode("utf-8")).hexdigest()

    def base_cluster_args(self):
        return list(self._cteno_ids), deepcopy(self._cteno_sim_scores)

    def get_data(self, data):
        if data == "cteno_panxs":
            return Sb.make_copy(self._cteno_panxs)
        elif data == "cteno_panxs_aln":
            return Alb.make_copy(self._cteno_panxs_aln)
        elif data == "cteno_ids":
            return deepcopy(self._cteno_ids)
        elif data == "cteno_sim_scores":
            return deepcopy(self._cteno_sim_scores)
        elif data == "ss2_dfs":
            psi_pred_ss2_dfs = Sb.OrderedDict()
            for rec in cteno_panxs.records:
                path = os.path.join(self.resource_path, "psi_pred", "%s.ss2" % rec.id)
                psi_pred_ss2_dfs[rec.id] = pd.read_csv(path, comment="#", header=None, delim_whitespace=True)
                psi_pred_ss2_dfs[rec.id].columns = ["indx", "aa", "ss", "coil_prob", "helix_prob", "sheet_prob"]
            return psi_pred_ss2_dfs
        elif data == "ss2_paths":
            psi_pred_ss2 = Sb.OrderedDict()
            for rec in cteno_panxs.records:
                psi_pred_ss2[rec.id] = os.path.join(self.resource_path, "psi_pred", "%s.ss2" % rec.id)
            return psi_pred_ss2
        else:
            raise AttributeError("Unknown data type: %s" % data)

    @staticmethod
    def get_db_graph(id_hash, broker):
        graph = broker.query("SELECT (graph) FROM data_table WHERE hash='%s'" % id_hash)
        if graph and graph[0][0]:
            graph = pd.read_csv(StringIO(graph[0][0]), index_col=False, header=None)
            graph.columns = ["seq1", "seq2", "subsmat", "psi", "raw_score", "score"]
        else:
            graph = pd.DataFrame(columns=["seq1", "seq2", "subsmat", "psi", "raw_score", "score"])
        return graph

    def get_sim_scores(self, id_subset):
        df = self._cteno_sim_scores.loc[self._cteno_sim_scores['seq1'].isin(id_subset)]
        df = df.loc[df["seq2"].isin(id_subset)]
        return df

    def get_test_clusters(self, broker, parent_sb, rdmcl):
        # psi_pred_ss2_paths = self.get_data("ss2_paths")
        # graph, alignbuddy = rdmcl.retrieve_all_by_all_scores(parent_sb, psi_pred_ss2_paths, broker)
        graph = self.get_db_graph("b4e3c45e680cdc39b293a778bbc60e4e", broker)
        parent_ids = [rec.id for rec in parent_sb.records]
        parent_cluster = rdmcl.Cluster(parent_ids, graph, collapse=False)

        cluster1 = ['BOL-PanxαB', 'Bab-PanxαA', 'Bch-PanxαA', 'Bfo-PanxαE', 'Bfr-PanxαA', "Oma-PanxαB"]
        # seqs = rdmcl.Sb.pull_recs(rdmcl.Sb.make_copy(parent_sb), "^%s$" % "$|^".join(cluster1))
        # graph, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqs, psi_pred_ss2_paths, broker)
        graph = self.get_db_graph("31777c0c0014ca26bfef56f33ad434bc", broker)
        cluster1 = rdmcl.Cluster(cluster1, graph, parent=parent_cluster)
        cluster1.set_name()

        cluster2 = ['Lla-PanxαA', 'Mle-Panxα11', 'Oma-PanxαD', 'Pba-PanxαB', 'Tin-PanxαF']
        # seqs = rdmcl.Sb.pull_recs(rdmcl.Sb.make_copy(parent_sb), "^%s$" % "$|^".join(cluster2))
        # graph, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqs, psi_pred_ss2_paths, broker)
        graph = self.get_db_graph("441c3610506fda8a6820da5f67fdc470", broker)
        cluster2 = rdmcl.Cluster(cluster2, graph, parent=parent_cluster)
        cluster2.set_name()

        cluster3 = ['Vpa-PanxαD']  # This should be placed in cluster 2
        # seqs = rdmcl.Sb.pull_recs(rdmcl.Sb.make_copy(parent_sb), "^%s$" % "$|^".join(cluster3))
        # graph, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqs, psi_pred_ss2_paths, broker)
        graph = self.get_db_graph("ded0bc087974589c24945bb36197d36f", broker)
        cluster3 = rdmcl.Cluster(cluster3, graph, parent=parent_cluster)
        cluster3.set_name()

        cluster4 = ['Hca-PanxαA', 'Lcr-PanxαG']  # This should be placed in cluster 1
        # seqs = rdmcl.Sb.pull_recs(rdmcl.Sb.make_copy(parent_sb), "^%s$" % "$|^".join(cluster4))
        # graph, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqs, psi_pred_ss2_paths, broker)
        graph = self.get_db_graph("035b3933770942c7e32e27c06e619825", broker)
        cluster4 = rdmcl.Cluster(cluster4, graph, parent=parent_cluster)
        cluster4.set_name()

        cluster5 = ['Hvu-PanxβA']  # This should not be placed in a cluster
        # seqs = rdmcl.Sb.pull_recs(rdmcl.Sb.make_copy(parent_sb), "^%s$" % "$|^".join(cluster5))
        # graph, alignbuddy = rdmcl.retrieve_all_by_all_scores(seqs, psi_pred_ss2_paths, broker)
        graph = self.get_db_graph("c62b6378d326c2479296c98f0f620d0f", broker)
        cluster5 = rdmcl.Cluster(cluster5, graph, parent=parent_cluster)
        cluster5.set_name()
        clusters = [cluster1, cluster2, cluster3, cluster4, cluster5]
        return clusters
