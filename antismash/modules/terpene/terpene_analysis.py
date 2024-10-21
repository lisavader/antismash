# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis section of the terpene module
"""

from typing import Dict, Iterable, List
from collections import defaultdict

from antismash.common import fasta, path, subprocessing
from antismash.common.hmmscan_refinement import refine_hmmscan_results, HMMResult
from antismash.common.utils import get_hmm_lengths
from antismash.common.secmet import Protocluster, CDSFeature

from .results import ProtoclusterPrediction, DomainPrediction


def run_terpene_hmmscan(cds_features: Iterable[CDSFeature]) -> Dict[str, List[HMMResult]]:
    """ Runs hmmscan for terpene proteins on the given CDSFeatures

        Arguments:
            cluster: Protocluster on which to run the terpene hmmscan

        Returns:
            a dictionary of key: cds and value: list of HMMResults, for hmmscan results of the cluster
    """
    cluster_fasta = fasta.get_fasta_from_features(cds_features)
    hmm_file = path.get_full_path(__file__, "data", "subprofiles_all.hmm")
    hmm_results = subprocessing.run_hmmscan(hmm_file, cluster_fasta, opts=['--cut_tc'])
    hmm_lengths = get_hmm_lengths(hmm_file)
    return refine_hmmscan_results(hmm_results, hmm_lengths)


def get_cds_predictions(subprofile_hits: Dict[str, List[HMMResult]]) -> Dict[str, List[DomainPrediction]]:
    """ Convert list of HMMResults in CDS mapping to a list of DomainPredictions

        Arguments:
            subprofile_hits: a mapping of CDS name to a list of HMMResults found by hmmscan

        Returns:
            cds_predictions: a mapping of CDS name to a list of DomainPredictions
    """
    cds_predictions : Dict[str, List[DomainPrediction]] = defaultdict(list)

    for cds_name, hmm_results in subprofile_hits.items():
        for hmm_result in hmm_results:
            domain_pred = DomainPrediction(type = "type_not_implemented", subtype = hmm_result.hit_id,
                                        start = hmm_result.query_start, end = hmm_result.query_end)
            cds_predictions[cds_name].append(domain_pred)
    return cds_predictions


def analyse_cluster(cluster: Protocluster) -> ProtoclusterPrediction:
    """ Analyse a terpene cluster

        Arguments:
            cluster: the Protocluster to analyse

        Returns:
            a single ProtoclusterPrediction instance with analysis results
    """
    assert cluster.product == "terpene"

    subprofile_hits = run_terpene_hmmscan(cluster.definition_cdses)
    cds_predictions = get_cds_predictions(subprofile_hits)

    return ProtoclusterPrediction(cds_predictions, cluster.location.start, cluster.location.end)
