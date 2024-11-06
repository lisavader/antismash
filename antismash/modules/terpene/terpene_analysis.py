# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis section of the terpene module
"""
import json

from typing import Dict, Iterable, List
from collections import defaultdict

from antismash.common import fasta, path
from antismash.common.hmmscan_refinement import HMMResult, HSP
from antismash.common.secmet import Protocluster, CDSFeature

from .results import ProtoclusterPrediction, DomainPrediction, TerpeneHMM
from .hmmscan import run_hmmscan


_HMM_PROPERTIES_CACHE: list[TerpeneHMM] = []


def get_hmm_properties() -> list[TerpeneHMM]:
    """ Generates the HMM properties from hmm_properties.json and compound_groups.json
        Only does the processing once per python invocation, future runs access
        existing properties
    """
    # if already called once, then just reuse the cached results
    if not _HMM_PROPERTIES_CACHE:

        loaded_json = []
        for file_name in ["hmm_properties", "compound_groups"]:
            file_path = path.get_full_path(__file__, "data", file_name + ".json")
            with open(file_path) as handle:
                try:
                    loaded_json.append(json.load(handle))
                except ValueError as error:
                    raise ValueError(f"{file_path!r} is not a valid json file: {error}")
        properties_json = loaded_json[0]
        compounds_json = loaded_json[1]

        terpene_hmms = []
        for hmm_name, hmm_data in properties_json.items():
            terpene_hmms.append(TerpeneHMM.from_json(hmm_name, hmm_data, compounds_json))

        _HMM_PROPERTIES_CACHE.extend(terpene_hmms)

    return _HMM_PROPERTIES_CACHE


def run_terpene_hmmscan(cds_features: Iterable[CDSFeature],
                        hmm_by_name: dict[str, TerpeneHMM]) -> Dict[str, set[HMMResult]]:
    """ Runs hmmscan for terpene proteins on the given CDSFeatures

        Arguments:
            cluster: Protocluster on which to run the terpene hmmscan

        Returns:
            a dictionary of key: cds and value: list of HMMResults, for hmmscan results of the cluster
    """
    cluster_fasta = fasta.get_fasta_from_features(cds_features)
    hmm_file = path.get_full_path(__file__, "data", "all_profiles.hmm")
    raw_results = run_hmmscan(hmm_file, cluster_fasta)

    results_by_id: Dict[str, set[HMMResult]] = defaultdict(set)
    for result in raw_results:
        # Store result if it is above cut-off
        for hsp in result.hsps:
            if hsp.hit_id in hmm_by_name:
                terpene_hmm = hmm_by_name[hsp.hit_id]
            else:
                raise ValueError(f"Failed to find signature for ID {hsp.hit_id}")
            if hsp.bitscore > terpene_hmm.cutoff:
                 results_by_id[hsp.query_id].add(HMMResult(hsp.hit_id, hsp.query_start,
                                                      hsp.query_end, hsp.evalue,
                                                      hsp.bitscore))
    return results_by_id


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

    hmm_properties_by_name = {terpene_hmm.name: terpene_hmm for terpene_hmm in get_hmm_properties()}
    hmm_hits = run_terpene_hmmscan(cluster.definition_cdses, hmm_properties_by_name)
    cds_predictions = get_cds_predictions(hmm_hits)

    return ProtoclusterPrediction(cds_predictions, cluster.location.start, cluster.location.end)
