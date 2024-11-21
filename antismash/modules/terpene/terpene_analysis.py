# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis section of the terpene module
"""
import json
import math

from typing import Dict, Iterable, List
from collections import defaultdict

from antismash.common import fasta, path
from antismash.common.hmmscan_refinement import HMMResult, refine_hmmscan_results, QueryResult
from antismash.common.secmet import Protocluster, CDSFeature
from antismash.common.subprocessing.hmmscan import run_hmmscan, OutputType

from .results import ProtoclusterPrediction, DomainPrediction, TerpeneHMM, CompoundGroup


_HMM_PROPERTIES_CACHE: dict[str, TerpeneHMM] = {}

_TYPE_MAPPINGS = {"phytoene_synt": "PT_phytoene_like"}

def load_json(json_file):
    file_path = path.get_full_path(__file__, "data", json_file + ".json")
    with open(file_path) as handle:
        try:
            return json.load(handle)
        except ValueError as error:
            raise ValueError(f"{file_path!r} is not a valid json file: {error}")


def get_hmm_properties() -> dict[str, TerpeneHMM]:
    """ Generates the HMM properties from hmm_properties.json and compound_groups.json
        Only does the processing once per python invocation, future runs access
        existing properties
    """
    # if already called once, then just reuse the cached results
    if not _HMM_PROPERTIES_CACHE:
        properties_json = load_json("hmm_properties")
        compounds_json = load_json("compound_groups")

        compound_groups = {}
        for group_data in compounds_json:
            compound_group = CompoundGroup.from_json(group_data)
            compound_groups[compound_group.name] = compound_group

        terpene_hmms = {}
        for hmm_data in properties_json:
            terpene_hmm = TerpeneHMM.from_json(hmm_data, terpene_hmms, compound_groups)
            terpene_hmms[terpene_hmm.name] = terpene_hmm

        # Retrieve the main type from parents
        def map_name(name):
            try:
                name_out = _TYPE_MAPPINGS[name]
            except KeyError:
                name_out = name
            return name_out

        for hmm_name, terpene_hmm in terpene_hmms.items():
            main_type = map_name(hmm_name)
            parent_hmm = terpene_hmm
            while parent_hmm.parents:
                parent_name = parent_hmm.parents[0]
                assert parent_name != hmm_name
                parent_hmm = terpene_hmms[parent_name]
                main_type = map_name(parent_name)
            terpene_hmm.add_main_type(main_type)

        _HMM_PROPERTIES_CACHE.update(terpene_hmms)

    return _HMM_PROPERTIES_CACHE


def run_terpene_hmmscan(cds_features: Iterable[CDSFeature]) -> List[QueryResult]:
    """ Runs hmmscan for terpene proteins on the given CDSFeatures

        Arguments:
            cluster: Protocluster on which to run the terpene hmmscan

        Returns:
            hmmscan_results: a list of QueryResult objects from Bio's SearchIO
    """
    cluster_fasta = fasta.get_fasta_from_features(cds_features)
    hmm_file = path.get_full_path(__file__, "data", "all_profiles.hmm")
    hmmscan_results = run_hmmscan(target_hmmfile = hmm_file, query_sequence = cluster_fasta,
                                  output_type = OutputType.TABULAR)

    return hmmscan_results


def filter_by_score(hmm_results: Dict[str, List[HMMResult]],
                    hmm_properties: dict[str, TerpeneHMM]) -> Dict[str, List[HMMResult]]:
    """ Removes hmm results with bitscores that are below the specified cutoff for the hmm
    """
    results_by_id: Dict[str, list[HMMResult]] = defaultdict(list)
    for cds_name, results in hmm_results.items():
        # Store result if it is above the cut-off
        for result in results:
            if result.hit_id in hmm_properties:
                terpene_hmm = hmm_properties[result.hit_id]
            else:
                raise ValueError(f"Failed to find signature for ID {result.hit_id}")
            if result.bitscore > terpene_hmm.cutoff:
                results_by_id[cds_name].append(result)
    return results_by_id


def remove_overlapping(results: List[HMMResult]) -> List[HMMResult]:
    """ Filters overlapping hmm results, keeping the hit with the lowest evalue.
        If the evalues are too similar, both hits are kept.
    """
    non_overlapping = []
    if results:
        non_overlapping = [results[0]]
        for result in results[1:]:
            previous = non_overlapping[-1]
            if result.overlaps_with(previous):
                # Replace the previous result if the current result scores significantly better
                if math.log10(result.evalue) < 1.5 * math.log10(previous.evalue):
                    non_overlapping[-1] = result
                # Keep the previous result if the previous result scores significantly better
                elif math.log10(previous.evalue) < 1.5 * math.log10(result.evalue):
                    pass
                # Keep both results if there is not substantial difference in evalues
                else:
                    non_overlapping.append(result)
            else:
                non_overlapping.append(result)
    return non_overlapping


def filter_overlaps(hmm_results: Dict[str, List[HMMResult]],
                   hmm_properties: dict[str, TerpeneHMM]) -> Dict[str, List[HMMResult]]:
    """ Removes overlaps within main profile hits and within subtype hits
    """
    results_by_id: Dict[str, list[HMMResult]] = defaultdict(list)
    for cds_name, results in hmm_results.items():
        subtypes: list[HMMResult] = []
        main_types: list[HMMResult] = []
        for result in results:
            if hmm_properties[result.hit_id].is_subtype():
                subtypes.append(result)
            else:
                main_types.append(result)
        results_by_id[cds_name].extend(remove_overlapping(subtypes) +
                                       remove_overlapping(main_types))
    return results_by_id


def get_cds_predictions(hmm_results: Dict[str, List[HMMResult]],
                        hmm_properties: dict[str, TerpeneHMM]) -> Dict[str, List[DomainPrediction]]:
    """ Convert list of HMMResults in CDS mapping to a list of DomainPredictions

        Arguments:
            hmm_results: a mapping of CDS name to a list of HMMResults

        Returns:
            cds_predictions: a mapping of CDS name to a list of DomainPredictions
    """
    cds_predictions : Dict[str, List[DomainPrediction]] = defaultdict(list)

    for cds_name, hmm_results in hmm_results.items():
        # Extract groups of overlapping hmm results
        groups = [[hmm_results[0]]]
        for result in hmm_results[1:]:
            overlapping = False
            for member in groups[-1]:
                if result.overlaps_with(member):
                    overlapping = True
            if overlapping:
                groups[-1].append(result)
            else:
                groups.append([result])

        # Add a domain prediction for each group
        for group in groups:
            start_locations = []
            end_locations = []
            main_types = set()
            subtypes = set()
            all_predictions = []
            merged_predictions = set()
            for hmm_result in group:
                start_locations.append(hmm_result.query_start)
                end_locations.append(hmm_result.query_end)
                terpene_hmm = hmm_properties[hmm_result.hit_id]
                main_types.add(terpene_hmm.main_type)
                if terpene_hmm.is_subtype():
                    subtypes.add(terpene_hmm.name)
                for pred in terpene_hmm.predictions:
                    print(pred)
                    all_predictions.append(pred)
            if all_predictions:
                merged = all_predictions[0]
                for other in all_predictions[1:]:
                    if merged.substrates == other.substrates:
                        merged = merged.merge(other)
                    else:
                        merged = other
                    merged_predictions.add(merged)
            assert len(main_types) == 1, "Overlapping hits cannot belong to different main types"
            domain_pred = DomainPrediction(type = main_types.pop(), subtypes = tuple(subtypes),
                                           start = min(start_locations), end = max(end_locations),
                                           predictions = tuple(merged_predictions))
            print(domain_pred)
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

    hmm_properties = get_hmm_properties()
    hmm_lengths = {hmm_name: hmm_obj.length for hmm_name, hmm_obj in hmm_properties.items()}

    hmmscan_results = run_terpene_hmmscan(cluster.cds_children)
    refined_results = refine_hmmscan_results(hmmscan_results, hmm_lengths, remove_incomplete_only = True)
    refined_results = filter_by_score(refined_results, hmm_properties)
    refined_results = filter_overlaps(refined_results, hmm_properties)
    print(refined_results)
    cds_predictions = get_cds_predictions(refined_results, hmm_properties)

    return ProtoclusterPrediction(cds_predictions)
