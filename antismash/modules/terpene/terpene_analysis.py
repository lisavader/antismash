# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis section of the terpene module
"""
import math
import logging

from typing import Iterable, Optional
from collections import defaultdict

from antismash.common import fasta, path
from antismash.common.hmmscan_refinement import HMMResult, refine_hmmscan_results, QueryResult
from antismash.common.secmet import Protocluster, CDSFeature
from antismash.common.subprocessing.hmmscan import run_hmmscan, OutputType

from .results import ProtoclusterPrediction, DomainPrediction, Reaction, TerpeneHMM, CompoundGroup
from .load_json import load_json

_HMM_PROPERTIES_CACHE: dict[str, TerpeneHMM] = {}

_TOP_LEVEL_HMMS = {"PT_phytoene_like", "phytoene_synt," "Lycopene_cycl", "Lycopene_cycl_fung",
                    "PT_FPPS_like", "T1TS", "T1TS_KS", "T2TS", "TS_UbiA", "TS_Pyr4"}

_MAIN_TYPES_ORDER = ("PT_FPPS_like",
                     "PT_phytoene_like",
                     "T2TS",
                     "HAD_2",
                     "T1TS",
                     "TS_UbiA",
                     "TS_Pyr4",
                     "Lycopene_cycl")

def _load_hmm_properties() -> dict[str, TerpeneHMM]:
    """ Generates the HMM properties from hmm_properties.json and compound_groups.json
        Only does the processing once per python invocation, future runs access
        existing properties
    """
    # if already called once, then just reuse the cached results
    if not _HMM_PROPERTIES_CACHE:
        properties_json = load_json("hmm_properties")
        compounds_json = load_json("compound_groups")

        compound_groups: dict[str, CompoundGroup] = {}
        for group_data in compounds_json:
            compound_group = CompoundGroup.from_json(group_data)
            compound_groups[compound_group.name] = compound_group

        terpene_hmms: dict[str, TerpeneHMM] = {}
        for hmm_data in properties_json:
            terpene_hmm = TerpeneHMM.from_json(hmm_data, terpene_hmms, compound_groups)
            terpene_hmms[terpene_hmm.name] = terpene_hmm

        _HMM_PROPERTIES_CACHE.update(terpene_hmms)

    return _HMM_PROPERTIES_CACHE


def run_terpene_hmmscan(cds_features: Iterable[CDSFeature]) -> list[QueryResult]:
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


def filter_by_score(hmm_results_per_cds: dict[str, list[HMMResult]],
                    hmm_properties: dict[str, TerpeneHMM]) -> dict[str, list[HMMResult]]:
    """ Removes hmm results with bitscores that are below the specified cutoff for the hmm
    """
    results_by_id: dict[str, list[HMMResult]] = defaultdict(list)
    for cds_name, results in hmm_results_per_cds.items():
        # Store result if it is above the cut-off
        for result in results:
            if result.hit_id in hmm_properties:
                terpene_hmm = hmm_properties[result.hit_id]
            else:
                raise ValueError(f"Failed to find signature for ID {result.hit_id}")
            if result.bitscore > terpene_hmm.cutoff:
                results_by_id[cds_name].append(result)
    return results_by_id


def group_hmm_results(hmm_results: list[HMMResult]) -> list[list[HMMResult]]:
    """ Groups overlapping hmm results.
    """
    hmm_results = sorted(hmm_results, key=lambda result: result.query_start)
    allowed_overlap = 20
    groups = [[hmm_results[0]]]
    end = hmm_results[0].query_end
    for hmm_result in hmm_results[1:]:
        if hmm_result.query_start < (end - allowed_overlap):
            groups[-1].append(hmm_result)
        else:
            groups.append([hmm_result])
        end = hmm_result.query_end
    return groups


def filter_overlaps(hmm_results: list[HMMResult],
                    hmm_properties: dict[str, TerpeneHMM]) -> list[HMMResult]:
    """ Finds the best results within overlapping hits
    """
    results_by_name = {}
    for result in hmm_results:
        results_by_name[result.hit_id] = []
    for name in results_by_name:
        if name in _TOP_LEVEL_HMMS:
            print(name)
    raise NotImplementedError()

def remove_overlapping(hmm_results: list[HMMResult],
                       hmm_properties: dict[str, TerpeneHMM]) -> list[HMMResult]:
    """ Filters overlapping hmm results, keeping the hit with the lowest evalue.
        If the evalues are too similar, both hits are kept.
        Domains with an overlap of 20% or less aren't considered to be overlapping
    """
    non_overlapping = []
    if hmm_results:
        non_overlapping = [hmm_results[0]]
        for result in hmm_results[1:]:
            previous = non_overlapping[-1]
            maxoverlap = 0.20 * max(hmm_properties[result.hit_id].length,
                                hmm_properties[previous.hit_id].length)
            if result.query_start < (previous.query_end - maxoverlap):
                # Keep the current result if the current result scores significantly better
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


def merge_reactions_by_substrate(profiles: list[TerpeneHMM] #have hmms as input
                      ) -> tuple[Reaction, ...]:


    if not profiles:
        return tuple()

    # each tuple is a set of reactions from a single profile/HMM #move to tests
    # of which none will have multiple reactions with the same substrate


    # gather all reactions from all groups with the same substrate
    reactions_by_substrate: dict[tuple[CompoundGroup, ...], list[Reaction]] = defaultdict(list)
    for profile in profiles:
        for reaction in profile.reactions:
            reactions_by_substrate[reaction.substrates].append(reaction)

    # find the substrates which are present in all groups
    intersecting_groups = {sub: reactions for sub, reactions in reactions_by_substrate.items() if len(reactions) == len(reaction_groups)}
    # if no substrates are present in all groups, return nothing
    if not intersecting_groups:
        return tuple()

    # for each matching substrate, compress the reaction down to the intersection of all products
    results = []
    for substrate_group in intersecting_groups.values():
        first = substrate_group[0]
        for reaction in substrate_group[1:]:
            first = first.merge(reaction)
        # if no products remain, don't include the substrate
        if first.products:
            results.append(first)

    # return all remaining merged reactions
    return tuple(results)


def get_domain_prediction(hmm_results: list[HMMResult],
                          hmm_properties: dict[str, TerpeneHMM]) -> DomainPrediction:
    """ Converts a list of HMMResults to a DomainPrediction

        Arguments:
            hmm_results: a list of HMMResults

        Returns:
            domain_prediction: a DomainPrediction
    """
    start = min(hmm_result.query_start for hmm_result in hmm_results)
    end = max(hmm_result.query_end for hmm_result in hmm_results)
    main_types = set()
    subtypes = set()
    reaction_groups = []
    profiles = []
    for hmm_result in hmm_results:
        hmm = hmm_properties[hmm_result.hit_id]
        main_types.add(hmm.type)
        subtype_names = [subtype.name for subtype in hmm_properties[hmm.type].subtypes]
        if hmm_result.hit_id in subtype_names:
            subtypes.add(hmm_result.hit_id)
        if hmm.reactions:
            reaction_groups.append(hmm.reactions)
            profiles.append(hmm)
    if len(main_types) > 1:
        logging.warning("Overlapping hits for different main types.")
        type = "ambiguous hit"
        subtypes = []
    else:
        type = main_types.pop()
    final_reactions = merge_reactions_by_substrate(profiles)
    return DomainPrediction(type = type, subtypes = tuple(subtypes),
                            start = start, end = end,
                            reactions = final_reactions)


def get_cds_predictions(hmm_results_per_cds: dict[str, list[HMMResult]],
                        hmm_properties: dict[str, TerpeneHMM]) -> dict[str, list[DomainPrediction]]:
    """ Convert list of HMMResults in CDS mapping to a list of DomainPredictions

        Arguments:
            hmm_results: a mapping of CDS name to a list of HMMResults

        Returns:
            cds_predictions: a mapping of CDS name to a list of DomainPredictions
    """
    cds_predictions : dict[str, list[DomainPrediction]] = defaultdict(list)

    for cds_name, hmm_results in hmm_results_per_cds.items():
        grouped_results = group_hmm_results(hmm_results)
        # Add a domain prediction for each group
        for group in grouped_results:
            cds_predictions[cds_name].append(get_domain_prediction(group, hmm_properties))
            print(get_domain_prediction(group, hmm_properties))
    return cds_predictions


def get_cluster_prediction(cds_predictions: dict[str, list[DomainPrediction]]) -> ProtoclusterPrediction:
    cluster_pred = ProtoclusterPrediction(cds_predictions)
    all_domains = [domain for domains in cds_predictions.values() for domain in domains]
    # List domains that don't have a subtype first
    ordered_domains = sorted(all_domains, key=lambda domain: not domain.subtypes)
    # Then order by the predefined main type order
    ordered_domains = sorted(all_domains, key=lambda domain: _MAIN_TYPES_ORDER.index(domain.type))
    # Find the products of the last domain
    for reaction in ordered_domains[-1].reactions:
        for product in reaction.products:
            cluster_pred.add_product(product)
    return cluster_pred


def analyse_cluster(cluster: Protocluster) -> ProtoclusterPrediction:
    """ Analyse a terpene cluster

        Arguments:
            cluster: the Protocluster to analyse

        Returns:
            a single ProtoclusterPrediction instance with analysis results
    """
    assert cluster.product_category == "terpene"
    hmm_properties = _load_hmm_properties()
    hmm_lengths = {hmm_name: hmm_obj.length for hmm_name, hmm_obj in hmm_properties.items()} #make class HMM_Property with this as atr

    hmmscan_results = run_terpene_hmmscan(cluster.cds_children)
    refined_results = refine_hmmscan_results(hmmscan_results, hmm_lengths, remove_incomplete_only = True) #replace with remove_incomplete
    refined_results = filter_by_score(refined_results, hmm_properties)
    print(refined_results)
    #refined_results = filter_overlaps(refined_results, hmm_properties)
    cds_predictions = get_cds_predictions(refined_results, hmm_properties)
    return get_cluster_prediction(cds_predictions)
