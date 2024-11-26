# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The analysis section of the terpene module
"""
import math

from typing import Iterable
from collections import defaultdict

from antismash.common import fasta, path
from antismash.common.hmmscan_refinement import HMMResult, refine_hmmscan_results, QueryResult
from antismash.common.secmet import Protocluster, CDSFeature
from antismash.common.subprocessing.hmmscan import run_hmmscan, OutputType

from .results import ProtoclusterPrediction, DomainPrediction, Reaction, TerpeneHMM, CompoundGroup
from .load_json import load_json

_HMM_PROPERTIES_CACHE: dict[str, TerpeneHMM] = {}

_TYPE_MAPPINGS = {"phytoene_synt": "PT_phytoene_like",
                  "Lycopene_cycl_fung": "Lycopene_cycl"}

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

        # Retrieve the main type from parents
        def map_name(name: str) -> str:
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
            assert main_type in _MAIN_TYPES_ORDER, f"{main_type} is not an accepted main type"
            terpene_hmm.add_main_type(main_type)

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


def filter_overlaps(hmm_results_per_cds: dict[str, list[HMMResult]],
                   hmm_properties: dict[str, TerpeneHMM]) -> dict[str, list[HMMResult]]:
    """ Removes overlaps within main profile hits and within subtype hits
    """
    filtered_results: dict[str, list[HMMResult]] = defaultdict(list)
    # Remove overlapping results within subtypes and within main types
    for cds_name, results in hmm_results_per_cds.items():
        subtypes: list[HMMResult] = []
        main_types: list[HMMResult] = []
        for result in results:
            if hmm_properties[result.hit_id].subtypes:
                main_types.append(result)
            else:
                subtypes.append(result)
        filtered_results[cds_name].extend(remove_overlapping(subtypes, hmm_properties) +
                                       remove_overlapping(main_types, hmm_properties))
    return filtered_results


def merge_reactions_by_substrate(reactions_per_hmm: list[tuple[Reaction, ...]]
                      ) -> tuple[Reaction, ...]:
    """ When one domain has multiple hmm hits,
        merge the associated reactions by taking the intersect of their products.
        Reactions are only merged if their substrate(s) is/are equal.
        If none of the reactions share substrates, return an empty tuple.
    """
    if not reactions_per_hmm:
        return tuple()
    final_reactions: list[tuple[Reaction, ...]] = [reactions_per_hmm[0]]
    if len(reactions_per_hmm) > 1:
        for reactions in reactions_per_hmm[1:]:
            merged_reactions = tuple(reaction.merge(final_reaction)
                            for reaction in reactions
                            for final_reaction in final_reactions[-1]
                            if reaction.has_equal_substrates(final_reaction))
            if merged_reactions:
                final_reactions[-1] = merged_reactions
            #If two hits have contrasting predictions, don't predict anything
            else:
                return tuple()
    return tuple(reaction for reactions in final_reactions for reaction in reactions)

"""
substrate_to_reaction: dict[substrates, (group count, Reaction)]
for reaction_group in reactions_per_hmm
    for reac in group
        if reaction.substrates not in substrate_to_reaction
            substrate_to_reaction{substrate} = reaction
            continue
    count, merged = substrate_to_reaction{substrate}
    substrate_to_reaction{substrate} = count += 1, reaction.merge(merged)

return tuple(sub for sub in substrate_to_reaction.values() if sub.products)
"""

def get_domain_prediction(hmm_results: list[HMMResult],
                          hmm_properties: dict[str, TerpeneHMM]) -> DomainPrediction:
    """ Converts a list of HMMResults to a DomainPrediction

        Arguments:
            hmm_results: a list of HMMResults

        Returns:
            domain_prediction: a DomainPrediction
    """
    start_locations = []
    end_locations = []
    main_types = set()
    subtypes = set()
    reactions_per_hmm = []
    for hmm_result in hmm_results:
        start_locations.append(hmm_result.query_start) #no list just int
        end_locations.append(hmm_result.query_end)
        terpene_hmm = hmm_properties[hmm_result.hit_id]
        main_types.add(terpene_hmm.main_type)
        if terpene_hmm.is_subtype():
            subtypes.add(terpene_hmm.name)
        if terpene_hmm.reactions:
            reactions_per_hmm.append(terpene_hmm.reactions)
    assert len(main_types) == 1, "Overlapping hits cannot belong to different main types" #main_types ambiguous and subtypes unknown
    final_reactions = merge_reactions_by_substrate(reactions_per_hmm)

    return DomainPrediction(type = main_types.pop(), subtypes = tuple(subtypes),
                                start = min(start_locations), end = max(end_locations),
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
        hmm_results = sorted(hmm_results, key=lambda result: result.query_start)
        # Extract groups of overlapping hmm results
        groups = [[hmm_results[0]]]
        for result in hmm_results[1:]:
            overlapping = False
            for member in groups[-1]:
                maxoverlap = 0.20 * max(hmm_properties[result.hit_id].length,
                                        hmm_properties[member.hit_id].length)
                if result.query_start < (member.query_end - maxoverlap):
                    overlapping = True
                    break
            if overlapping:
                groups[-1].append(result)
            else:
                groups.append([result])

        # Add a domain prediction for each group
        for group in groups:
            cds_predictions[cds_name].append(get_domain_prediction(group, hmm_properties))
    return cds_predictions


def get_cluster_prediction(cds_predictions: dict[str, list[DomainPrediction]]) -> ProtoclusterPrediction:
    cluster_pred = ProtoclusterPrediction(cds_predictions)
    all_domains = [domain for domains in cds_predictions.values() for domain in domains]
    # List domains that don't have a subtype first
    ordered_domains = sorted(all_domains, key=lambda domain: not domain.subtypes)
    # Then order by the predefined main type order
    ordered_domains = sorted(all_domains, key=lambda domain: _MAIN_TYPES_ORDER.index(domain.type))
    # Find the products of the last domain
    cluster_products = [reaction.products for reaction in ordered_domains[-1].reactions]
    if cluster_products:
        for product in cluster_products:
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
    refined_results = filter_overlaps(refined_results, hmm_properties)
    cds_predictions = get_cds_predictions(refined_results, hmm_properties)
    return get_cluster_prediction(cds_predictions)
