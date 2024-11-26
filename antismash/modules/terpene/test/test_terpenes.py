# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import unittest
from unittest.mock import Mock, patch

from antismash.common import json
from antismash.common.secmet.test.helpers import DummyRecord
from antismash.modules.terpene.results import (
    CompoundGroup,
    Reaction,
    DomainPrediction,
    ProtoclusterPrediction,
    TerpeneResults
)

existing_compound_groups = {
    "GFPP": CompoundGroup(
        name = "GFPP",
        extended_name = "geranylfarnesyl diphosphate",
        single_compound = True,
        biosynthetic_class = "sesterterpene",
        biosynthetic_subclass = None,
        chain_length = 25,
        initial_cyclisations = tuple(),
        functional_groups = ("PP",)),

    "C40_carotenoid": CompoundGroup(
        name = "C40_carotenoid",
        extended_name = None,
        single_compound = False,
        biosynthetic_class = "tetraterpene",
        biosynthetic_subclass = "C40 carotenoid",
        chain_length = 40,
        initial_cyclisations = ("unknown",),
        functional_groups = ("unknown",))
}

dummy_compound_group = CompoundGroup(
    name = "Compound",
    extended_name = "Longname",
    single_compound = True,
    biosynthetic_class = "diterpene",
    biosynthetic_subclass = "indole diterpenoid",
    chain_length = 20,
    initial_cyclisations = ["C6-C11", "C10-C15"],
    functional_groups = ["OH", "PP", "epoxy", "indole"])

dummy_reaction = Reaction(substrates = (existing_compound_groups["GFPP"],),
                          products = (existing_compound_groups["C40_carotenoid"],))

dummy_domain1 = DomainPrediction(type = "T1TS", subtypes = ("T1TS_Bas_a",),
                         start = 1, end = 200, reactions = (dummy_reaction,))
dummy_domain2 = DomainPrediction(type = "PT_phytoene_like", subtypes = tuple(),
                               start = 220, end = 400, reactions = tuple())

dummy_cds_preds = {
"cds1": [dummy_domain1, dummy_domain2],
"cds2": []
}

dummy_cluster_pred = ProtoclusterPrediction(cds_predictions = dummy_cds_preds)

class TestJSONConversion(unittest.TestCase):
    def test_regeneration(self):
        record = DummyRecord()
        results = TerpeneResults(record.id)
        results.cluster_predictions = {1: dummy_cluster_pred}
        regenerated = TerpeneResults.from_json(json.loads(json.dumps(results.to_json())), record)
        assert regenerated.cluster_predictions == results.cluster_predictions
"""
    def test_cluster_prediction(self):
        self.fail()


    def test_missing_compounds(self):
        pred = Predction()
        with self.assertRaises(MissingCompoundError):
            pred.from_json()
        full_results = Results(pred)

        res = Results.from_json(pred.to_json())
        assert res is None
"""

"""
class TestAnalysis(unittest.TestCase):
    def test_overlap_detection(self):
        self.fail()

    def test_cds_pred(self):
        with patch.object(terpene, "run_terpene_hmmscan", return_value=["stuff", "bob"]):
            res = terpene.terpene_analysis.get_cds_predictions(..)
        assert res.blah == ksjlkfdj

"""