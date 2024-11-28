# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

# for test files, silence irrelevant and noisy pylint warnings
# pylint: disable=use-implicit-booleaness-not-comparison,protected-access,missing-docstring

import dataclasses
import unittest
from unittest.mock import Mock, patch

from antismash.common import json
from antismash.common.secmet.test.helpers import DummyRecord
from antismash.modules.terpene.results import (
    CompoundGroup,
    Reaction,
    TerpeneHMM,
    DomainPrediction,
    ProtoclusterPrediction,
    TerpeneResults,
    MissingCompoundError,
    MissingHmmError
)

## BUILD DUMMY DATA CLASSES ##

def build_dummy_compound_groups():
    return {
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

def build_dummy_reaction():
    return Reaction(substrates = (build_dummy_compound_groups()["GFPP"],),
                    products = (build_dummy_compound_groups()["C40_carotenoid"],))

def build_dummy_domain():
    return DomainPrediction(type = "T1TS", subtypes = ("T1TS_Bas_a",),
                         start = 1, end = 200, reactions = (build_dummy_reaction(),))

def build_dummy_cds_preds():
    return {"cds1": [build_dummy_domain()],
            "cds2": []}

def build_dummy_cluster_pred():
    return ProtoclusterPrediction(cds_predictions = build_dummy_cds_preds(),)

def build_dummy_terpene_hmms():
    return {"PT_FPP_bact": TerpeneHMM(
            name = "PT_FPP_bact",
            description = "Prenyltransferase; Farnesyl diphosphate synthase, bacterial",
            length = 265,
            cutoff = 250,
            subtypes = tuple(),
            reactions = (build_dummy_reaction(),))
            }

## FAKE DATA ##

fake_reaction_data = {"substrates" : ["Fake1","Fake2"],
                      "products" : ["Fake3"]}

fake_hmm_data = {"name": "PT_FPP_bact",
                "description": "Prenyltransferase; Farnesyl diphosphate synthase, bacterial",
                "length": 265,
                "cutoff": 250,
                "subtypes": ["Fake_subtype"],
                "reactions": [fake_reaction_data]}

fake_cluster_data = {"cds_preds": {},
                     "products": ["Fake1","Fake2"]}

class TestJSONConversion(unittest.TestCase):
    def test_compound_regeneration(self):
        compound = build_dummy_compound_groups()["GFPP"]
        regenerated = CompoundGroup.from_json(json.loads(json.dumps(compound.to_json())))
        assert regenerated == compound

    def test_reaction_regeneration(self):
        reaction = build_dummy_reaction()
        regenerated = Reaction.from_json(json.loads(json.dumps(reaction.to_json())), build_dummy_compound_groups())
        assert regenerated == reaction

    def test_domain_pred_regeneration(self):
        pred = build_dummy_domain()
        regenerated = DomainPrediction.from_json(json.loads(json.dumps(pred.to_json())), build_dummy_compound_groups())
        assert regenerated == pred

    def test_cluster_pred_regeneration(self):
        pred = build_dummy_cluster_pred()
        pred.add_product(build_dummy_compound_groups()["GFPP"])
        regenerated = ProtoclusterPrediction.from_json(json.loads(json.dumps(pred.to_json())), build_dummy_compound_groups())
        assert regenerated == pred

    def test_regeneration(self):
        record = DummyRecord()
        results = TerpeneResults(record.id)
        results.cluster_predictions = {1: build_dummy_cluster_pred()}
        regenerated = TerpeneResults.from_json(json.loads(json.dumps(results.to_json())), record)
        assert regenerated.cluster_predictions == results.cluster_predictions

    def test_missing_compounds(self):
        with self.assertRaises(MissingCompoundError):
            Reaction.from_json(fake_reaction_data, build_dummy_compound_groups())
        with self.assertRaises(MissingCompoundError):
            ProtoclusterPrediction.from_json(fake_cluster_data, build_dummy_compound_groups())
        with self.assertRaises(MissingHmmError):
            TerpeneHMM.from_json(fake_hmm_data, build_dummy_terpene_hmms(), build_dummy_compound_groups())

"""
class TestAnalysis(unittest.TestCase):
    def test_overlap_detection(self):
        self.fail()

    def test_cds_pred(self):
        with patch.object(terpene, "run_terpene_hmmscan", return_value=["stuff", "bob"]):
            res = terpene.terpene_analysis.get_cds_predictions(..)
        assert res.blah == ksjlkfdj

"""