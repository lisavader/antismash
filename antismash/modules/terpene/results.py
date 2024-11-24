# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the terpene module """

import logging
from typing import Any, Dict, List, Optional, Self, Union
from typing import Tuple
from dataclasses import dataclass, field
from collections import defaultdict

from antismash.common.json import JSONBase
from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record


@dataclass(frozen = True)
class CompoundGroup():
    """ Biosynthetic and chemical properties for a group of compounds.
    """
    name: str
    extended_name: str
    single_compound: bool
    biosynthetic_class: str
    biosynthetic_subclass: str
    chain_length: int
    initial_cyclisations: tuple[str]
    functional_groups: tuple[str]

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstructs an instance from a JSON representation """
        return cls(data["name"], data["extended_name"], data["single_compound"],
                   data["biosynthetic_class"], data["biosynthetic_subclass"],
                   data["chain_length"], tuple(data["initial_cyclisations"]),
                   tuple(data["functional_groups"]))


@dataclass(frozen = True)
class ReactionPrediction():
    """ Prediction of the substrates and products of a reaction.
    """
    substrates: tuple[CompoundGroup, ...]
    products: tuple[CompoundGroup, ...]

    def has_equal_substrates(self, other: "ReactionPrediction") -> bool:
        """ Returns wether this instance and the provided instance
            have the same substrates.
        """
        return self.substrates == other.substrates

    def merge(self, other: "ReactionPrediction") -> "ReactionPrediction":
        """ Creates a new ReactionPrediction instance that contains the
            intersections of the substrates and products of this instance
            and the provided instance.
        """
        return ReactionPrediction(tuple(set(self.substrates) & set(other.substrates)),
                                  tuple(set(self.products) & set(other.products)))

    def __str__(self) -> str:
        return (f"substrates={[compound_group.name for compound_group in self.substrates]}, "
                f"products={[compound_group.name for compound_group in self.products]}")

    @classmethod
    def from_json(cls, data: dict[str, list[str]], compound_groups: dict[str, CompoundGroup]) -> Self:
        """ Reconstructs an instance from a JSON representation """
        try:
            substrates = data["substrates"]
            products = data["products"]
        except KeyError as key:
            raise ValueError(f"{data}: Field {key} is missing")
        try:
            return cls(tuple(compound_groups[name] for name in substrates),
                       tuple(compound_groups[name] for name in products))
        except KeyError as key:
            raise ValueError(f"Compound group {key} is not defined")


@dataclass
class TerpeneHMM:
    """ Properties associated with a terpene hmm profile
    """
    name: str
    description: str
    length: int
    cutoff: int
    subtypes: tuple["TerpeneHMM", ...]
    predictions: tuple[ReactionPrediction, ...]
    __parents: list[str] = field(default_factory = list)
    __main_type: str = ""

    @property
    def parents(self) -> list:
        """ Returns the parents of the profile """
        return self.__parents

    @property
    def main_type(self) -> str:
        """ Returns the main type of the profile """
        return self.__main_type

    def is_subtype(self) -> bool:
        """ Returns whether the profile is a subtype of a another profile """
        return bool(self.__parents)

    def add_parent(self, parent: str) -> None:
        """ Adds a parent to parents """
        self.__parents.append(parent)

    def add_main_type(self, main_type: str) -> None:
        """ Adds the main type of the profile (if not defined yet) """
        if not self.__main_type:
            self.__main_type = main_type

    @classmethod
    def from_json(cls, hmm_json: Dict[str, Any], terpene_hmms: dict[str, "TerpeneHMM"],
                  compound_groups: dict[str, CompoundGroup]) -> Self:
        """ Reconstructs an instance from a JSON representation """
        try:
            subtypes = tuple(terpene_hmms[name] for name in hmm_json["subtypes"])
        except KeyError as key:
            raise ValueError(f"'{hmm_json['name']}': Subtype {key} not defined yet")
        for subtype in subtypes:
            subtype.add_parent(hmm_json["name"])
        return cls(hmm_json["name"], hmm_json["description"], hmm_json["length"], hmm_json["cutoff"], subtypes,
                   tuple(ReactionPrediction.from_json(pred, compound_groups) for pred in hmm_json["predictions"]))

@dataclass
class DomainPrediction:
    """ A prediction for a terpene biosynthetic domain
    """

    type: str
    subtypes: tuple[str, ...]
    start: int
    end: int
    predictions: tuple[ReactionPrediction, ...]

    def __str__(self) -> str:
        return(f"DomainPrediction(type={self.type}, subtypes={self.subtypes}, start={self.start}, "
               f"end={self.end}, predictions={[str(prediction) for prediction in self.predictions]}")

    def to_json(self) -> Tuple[str, Optional[str], int, int]:
        """ Returns a JSON-friendly representation of the DomainPrediction """
        data = {key: getattr(self, key) for key in self.__slots__}
        predictions = data.pop("predictions")
        if predictions:
            for pred in predictions:
                data["predictions"].append({"substrates": [compound_group.name for compound_group in pred.substrates],
                                            "products": [compound_group.name for compound_group in pred.products]})
        return data

    @staticmethod
    def from_json(data: Dict[str, Any]) -> "DomainPrediction":
        """ Reconstructs a DomainPrediction from a JSON representation """
        predictions = [ReactionPrediction.from_json(pred) for pred in data.get("predictions", [])]
        return DomainPrediction(str(json[0]), str(json[1]), int(json[2]), int(json[3]))

@dataclass
class ProtoclusterPrediction:
    """ A prediction for a terpene protocluster
    """

    def __init__(self, cds_predictions: Dict[str, List[DomainPrediction]]) -> None:
        self.cds_predictions = cds_predictions

    def __str__(self) -> str:
        parts = [
            f"CDSs: {len(self.cds_predictions)}",
        ]

        for cds, predictions in self.cds_predictions.items():
            parts.append(str(cds))
            parts.append(" " + ("\n ".join(map(str, predictions))))
        return "\n".join(parts)

    def to_json(self) -> Dict[str, Any]:
        """ Converts a ProtoclusterPrediction into a JSON friendly format """
        cds_preds = {}
        for name, preds in self.cds_predictions.items():
            cds_preds[name] = [pred.to_json() for pred in preds]

        return {"cds_preds": cds_preds,
                }

    @staticmethod
    def from_json(json: Dict[str, Any]) -> "ProtoclusterPrediction":
        """ Rebuilds a ProtoclusterPrediction from JSON """
        assert isinstance(json, dict), json
        cds_predictions = {name: list(map(DomainPrediction.from_json, preds)) for name, preds in json["cds_preds"].items()}
        return ProtoclusterPrediction(cds_predictions, int(json["start"]), int(json["end"]))

class TerpeneResults(ModuleResults):
    """ The combined results of the terpene module """
    _schema_version = 1
    __slots__ = ["cluster_predictions"]

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self.cluster_predictions: Dict[int, ProtoclusterPrediction] = {}

    def __repr__(self) -> str:
        return f"TerpeneResults(clusters={list(self.cluster_predictions)})"

    def __str__(self) -> str:
        parts = []
        for cluster_id, prediction in self.cluster_predictions.items():
            parts.append(f"Protocluster {cluster_id}\n")
            parts.append(str(prediction))
        return "".join(parts)

    def to_json(self) -> Dict[str, Any]:
        clusters = {cluster_number: pred.to_json() for cluster_number, pred in self.cluster_predictions.items()}
        results = {"schema_version": self._schema_version,
                   "record_id": self.record_id,
                   "protocluster_predictions": clusters}
        return results

    @staticmethod
    def from_json(json: Dict[str, Any], _record: Record) -> Optional["TerpeneResults"]:
        assert "record_id" in json
        if json.get("schema_version") != TerpeneResults._schema_version:
            logging.warning("Mismatching schema version, dropping terpene results")
            return None
        results = TerpeneResults(json["record_id"])

        return results

    def add_to_record(self, record: Record) -> None:
        """ Save terpene prediction in record.

            Cluster predictions are saved the relevant Cluster feature.
            Gene functions are added to each CDSFeature.
        """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
