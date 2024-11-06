# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the terpene module """

import logging
from typing import Any, Dict, List, Optional, Union
from typing import Tuple
from dataclasses import dataclass

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
    initial_cyclisations: list[str]
    functional_groups: list[str]

    @staticmethod
    def from_json(name: str, data: dict[str, Any]) -> "CompoundGroup":
        """ Reconstructs an instance from a JSON representation """
        return CompoundGroup(name, **data)


@dataclass(frozen = True)
class TerpeneHMM:
    """ Properties associated with a terpene hmm profile
    """
    name: str
    description: str
    cutoff: int
    main_profile: bool
    predictions: dict[str, list[CompoundGroup]]

    def from_json(name: str, hmm_json: Dict[str, Any], compounds_json: dict[dict[str, Any]]) -> "TerpeneHMM":
        """ Reconstructs an instance from a JSON representation """
        predictions = {}
        compound_groups = []
        for pred in hmm_json["predictions"]:
            for field, group_names in pred.items():
                for group_name in group_names:
                    try:
                        compound_groups.append(CompoundGroup.from_json(group_name, compounds_json[group_name]))
                    except KeyError as err:
                        raise ValueError(f"Compound group {err} does not exist.")
            predictions[field] = compound_groups
        return TerpeneHMM(name, hmm_json["description"], hmm_json["cutoff"], hmm_json["main_profile"], predictions)


class DomainPrediction:
    """ A prediction for a terpene biosynthetic domain
    """
    def __init__(self, type: str, subtype: Optional[str], start: int, end: int) -> None:
        self.type = type
        self.subtype = subtype
        self.start = start
        self.end = end

    def __repr__(self) -> str:
        return f"DomainPrediction({self.type}, {self.subtype}, {self.start}, {self.end})"

    def __str__(self) -> str:
        parts = [
            f"Type: {self.type}",
            f"Subtype: {self.subtype}",
            f"Start: {self.start}",
            f"End: {self.end}"
        ]
        return "\n".join(parts)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, DomainPrediction):
            return False
        return (self.subtype == other.subtype
                and self.start == other.start
                and self.end == other.end)

    def to_json(self) -> Tuple[str, Optional[str], int, int]:
        """ Returns a JSON-friendly representation of the DomainPrediction """
        return self.type, self.subtype, self.start, self.end

    @staticmethod
    def from_json(json: List[Union[str, int]]) -> "DomainPrediction":
        """ Reconstructs a Prediction from a JSON representation """
        return DomainPrediction(str(json[0]), str(json[1]), int(json[2]), int(json[3]))

class ProtoclusterPrediction:
    """ A prediction for a terpene protocluster
    """

    def __init__(self, cds_predictions: Dict[str, List[DomainPrediction]],
                 start: int, end: int) -> None:
        self.cds_predictions = cds_predictions
        self.start = start
        self.end = end

    def __repr__(self) -> str:
        return f"Prediction({self.cds_predictions}, {self.start}, {self.end})"

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
                "start": self.start,
                "end": self.end
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
