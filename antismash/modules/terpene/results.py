# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the terpene module """

import logging
from typing import Any, Optional, Self
from dataclasses import dataclass, field
from functools import cached_property

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record

from .load_json import load_json

class MissingCompoundError(ValueError):
    pass

class MissingHmmError(ValueError):
    pass


@dataclass(frozen = True, slots = True)
class CompoundGroup():
    """ Biosynthetic and chemical properties for a group of compounds.
    """
    name: str
    extended_name: Optional[str]
    single_compound: bool
    biosynthetic_class: str
    biosynthetic_subclass: Optional[str]
    chain_length: int
    initial_cyclisations: tuple[str]
    functional_groups: tuple[str]

    def has_equal_skeleton(self, other: "CompoundGroup") -> bool:
        """ Returns whether this instance and the provided instance
            have the same chain length and initial cyclisations
        """
        return any([
            self.name == other.name,
            self.chain_length == other.chain_length \
            and (self.initial_cyclisations == other.initial_cyclisations
            or "unknown" in self.initial_cyclisations + other.initial_cyclisations)])

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-friendly representation """
        return {key: getattr(self, key) for key in self.__slots__}

    @classmethod
    def from_json(cls, data: dict[str, Any]) -> Self:
        """ Reconstructs an instance from a JSON representation """
        return cls(str(data["name"]), data["extended_name"], bool(data["single_compound"]),
                   str(data["biosynthetic_class"]), data["biosynthetic_subclass"],
                   int(data["chain_length"]), tuple(data["initial_cyclisations"]),
                   tuple(data["functional_groups"]))


@dataclass(frozen = True)
class Reaction():
    """ Contains the substrates and products of a chemical reaction.
    """
    substrates: tuple[CompoundGroup, ...]
    products: tuple[CompoundGroup, ...]

    def has_equal_substrates(self, other: "Reaction") -> bool:
        """ Returns wether this instance and the provided instance
            have the same substrates.
        """
        return self.substrates == other.substrates

    def merge(self, other: "Reaction") -> "Reaction":
        """ Creates a new Reaction instance that contains the
            intersections of the substrates and products of this instance
            and the provided instance.
        """
        return Reaction(tuple(set(self.substrates) & set(other.substrates)),
                                  tuple(set(self.products) & set(other.products)))

    def __str__(self) -> str:
        return (f"substrates={[compound.name for compound in self.substrates]}, "
                f"products={[compound.name for compound in self.products]}")

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-friendly representation  """
        return {"substrates": [compound.name for compound in self.substrates],
                "products": [compound.name for compound in self.products]}

    @classmethod
    def from_json(cls, data: dict[str, list[str]], compound_groups: dict[str, CompoundGroup]) -> Self:
        """ Reconstructs an instance from a JSON representation """
        try:
            substrates = data["substrates"]
            products = data["products"]
        except KeyError as key:
            raise MissingCompoundError(f"{data}: Field {key} is missing.")
        try:
            return cls(tuple(compound_groups[name] for name in substrates),
                       tuple(compound_groups[name] for name in products))
        except KeyError as key:
            raise MissingCompoundError(f"Compound group {key} is not defined.")


@dataclass
class TerpeneHMM:
    """ Properties associated with a terpene hmm profile
    """
    name: str
    description: str
    length: int
    cutoff: int
    subtypes: tuple["TerpeneHMM", ...]
    reactions: tuple[Reaction, ...]
    __parents: list[str] = field(default_factory = list)
    __main_type: str = ""

    @property
    def parents(self) -> list[str]:
        """ Returns the parents of the profile """
        return self.__parents

    @cached_property
    def main_type(self) -> str:
        """ Returns the main type of the profile """ #lookup name in name_to_type dict
        return self.__main_type

    def is_subtype(self) -> bool:
        """ Returns whether the profile is a subtype of a another profile """
        return bool(self.__parents)

    def add_parent(self, parent: str) -> None:
        """ Adds a parent to parents """
        self.__parents.append(parent)

    def add_main_type(self, main_type: str) -> None:
        """ Adds the main type of the profile, if it doesn't exist """
        if not self.__main_type:
            self.__main_type = main_type

    @classmethod
    def from_json(cls, hmm_json: dict[str, Any], terpene_hmms: dict[str, "TerpeneHMM"],
                  compound_groups: dict[str, CompoundGroup]) -> Self:
        """ Reconstructs an instance from a JSON representation """
        try:
            subtypes = tuple(terpene_hmms[name] for name in hmm_json["subtypes"])
        except KeyError as key:
            raise MissingHmmError(f"'{hmm_json['name']}': Subtype {key} not defined yet")
        for subtype in subtypes:
            subtype.add_parent(hmm_json["name"])
        return cls(str(hmm_json["name"]), str(hmm_json["description"]), int(hmm_json["length"]),
                   int(hmm_json["cutoff"]), subtypes,
                   tuple(Reaction.from_json(reaction, compound_groups)
                         for reaction in hmm_json["reactions"]))


@dataclass(slots = True)
class DomainPrediction:
    """ A prediction for a terpene biosynthetic domain
    """
    type: str
    subtypes: tuple[str, ...]
    start: int
    end: int
    reactions: tuple[Reaction, ...]

    def __str__(self) -> str:
        return(f"DomainPrediction(type={self.type}, subtypes={self.subtypes}, start={self.start}, "
               f"end={self.end}, reactions={[str(reaction) for reaction in self.reactions]}")

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-friendly representation """
        data = {key: getattr(self, key) for key in self.__slots__}
        reactions = data.pop("reactions")
        if reactions:
            data["reactions"] = [reaction.to_json() for reaction in reactions]
        return data

    @staticmethod
    def from_json(data: dict[str, Any], compound_groups: dict[str, CompoundGroup]) -> "DomainPrediction":
        """ Reconstructs an instance from a JSON representation """
        reactions = tuple(Reaction.from_json(reaction, compound_groups) for reaction in data.get("reactions", []))
        return DomainPrediction(str(data["type"]), tuple(data["subtypes"]), int(data["start"]),
                                int(data["end"]), reactions)

@dataclass
class ProtoclusterPrediction:
    """ A prediction for a terpene protocluster
    """
    cds_predictions: dict[str, list[DomainPrediction]]
    __products: list[CompoundGroup] = field(default_factory = list)

    @property
    def products(self) -> list[CompoundGroup]:
        """ Returns the predicted products of the cluster """
        return self.__products

    def add_product(self, product: CompoundGroup) -> None:
        """ Adds a product to products """
        self.__products.append(product)

    def __str__(self) -> str:
        parts = [
            f"CDSs: {len(self.cds_predictions)}",
        ]
        for cds, predictions in self.cds_predictions.items():
            parts.append(str(cds))
            parts.append(("\n ".join(map(str, predictions))))
        return "\n".join(parts)

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-friendly representation """
        cds_preds = {}
        for name, preds in self.cds_predictions.items():
            cds_preds[name] = [pred.to_json() for pred in preds]
        return {"cds_preds": cds_preds}

    @staticmethod
    def from_json(json: dict[str, Any], compound_groups: dict[str, CompoundGroup]
                  ) -> "ProtoclusterPrediction":
        """ Reconstructs an instance from a JSON representation """
        cds_preds: dict[str, list[DomainPrediction]] = {}
        for name, preds in json["cds_preds"].items():
            domains = [DomainPrediction.from_json(pred, compound_groups) for pred in preds]
            cds_preds[name] = domains
        return ProtoclusterPrediction(cds_preds)

class TerpeneResults(ModuleResults):
    """ The combined results of the terpene module """
    _schema_version = 1
    __slots__ = ["cluster_predictions"]

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)
        self.cluster_predictions: dict[int, ProtoclusterPrediction] = {}

    def __repr__(self) -> str:
        return f"TerpeneResults(clusters={list(self.cluster_predictions)})"

    def __str__(self) -> str:
        parts = []
        for cluster_id, prediction in self.cluster_predictions.items():
            parts.append(f"Protocluster {cluster_id}\n")
            parts.append(str(prediction))
        return "".join(parts)

    def to_json(self) -> dict[str, Any]:
        """ Returns a JSON-friendly representation """
        clusters = {cluster_number: pred.to_json() for cluster_number, pred in self.cluster_predictions.items()}
        results = {"schema_version": self._schema_version,
                   "record_id": self.record_id,
                   "protocluster_predictions": clusters}
        return results

    @staticmethod
    def from_json(json: dict[str, Any], record: Record,
                  ) -> Optional["TerpeneResults"]:
        """ Reconstructs an instance from a JSON representation """
        assert "record_id" in json
        if json.get("schema_version") != TerpeneResults._schema_version:
            logging.warning("Mismatching schema version, dropping terpene results")
            return None

        compound_groups: dict[str, CompoundGroup] = {}
        compounds_json = load_json("compound_groups")
        for group_data in compounds_json:
            compound_group = CompoundGroup.from_json(group_data)
            compound_groups[compound_group.name] = compound_group

        results = TerpeneResults(json["record_id"])
        try:
            for cluster_id, prediction in json["protocluster_predictions"].items():
                cluster_prediction = ProtoclusterPrediction.from_json(prediction, compound_groups)
                results.cluster_predictions[int(cluster_id)] = cluster_prediction
        except (MissingCompoundError, MissingHmmError) as err:
            logging.warning(f"Discarding results, missing referenced terpene data: {err}")
            return None
        return results

    def add_to_record(self, record: Record) -> None:
        """ Save terpene prediction in record.

            Cluster predictions are saved the relevant Cluster feature.
            Gene functions are added to each CDSFeature.
        """
        if record.id != self.record_id:
            raise ValueError("Record to store in and record analysed don't match")
