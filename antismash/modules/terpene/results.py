# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Contains the results classes for the terpene module """

import logging
from typing import Any, Dict, Optional

from antismash.common.module_results import ModuleResults
from antismash.common.secmet import Record


class TerpeneResults(ModuleResults):
    """ The combined results of the terpene module """
    _schema_version = 1

    def __init__(self, record_id: str) -> None:
        super().__init__(record_id)

    def to_json(self) -> Dict[str, Any]:
        results = {"schema_version": self._schema_version,
                   "record_id": self.record_id}
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