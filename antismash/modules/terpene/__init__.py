# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Terpene analysis module
"""

from typing import Any, Dict, List, Optional

from antismash.common.secmet import Record
from antismash.config import ConfigType
from antismash.config.args import ModuleArgs

from .results import TerpeneResults

NAME = "terpene"
SHORT_DESCRIPTION = "Terpene analysis"


def get_arguments() -> ModuleArgs:
    """ Return the args for the terpene module """
    args = ModuleArgs('Advanced options', 'terpene', enabled_by_default=True)
    return args


def is_enabled(options: ConfigType) -> bool:
    """ Whether the module is enabled """
    return not options.minimal or options.terpene_enabled


def check_options(_options: ConfigType) -> List[str]:
    """ No options to check at the moment """
    return []


def check_prereqs(options: ConfigType) -> List[str]:
    """ Check the prerequisites.
            hmmscan: domain detection
            hmmpress: compressing the hmm files

        Returns:
            a list of strings describing any errors, if they occurred
    """
    failure_messages = []
    for binary_name in ["hmmsearch", "hmmpress"]:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate executable for {binary_name!r}")

    return failure_messages


def regenerate_previous_results(previous: Dict[str, Any], record: Record,
                                _options: ConfigType) -> Optional[TerpeneResults]:
    """ Regenerate the previous results from JSON format.

        Arguments:
            previous: the previous results as from JSON
            record: the Record these previous results were originally created from
            options: the current antismash config object

        Returns:
            an instance of the module's ModuleResults implementation,
            or None if the current options require the analysis to be rerun or cannot be regenerated
    """
    # if there isn't anything to work with, just return None
    if not previous:
        return None
    return TerpeneResults.from_json(previous, record)


def run_on_record(record: Record, results: TerpeneResults, options: ConfigType) -> TerpeneResults:
    """ Run the analysis, unless the previous results apply to the given record

        Arguments:
            record: the Record being analysed
            results: an existing instance of the module's ModuleResults implementation (or None)
            options: the current antismash config object

        Returns:
            a TerpeneResults instance
    """
    if isinstance(results, TerpeneResults) and results.record_id == record.id:
        return results

    results = TerpeneResults(record.id)

    return results
