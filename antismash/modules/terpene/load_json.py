
from typing import Any

from antismash.common import json, path

def load_json(json_file: str) -> list[dict[str, Any]]:
    """ Load a json file that's present in data """
    file_path = path.get_full_path(__file__, "data", json_file + ".json")
    with open(file_path) as handle:
        try:
            return json.load(handle)
        except ValueError as error:
            raise ValueError(f"{file_path!r} is not a valid json file: {error}")