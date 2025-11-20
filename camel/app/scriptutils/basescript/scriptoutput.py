import dataclasses
from pathlib import Path
from typing import Any

from camel.app.scriptutils import model


@dataclasses.dataclass(frozen=True)
class ScriptOutput(model.BaseOutput):
    """
    Contains the script output.
    """

    html: Path
    dir: Path
    tsv: Path | None = None
    json: Path | None = None
    fasta: Path | None = None

    def to_dict(self) -> dict[str, Any]:
        """
        Returns the output as a dictionary.
        :return: Dictionary with the output files.
        """
        data = {
            "html": str(self.html.absolute()),
            "dir": str(self.dir.absolute()),
        }
        if self.tsv is not None:
            data["tsv"] = str(self.tsv)
        if self.json is not None:
            data["json"] = str(self.json)
        return data
