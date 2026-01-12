import dataclasses
from pathlib import Path

import click

from camel.app.scriptutils.model import BaseOptions


@dataclasses.dataclass(frozen=True)
class ScriptOptions(BaseOptions):
    """
    Contains the script options.
    """

    working_dir: Path = dataclasses.field(default=Path.cwd(), metadata={"help": "Working directory"})
    typing_method: str | None = dataclasses.field(
        default="mist",
        metadata={"choices": ["blast", "mist", "kma"], "show_default": True},
    )
    gene_detection_method: str | None = dataclasses.field(
        default="blast",
        metadata={"choices": ["blast", "kma"], "show_default": True},
    )
    trimming_method: str | None = dataclasses.field(
        default="fastp",
        metadata={"choices": ["fastp", "trimmomatic"], "show_default": True},
    )
    cov_max: int | None = dataclasses.field(
        default=100,
        metadata={
            "help": "Maximum coverage (datasets with higher estimated coverage are downsampled)",
            "show_default": True,
        },
    )
    ont_min_qual: int | None = dataclasses.field(
        default=10,
        metadata={
            "help": "Minimum quality of ONT reads",
            "show_default": True,
            "type": click.IntRange(0, 40),
        },
    )
    ont_min_length: int | None = dataclasses.field(
        default=1000,
        metadata={
            "help": "Minimum length of ONT reads",
            "show_default": True,
            "type": click.IntRange(min=0),
        },
    )
    threads: int | None = dataclasses.field(
        default=4, metadata={"help": "Nb. of threads", "show_default": True}
    )
    include_bam: bool = dataclasses.field(
        default=False, metadata={"help": "Include the BAM file in the output report"}
    )
    log: bool = dataclasses.field(
        default=False,
        metadata={"help": "Save logs to the directory specified in the config file"},
    )
    kraken2_small_db: bool = dataclasses.field(
        default=False, metadata={"help": "Use the small K2 database"}
    )
