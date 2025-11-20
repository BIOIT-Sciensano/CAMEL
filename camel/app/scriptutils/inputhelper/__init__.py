import dataclasses
from typing import Callable

from camel.app.cli import cliutils
from camel.app.scriptutils import model
from camel.app.scriptutils.inputhelper.fastahelper import FastaHelper
from camel.app.scriptutils.inputhelper.illuminahelper import (
    IlluminaHelper,
    IlluminaTrimmingOpts,
    IlluminaAssemblyOpts,
)
from camel.app.scriptutils.inputhelper.inputhelperbase import TrimmingOpts, AssemblyOpts
from camel.app.scriptutils.inputhelper.onthelper import (
    ONTHelper,
    ONTTrimmingOpts,
    ONTAssemblyOpts,
)

helper_by_input_type = {
    model.InputType.FASTA: FastaHelper,
    model.InputType.ILLUMINA: IlluminaHelper,
    model.InputType.ONT: ONTHelper,
}

def add_helper_opts(f: Callable) -> Callable:
    """
    Adds the CLI options for the input helpers.
    :param f: Input function
    """
    f = cliutils.add_click_options_from_dataclass(IlluminaTrimmingOpts)(f)
    f = cliutils.add_click_options_from_dataclass(
        ONTTrimmingOpts, skip=[f.name for f in dataclasses.fields(TrimmingOpts)])(f)
    f = cliutils.add_click_options_from_dataclass(IlluminaAssemblyOpts)(f)
    f = cliutils.add_click_options_from_dataclass(
        ONTAssemblyOpts, skip=[f.name for f in dataclasses.fields(AssemblyOpts)])(f)
    return f
