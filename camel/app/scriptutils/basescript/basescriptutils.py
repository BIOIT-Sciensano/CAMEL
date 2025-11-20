from pathlib import Path
from typing import Callable, Any

import click

from camel.app.cli import cliutils
from camel.app.core.utils import fastqutils, fileutils
from camel.app.scriptutils import model
from camel.app.scriptutils.basescript.scriptinput import ScriptInput
from camel.app.scriptutils.basescript.scriptoptions import ScriptOptions
from camel.app.scriptutils.basescript.scriptoutput import ScriptOutput


def add_input_opts(*, supported: list[model.InputType] | None = None) -> Callable:
    """
    Creates a decorator to add pipeline input options to a Click command.
    :param supported: List of supported input types.
    :return: Decorator
    """

    def decorator(f: Callable) -> Callable:
        """
        Decorator to add pipeline input options to a Click command.
        :return: Decorated function.
        """
        options = [click.option("--sample-name", type=str, help="Sample name")]

        # Input type option is only needed if there are multiple input types
        if supported is None or len(supported) != 1:
            input_type_opts = [
                t.value
                for t in (supported if supported is not None else model.InputType)
            ]
            options.append(
                click.option(
                    "--input-type",
                    type=click.Choice(input_type_opts),
                    required=True,
                    help="Type of input data",
                )
            )

        # FASTA input
        if supported is None or model.InputType.FASTA in supported:
            options.append(
                click.option(
                    "--fasta",
                    type=click.Path(exists=True, path_type=Path),
                    help="FASTA input",
                )
            )
            options.append(
                click.option(
                    "--fasta-name", type=str, help="FASTA input name"
                )
            )

        # VCF file
        if supported is None or model.InputType.FASTA_WITH_VCF in supported:
            options.append(
                click.option(
                    "--vcf-unfiltered",
                    type=click.Path(exists=True, path_type=Path),
                    help="VCF input",
                )
            )

        # ONT input
        if supported is None or model.InputType.ONT in supported:
            options.append(
                click.option(
                    "--fastq-se",
                    type=click.Path(exists=True, path_type=Path),
                    help="ONT FASTQ input",
                )
            )
            options.append(
                click.option(
                    "--fastq-se-name", type=str, help="ONT FASTQ input name"
                )
            )

        # Illumina input
        if supported is None or model.InputType.ILLUMINA in supported:
            options.append(
                click.option(
                    "--fastq-pe",
                    nargs=2,
                    type=click.Path(exists=True, path_type=Path),
                    help="Paired-end FASTQ input (R1 R2)",
                )
            )
            options.append(
                click.option(
                    "--fastq-pe-names",
                    nargs=2,
                    type=str,
                    help="Paired-end FASTQ input names (R1_name R2_name)",
                )

            )

        for opt in reversed(options):
            f = opt(f)
        return f

    return decorator


def add_output_opts(f: Callable) -> Callable:
    """
    Creates a decorator to add pipeline output options to a Click command.
    :param f: Input function
    :return: Decorator
    """
    options = [
        click.option(
            "--output-html",
            type=click.Path(path_type=Path),
            default=Path("out/report.html"),
            help="Output report (HTML)",
            show_default=True,
        ),
        click.option(
            "--output-dir",
            type=click.Path(path_type=Path),
            default=Path("out"),
            help="Output directory",
            show_default=True,
        ),
        click.option(
            "--output-tsv",
            type=click.Path(path_type=Path),
            default=Path("out/summary.tsv"),
            help="Summary output (TSV)"
        ),
        click.option(
            "--output-json",
            type=click.Path(path_type=Path),
            help="Summary output (JSON)",
        ),
        click.option(
            "--output-fasta",
            type=click.Path(path_type=Path),
            help="Assembly output (FASTA)",
        ),
    ]
    for opt in reversed(options):
        f = opt(f)
    return f


def add_general_opts(f: Callable) -> Callable:
    """
    Adds the general script options.
    :param f: Input function
    :return: Decorator
    """
    return cliutils.add_click_options_from_dataclass(ScriptOptions)(f)


def extract_sample_name(kwargs: Any) -> str:
    """
    Extracts the sample name from the given arguments.
    :param kwargs: Arguments
    :return: Sample name
    """
    if kwargs.get("sample_name") is not None:
        return fileutils.make_valid(kwargs["sample_name"])
    # FASTA input
    input_type = model.InputType(kwargs["input_type"])
    if input_type in (model.InputType.FASTA, model.InputType.FASTA_WITH_VCF):
        if kwargs.get("fasta_name") is not None:
            return fileutils.make_valid(Path(kwargs["fasta_name"]).stem)
        return fileutils.make_valid(Path(kwargs["fasta"]).stem)

    # PE reads (illumina / hybrid)
    elif input_type in (model.InputType.ILLUMINA, model.InputType.HYBRID):
        if kwargs.get("fastq_pe_names") is not None:
            return fastqutils.get_sample_name(
                kwargs["fastq_pe_names"][0], fastqutils.PATTERN_FQ_PE
            )
        return fastqutils.get_sample_name(
            kwargs["fastq_pe"][0], fastqutils.PATTERN_FQ_PE
        )

    # SE reads (ONT)
    elif input_type is model.InputType.ONT:
        if kwargs.get("fastq_se_name") is not None:
            return fastqutils.get_sample_name(
                kwargs["fastq_se_name"], fastqutils.PATTERN_FQ_ONT
            )
        return fastqutils.get_sample_name(kwargs["fastq_se"], fastqutils.PATTERN_FQ_ONT)
    raise ValueError("Cannot determine sample name")


def parse_script_input(kwargs: Any) -> ScriptInput:
    """
    Parses the script input from the given arguments (provided by click).
    :param kwargs: Arguments
    :return: Parsed script input options
    """
    input_type = model.InputType(kwargs["input_type"])
    if input_type in (model.InputType.FASTA, model.InputType.FASTA_WITH_VCF):
        if kwargs.get("fasta") is None:
            raise click.UsageError(
                f"FASTA input is required for {input_type.value} input (--fasta)."
            )
    if input_type is model.InputType.FASTA_WITH_VCF:
        if kwargs.get("vcf_unfiltered") is None:
            raise click.UsageError(
                f"VCF input is required for {input_type.value} input (--vcf-unfiltered)."
            )
    if input_type is model.InputType.ILLUMINA:
        if kwargs.get("fastq_pe") is None:
            raise click.UsageError(
                f"Paired-end FASTQs are required for {input_type.value} input (--fastq-pe)."
            )
    if input_type is model.InputType.ONT:
        if kwargs.get("fastq_se") is None:
            raise click.UsageError(
                f"Single-end FASTQs are required for {input_type.value} input (--fastq-se)."
            )

    return ScriptInput(
        type_=model.InputType(input_type if input_type else kwargs["input_type"]),
        sample_name=extract_sample_name(kwargs),
        fasta=kwargs.get("fasta"),
        fasta_name=kwargs.get("fasta_name"),
        fastq_pe=kwargs.get("fastq_pe"),
        fastq_pe_names=kwargs.get("fastq_pe_names"),
        fastq_se=kwargs.get("fastq_se"),
        fastq_se_name=kwargs.get("fastq_se_name"),
        vcf_unfiltered=kwargs.get("vcf_unfiltered"),
    )


def parse_script_output(kwargs: Any) -> ScriptOutput:
    """
    Parses the script output arguments.
    :param kwargs: Arguments
    :return: Parsed script output options
    """
    return ScriptOutput(
        dir=Path(kwargs["output_dir"]).absolute(),
        html=Path(kwargs["output_html"]).absolute(),
        tsv=Path(kwargs["output_tsv"]).absolute()  if kwargs.get("output_tsv") else None,
        json=Path(kwargs["output_json"]).absolute() if kwargs.get("output_json") else None,
        fasta=Path(kwargs["output_fasta"]).absolute() if kwargs.get("output_fasta") else None,
    )


def parse_script_opts(kwargs) -> ScriptOptions:
    """
    Parses the script options.
    """
    return ScriptOptions(
        detection_method=kwargs["detection_method"],
        working_dir=(
            Path(kwargs["working_dir"]).absolute()
            if kwargs.get("working_dir")
            else None
        ),
        cov_max=kwargs["cov_max"],
        threads=kwargs.get("threads", 1),
    )
