import dataclasses
from pathlib import Path

import pytest

from camel.app.scriptutils import model
from camel.app.scriptutils.basepipe import basepipeutils
from camel.app.scriptutils.basescript.scriptinput import ScriptInput


@pytest.fixture
def fasta_input() -> ScriptInput:
    """
    Creates a FASTA input.
    :return: FASTA input
    """
    path_fasta = Path("path/to/assembly.fasta")
    return ScriptInput(type_=model.InputType.FASTA, fasta=path_fasta)


@pytest.fixture
def ont_input() -> ScriptInput:
    """
    Creates an ONT input.
    :return: ONT input
    """
    path_fastq = Path("path/to/reads.fastq.gz")
    return ScriptInput(
        type_=model.InputType.ONT, fastq_se=path_fastq, sample_name="ONT_sample"
    )


@pytest.fixture
def illumina_input() -> ScriptInput:
    """
    Creates an Illumina input.
    :return: Illumina input
    """
    paths_fastq = (Path("path/to/reads_1.fastq.gz"), Path("path/to/reads_2.fastq.gz"))
    return ScriptInput(
        type_=model.InputType.ILLUMINA,
        fastq_pe=paths_fastq,
        sample_name="Illumina_sample",
    )


def test_input_str_fasta(fasta_input: ScriptInput) -> None:
    """
    Tests the input string for FASTA input.
    :return: None
    """
    assert fasta_input.input_str == "assembly.fasta"


def test_input_str_fasta_with_name(fasta_input: ScriptInput) -> None:
    """
    Tests the input string for FASTA input with a custom name.
    :return: None
    """
    updated = dataclasses.replace(fasta_input, fasta_name="custom_name.fna")
    assert updated.input_str == "custom_name.fna"


def test_input_str_ont(ont_input: ScriptInput) -> None:
    """
    Tests the input string for ONT input.
    :return: None
    """
    assert ont_input.input_str == "reads.fastq.gz"


def test_input_str_ont_with_name(ont_input: ScriptInput) -> None:
    """
    Tests the input string for ONT input with a custom name.
    :return: None
    """
    updated = dataclasses.replace(ont_input, fastq_se_name="custom_name.fastq")
    assert updated.input_str == "custom_name.fastq"


def test_input_str_illumina(illumina_input: ScriptInput) -> None:
    """
    Tests the input string for Illumina input.
    :return: None
    """
    assert illumina_input.input_str == "reads_1.fastq.gz, reads_2.fastq.gz"


def test_name_returns_sample_name(ont_input: ScriptInput) -> None:
    """
    Tests that the name property returns the sample name.
    :return: None
    """
    assert ont_input.name == "ONT_sample"


def test_from_dict_roundtrip_fasta(fasta_input: ScriptInput) -> None:
    """
    Checks to_dict followed by from_dict.
    """
    data = fasta_input.to_dict()
    restored = ScriptInput.from_dict(data)

    for field in dataclasses.fields(ScriptInput):
        if field.name.endswith('_name'):
            continue
        val_orig = getattr(fasta_input, field.name)
        val_restored = getattr(restored, field.name)
        assert val_orig == val_restored, f"mismatched field: {field.name}"


@pytest.mark.parametrize("fixture_name", ["ont_input", "illumina_input", "fasta_input"])
def test_from_dict_roundtrip(request, fixture_name) -> None:
    """
    Checks to_dict followed by from_dict.
    :param request: Pytest request object
    :param fixture_name: Fixture name
    :return: None
    """
    instance: ScriptInput = request.getfixturevalue(fixture_name)
    data = instance.to_dict()
    restored = ScriptInput.from_dict(data)

    for field in dataclasses.fields(ScriptInput):
        if field.name.endswith('_name') or field.name.endswith("_names"):
            continue
        val_orig = getattr(instance, field.name)
        val_restored = getattr(restored, field.name)
        assert val_orig == val_restored, f"mismatched field: {field.name}"


def test_resolve_species_basic() -> None:
    """
    Tests basic resolution of a species-specific value.
    :return: None
    """
    cfg = {
        "species": {"faecalis": {"mlst_db": "/path/to/mlst"}},
        "tool": {"db": {"resolve": "species", "key": "mlst_db"}},
    }

    resolved = basepipeutils.resolve_config(cfg, "faecalis")
    assert resolved["tool"]["db"] == "/path/to/mlst"


def test_resolve_nested_dict() -> None:
    """
    Tests resolution within nested dictionary structures.
    :return: None
    """
    cfg = {
        "species": {"faecalis": {"ref_fasta": "/ref.fa"}},
        "reference": {"genome": {"fasta": {"resolve": "species", "key": "ref_fasta"}}},
    }

    resolved = basepipeutils.resolve_config(cfg, "faecalis")

    assert resolved["reference"]["genome"]["fasta"] == "/ref.fa"


def test_resolve_in_list() -> None:
    """
    Tests resolution of multiple values inside a list.
    :return: None
    """
    cfg = {
        "species": {"faecalis": {"a": "A", "b": "B"}},
        "values": [
            {"resolve": "species", "key": "a"},
            {"resolve": "species", "key": "b"},
        ],
    }
    resolved = basepipeutils.resolve_config(cfg, "faecalis")
    assert resolved["values"] == ["A", "B"]


def test_unknown_species() -> None:
    """
    Tests that resolving with an unknown species raises an error.
    :return: None
    """
    cfg = {"species": {}, "tool": {"db": {"resolve": "species", "key": "mlst_db"}}}

    with pytest.raises(ValueError, match="Unknown species"):
        basepipeutils.resolve_config(cfg, "faecalis")


def test_missing_key() -> None:
    """
    Tests that resolving a missing key for a species raises an error.
    :return: None
    """
    cfg = {
        "species": {"faecalis": {}},
        "tool": {"db": {"resolve": "species", "key": "mlst_db"}},
    }

    with pytest.raises(ValueError, match="Missing key"):
        basepipeutils.resolve_config(cfg, "faecalis")


def test_unknown_resolver() -> None:
    """
    Tests that an unknown resolver type raises an error.
    :return: None
    """
    cfg = {
        "species": {"faecalis": {"mlst_db": "/path"}},
        "tool": {"db": {"resolve": "unknown", "key": "mlst_db"}},
    }

    with pytest.raises(ValueError, match="Unknown resolver"):
        basepipeutils.resolve_config(cfg, "faecalis")


def test_mixed_values() -> None:
    """
    Tests resolution alongside static values in the same structure.
    :return: None
    """
    cfg = {
        "species": {"faecalis": {"mlst_db": "/mlst"}},
        "tool": {"db": {"resolve": "species", "key": "mlst_db"}, "threads": 4},
    }

    resolved = basepipeutils.resolve_config(cfg, "faecalis")

    assert resolved["tool"]["db"] == "/mlst"
    assert resolved["tool"]["threads"] == 4
