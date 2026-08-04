from pathlib import Path

_current_dir = Path(__file__).parent
CONFIG_DATA = _current_dir / 'config_base.yml'
SNAKEFILE_MAIN = _current_dir / 'snakefile' / 'main.smk'
TSV_BASECALLING_MODELS = _current_dir / 'basecalling_models.tsv'
