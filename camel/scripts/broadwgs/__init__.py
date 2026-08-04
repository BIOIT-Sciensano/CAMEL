from pathlib import Path

_current_dir = Path(__file__).parent
CONFIG_DATA = _current_dir / 'config' / 'config_data.yml'
REFERENCES = _current_dir / 'config' / 'references.yml'
SLURM_SUBMIT = _current_dir / 'resources' / 'submit.py'
INTERVALS = _current_dir / 'resources' / 'intervalfiles'
SNAKEFILE_MAIN = _current_dir / 'snakefile' / 'main.smk'
TOOL_DATA = _current_dir / 'tool_data' / 'tool_data.yml'
