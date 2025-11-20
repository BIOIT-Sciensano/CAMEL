from pathlib import Path

_current_dir = Path(__file__).parent.resolve()
CONFIG_DATA = _current_dir / 'config_data.yml'
SNAKEFILE_MAIN = _current_dir / 'snakefile' / 'main.smk'
