from pathlib import Path

_current_dir = Path(__file__).parent
CONFIG_DATA = str(_current_dir / 'config_data.yml')
SNAKEFILE_MAIN = str(_current_dir / 'snakefile' / 'main.smk')
