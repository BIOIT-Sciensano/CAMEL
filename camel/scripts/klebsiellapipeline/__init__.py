from pathlib import Path

_current_dir = Path(__file__).parent

# Snakemake
SNAKEFILE_MAIN = _current_dir / 'snakefile' / 'main.smk'

# Other files
CONFIG_DATA = _current_dir / 'config_data.yml'
CITATIONS_HTML = _current_dir / 'citations.html'
