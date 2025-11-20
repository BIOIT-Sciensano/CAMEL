from pathlib import Path

_current_dir = Path(__file__).parent
SNAKEFILE_MAIN = _current_dir / 'snakefile' / 'main.smk'
CONFIG_DATA = _current_dir / 'config_data.yml'
AMR_CIRCOS_TEMPLATE = _current_dir / 'static' / 'amr_circos_template.txt'
