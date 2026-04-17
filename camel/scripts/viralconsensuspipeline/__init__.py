from importlib.resources import files
from pathlib import Path

_current_dir = Path(__file__).parent
CONFIG_DATA = _current_dir / 'config_data.yml'
SNAKEFILE_MAIN = Path(
    str(files('camel').joinpath('scripts/viralconsensuspipeline/snakefile/main.smk'))
)
