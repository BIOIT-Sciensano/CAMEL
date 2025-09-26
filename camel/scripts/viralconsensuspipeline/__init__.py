from importlib.resources import files
from pathlib import Path

SNAKEFILE_MAIN = Path(str(files('camel').joinpath('scripts/viralconsensuspipeline/snakefile/main.smk')))
