from pathlib import Path

_current_dir = Path(__file__).parent
SNAKEFILE_SAMTOOLS_CALLING_ALL = _current_dir / 'snakefile'/ 'samtools_calling_all.smk'
SNAKEFILE_SAMTOOLS_FILTERING_ALL = _current_dir / 'snakefile'/ 'samtools_filtering_all.smk'
SNAKEFILE_TRIMMING_ALL = _current_dir / 'snakefile'/ 'trimming_all.smk'
