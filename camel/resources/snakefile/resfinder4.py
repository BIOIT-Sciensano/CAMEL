from pathlib import Path

DIR = Path('resfinder4')
SNAKEFILE = f'{Path(__file__).parent / Path(__file__).stem}.smk'

# Report and summary
OUTPUT_INFORMS = DIR / 'informs.io'
OUTPUT_REPORT = DIR / 'html.iob'
OUTPUT_REPORT_EMPTY = DIR / 'html-empty.iob'
OUTPUT_SUMMARY = DIR / 'summary_resfinder.{ext}'
