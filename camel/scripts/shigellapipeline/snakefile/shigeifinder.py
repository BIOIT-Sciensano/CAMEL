from pathlib import Path


SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

_dir_shigeifinder = Path('shigeifinder')
OUTPUT_REPORT = 'shigeifinder/report/html.iob'
OUTPUT_REPORT_EMPTY = 'shigeifinder/report/html-empty.iob'
OUTPUT_SUMMARY = 'shigeifinder/summary_shigeifinder.{ext}'
OUTPUT_INFORMS = 'shigeifinder/tool/informs.io'
