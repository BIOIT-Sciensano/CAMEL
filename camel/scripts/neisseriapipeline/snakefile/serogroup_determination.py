from pathlib import Path

SNAKEFILE = Path(__file__).parent / f'{Path(__file__).stem}.smk'

# Legacy detection
OUTPUT_LEGACY_REPORT = 'serogroup_determination/legacy/report/html.iob'
OUTPUT_LEGACY_REPORT_EMPTY = 'serogroup_determination/legacy/report/html-empty.iob'

# Capsule typing tool
OUTPUT_REPORT = 'serogroup_determination/capsule/report/html.iob'
OUTPUT_REPORT_EMPTY = 'serogroup_determination/capsule/report/html-empty.iob'
OUTPUT_INFORMS = 'serogroup_determination/capsule/tool/informs.io'

# Summary
OUTPUT_SUMMARY = 'serogroup_determination/summary/summary_out.{ext}'
