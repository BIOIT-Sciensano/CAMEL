from pathlib import Path

SNAKEFILE_VARIANTCALLING = f'{Path(__file__).parent / Path(__file__).stem}.smk'

_dir_varcalling = Path("variant_calling")
OUTPUT_gVCF = _dir_varcalling / "merge_vcf" / "output.gVCF.gz.io"
OUTPUT_gVCF_index = _dir_varcalling / "merge_vcf" / 'output.vcf.gz.tbi'
