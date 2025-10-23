import json
from pathlib import Path

from camel.app.core.command import Command
from camel.app.scriptutils.fastqinput import FastqInput
from camel.app.core.io.tooliofile import ToolIOFile
from camel.app.core.snakemake import snakemakeutils


rule samtools_ampliconclip:
    """
    Clips the primer sequences.
    """
    input:
        BAM = config['input']['bam']
    output:
        BAM = 'ampliconclip/alignment_clipped.bam',
        INFORMS = 'ampliconclip/informs.json'
    params:
        dir_ = lambda wildcards: f'ampliconclip',
        path_bed = config.get('bed_primers')
    run:
        import re

        # Create working directory
        dir_ = Path(str(params.dir_)).absolute()
        dir_.mkdir(exist_ok=True, parents=True)

        # Run command
        command = Command(' '.join([
            'module load samtools/1.17;',
            f'samtools ampliconclip --both-ends -b {Path(params.path_bed).absolute()} {Path(input.BAM).absolute()} |',
            f'samtools sort - -o {Path(output.BAM).absolute()};',
            f'samtools index {Path(output.BAM).absolute()};'
        ]))
        command.run(dir_)
        if not command.returncode == 0:
            raise RuntimeError(f'Error removing primers: {command.stderr}')

        # Collect amplicon clipping steps (if performed)
        stats_out = {}
        for line in command.stderr.splitlines():
            m = re.match('^([A-Z ]+): (\d+)', line)
            if not m:
                continue
            stats_out[m.group(1)] = int(m.group(2))

        # Save the informs
        with open(output.INFORMS, 'w') as handle:
            json.dump({
                '_name': 'samtools ampliconclip',
                '_name_full': 'samtools ampliconclip 1.17',
                '_version': '1.17',
                '_command': command.command,
                'stats': stats_out
            }, handle, indent=2)

rule create_bam_for_seq_id:
    """
    Creates a BAM file for the given sequence id.
    """
    input:
        BAM = config['input']['bam'] if config.get('primers_bed') is None else rules.samtools_ampliconclip.output.BAM
    output:
        BAM = 'by_seq_id/{seq_id}/{seq_id}.bam',
        INFORMS = 'by_seq_id/{seq_id}/informs.json'
    params:
        dir_ = lambda wildcards: f'by_seq_id/{wildcards.seq_id}',
        seq_id = lambda wildcards: wildcards.seq_id
    run:
        # Create working directory
        dir_ = Path(str(params.dir_)).absolute()
        dir_.mkdir(exist_ok=True, parents=True)

        # Execute the command
        command = Command(' '.join([
            'module load samtools/1.17;',
            f'samtools view -b {Path(input.BAM).absolute()} "{params.seq_id}" > {Path(output.BAM).absolute()};',
            f'samtools index {Path(output.BAM).absolute()};'
        ]))
        command.run(dir_)
        if not command.returncode == 0:
            raise RuntimeError(f'Error extracting BAM file: {command.stderr}')
        with open(output.INFORMS, 'w') as handle:
            json.dump({
                '_name': 'samtools view',
                '_name_full': 'samtools view 1.17',
                '_version': '1.17',
                '_command': command.command
            }, handle, indent=2)

rule downsample_bam:
    """
    Down-samples the BAM file to the required coverage (if needed).
    """
    input:
        BAM = rules.create_bam_for_seq_id.output.BAM
    output:
        BAM = 'by_seq_id/{seq_id}/{seq_id}-ds.bam'
    params:
        dir_ = lambda wildcards: f'by_seq_id/{wildcards.seq_id}',
        ratio = lambda wildcards: config['downsampling'][wildcards.seq_id]
    run:
        dir_ = Path(str(params.dir_)).absolute()
        dir_.mkdir(exist_ok=True, parents=True)
        if params.ratio is not None:
            command = Command(' '.join([
                'module load samtools/1.17;',
                f'samtools view -b -s {params.ratio} {Path(input.BAM).absolute()} > {Path(output.BAM).absolute()};'
            ]))
            command.run(dir_)
            if not command.returncode == 0:
                raise RuntimeError(f'Error downsampling BAM file: {command.stderr}')
        else:
            Path(output.BAM).absolute().symlink_to(Path(input.BAM).absolute())

rule extract_fq_ont:
    """
    Extracts FASTQ files for single-end data.
    """
    input:
        BAM = rules.downsample_bam.output.BAM
    output:
        FASTQ = 'by_seq_id/{seq_id}/fastq/ont/{seq_id}.fastq',
        INFORMS = 'by_seq_id/{seq_id}/fastq/ont/informs_{seq_id}.json'
    params:
        dir_ = lambda wildcards: f'by_seq_id/{wildcards.seq_id}/fastq/ont'
    run:
        dir_ = Path(str(params.dir_)).absolute()
        dir_.mkdir(exist_ok=True, parents=True)
        command = Command(' '.join([
            'module load samtools/1.17;',
            f'samtools fastq -0 {Path(output.FASTQ).absolute()} {Path(input.BAM).absolute()};'
        ]))
        command.run(dir_)
        if not command.returncode == 0:
            raise RuntimeError(f'Error extracting FASTQ files: {command.stderr}')
        with open(output.INFORMS, 'w') as handle:
            json.dump({
                '_name': 'samtools fastq',
                '_name_full': 'samtools fastq 1.17;',
                '_version': '1.17',
                '_command': command.command
            }, handle, indent=2)

rule extract_fq_ont_merge:
    """
    Merges the extracted SE FASTQ files. 
    """
    input:
        FASTQ = [str(rules.extract_fq_ont.output.FASTQ).format(seq_id=seq_id) for seq_id in config['downsampling'].keys()]
    output:
        FQ_dict = 'merged/ont/fq_dict.io'
    params:
        dir_ = 'merged/ont'
    run:
        dir_ = Path(params.dir_).absolute()
        dir_.mkdir(parents=True, exist_ok=True)
        fq_out = dir_ / 'merged.fastq.gz'
        command = Command(f"cat {' '.join(str(Path(x).absolute()) for x in input.FASTQ)} | gzip > {fq_out};")
        command.run(dir_)
        if not command.returncode == 0:
            raise RuntimeError(f'Error merged FASTQ files: {command.returncode}')
        fq_dict = FastqInput('nanopore', se=[ToolIOFile(fq_out)], is_trimmed=True, is_pe=False)
        snakemakeutils.dump_object(fq_dict.to_fq_dict(), Path(output.FQ_dict))

rule extract_fq_illumina:
    """
    Extracts FASTQ files for paired-end Illumina data.
    """
    input:
        BAM = rules.downsample_bam.output.BAM
    output:
        FASTQ = 'by_seq_id/{seq_id}/fastq/illumina/{seq_id}_fq_dict.io',
        INFORMS = 'by_seq_id/{seq_id}/fastq/illumina/informs_{seq_id}.json'
    params:
        dir_ = lambda wildcards: f'by_seq_id/{wildcards.seq_id}/fastq/illumina'
    run:
        dir_ = Path(str(params.dir_)).absolute()
        dir_.mkdir(exist_ok=True, parents=True)

        ################################################
        # Output files                                 #
        # 1P - Paired forward reads                    #
        # 2P - Paired reverse reads                    #
        # single - Reads with mate not mapped          #
        # orphaned - Mapped reads orphaned by trimming #
        ################################################
        fq_1p_out = dir_ / 'reads_1P.fastq'
        fq_2p_out = dir_ / 'reads_2P.fastq'
        fq_single_out = dir_ / 'reads_single.fastq'
        fq_orphaned_out = dir_ / 'reads_orphaned.fastq'

        # Extract reads
        command = Command(' '.join([
            'module load samtools/1.17;',
            f'samtools fastq -N -1 {fq_1p_out} -2 {fq_2p_out} -s {fq_single_out} -0 {fq_orphaned_out} \
            {Path(input.BAM).absolute()};'
        ]))
        command.run(dir_)
        if not command.returncode == 0:
            raise RuntimeError(f'Error extracting FASTQ files: {command.stderr}')
        with open(output.INFORMS, 'w') as handle:
            json.dump({
                '_name': 'samtools fastq',
                '_name_full': 'samtools fastq 1.17;',
                '_version': '1.17',
                '_command': command.command
            }, handle, indent=2)

        # Split singleton reads
        command = Command(' '.join([
            'module load java/18.0.1.1; module load bbtools/39.15;'
            f'demuxbyname.sh in={fq_single_out} out=reads_%U.fastq delimiter="/" prefixmode=f'
        ]))
        command.run(dir_)
        if not command.returncode == 0:
            raise RuntimeError(f'Error splitting singletons: {command.stderr}')

        fq_1u_out = dir_ / 'reads_1U.fastq'
        fq_2u_out = dir_ / 'reads_2U.fastq'

        # Create empty files if they do not exist
        for path_fq in (fq_1u_out, fq_2u_out):
            if not path_fq.exists():
                path_fq.touch()

        # Concatenate orphaned reads to forward unpaired reads
        fq_1u_out_merged = dir_ / 'reads-merged_1U.fastq'
        command = Command(f'cat {fq_1u_out} {fq_orphaned_out} > {fq_1u_out_merged}')
        command.run(dir_)
        if not command.returncode == 0:
            raise RuntimeError(f'Error merging FASTQ files: {command.returncode}')

        # Save output dictionary
        fq_dict_out = FastqInput('illumina', [ToolIOFile(fq_1p_out), ToolIOFile(fq_2p_out)],
            se_fwd=[ToolIOFile(fq_1u_out_merged)], se_rev=[ToolIOFile(fq_2u_out)])
        snakemakeutils.dump_object(fq_dict_out.to_fq_dict(), Path(output.FASTQ))

rule extract_fq_illumina_merge:
    """
    Merges the extracted SE FASTQ files. 
    """
    input:
        FASTQ = [str(rules.extract_fq_illumina.output.FASTQ).format(seq_id=seq_id) for seq_id in config['downsampling'].keys()]
    output:
        FQ_dict = 'merged/illumina/fq_dict.io'
    params:
        dir_ = 'merged/illumina'
    run:
        # Create working directory
        dir_ = Path(params.dir_).absolute()
        dir_.mkdir(parents=True, exist_ok=True)

        # Parse input FASTQ dictionaries
        fq_dicts = [FastqInput.from_fq_dict(Path(p), 'illumina') for p in input.FASTQ]

        output_files = {
            'reads_1P.fastq.gz': [fq_dict.pe[0] for fq_dict in fq_dicts],
            'reads_2P.fastq.gz': [fq_dict.pe[1] for fq_dict in fq_dicts],
            'reads_1U.fastq.gz': [fq_dict.se_fwd[0] for fq_dict in fq_dicts],
            'reads_2U.fastq.gz': [fq_dict.se_rev[0] for fq_dict in fq_dicts],
        }
        for path_out, io_objs in output_files.items():
            command = Command(
                f"cat {' '.join(str(io.path.absolute()) for io in io_objs)} | gzip > {dir_ / path_out};")
            command.run(dir_)
            if not command.returncode == 0:
                raise RuntimeError(f'Error merged FASTQ files: {command.returncode}')
        fq_dict = FastqInput(
            'illumina',
            pe=[ToolIOFile(dir_ / 'reads_1P.fastq.gz'), ToolIOFile(dir_ / 'reads_2P.fastq.gz')],
            se_fwd=[ToolIOFile(dir_ / 'reads_1U.fastq.gz')],
            se_rev=[ToolIOFile(dir_ / 'reads_2U.fastq.gz')],
            is_trimmed=True, is_pe=True)
        snakemakeutils.dump_object(fq_dict.to_fq_dict(), Path(output.FQ_dict))

rule collect_informs:
    """
    Collects the informs for the workflow (limited to the first seq_id).
    """
    input:
        INFORMS_clipping = rules.samtools_ampliconclip.output.INFORMS if config.get('bed_primers') is not None else [],
        INFORMS_create_bam = str(rules.create_bam_for_seq_id.output.INFORMS).format(seq_id=next(iter(config['downsampling'].keys()))),
        INFORMS_fastq = str('by_seq_id/{seq_id}/fastq/{input_type}/informs_{seq_id}.json').format(
            input_type=config['input']['input_type'], seq_id=next(iter(config['downsampling'].keys())))
    output:
        JSON = 'informs_all.json'
    run:
        data_out = []
        for input_str in (input.INFORMS_clipping, input.INFORMS_create_bam, input.INFORMS_fastq):
            if len(str(input_str)) == 0:
                continue
            with Path(input_str).open() as handle:
                data_out.append(json.load(handle))
        with open(output.JSON, 'w') as handle:
            json.dump(data_out, handle, indent=2)
