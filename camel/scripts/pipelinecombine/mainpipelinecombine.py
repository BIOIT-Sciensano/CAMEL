#!/usr/bin/env python
import ast
import dataclasses
import re
from pathlib import Path
from typing import Any, Union

import click
import humanize
import pandas as pd

from camel.app.cli import cliutils
from camel.app.loggers import initialize_logging, logger
from camel.app.scriptutils.basescript.basescript import BaseScript
from camel.app.scriptutils.model import BaseInput, BaseOptions, BaseOutput
from camel.version import __VERSION__

GENE_FORMATS = {
    'simple': '{hit[1]}',
    'locus_with_id': '{hit[1]} ({hit[2]}%)',
    'locus_with_id_len': '{hit[1]} (id={hit[2]}%, len={hit[3]})',
}


@dataclasses.dataclass(frozen=True)
class PipelineCombineInput(BaseInput):
    """
    Inputs for the pipeline combine script.
    """
    tsv: list[Path]

@dataclasses.dataclass(frozen=True)
class PipelineCombineOutput(BaseOutput):
    """
    Outputs for the pipeline combine script.
    """
    output: Path = dataclasses.field(metadata={'help': 'TSV output file'})

@dataclasses.dataclass(frozen=True)
class PipelineCombineOpts(BaseOptions):
    """
    Options for the pipeline combine script.
    """
    exclude: str | None = dataclasses.field(default=None, metadata={'help': 'Comma separated list of keys to exclude'})
    include: str | None = dataclasses.field(default=None, metadata={'help': 'Comma separated list of keys to include'})
    gene_format: str | None = dataclasses.field(default=None, metadata={
        'help': 'Format for genes',
        'choices': list(GENE_FORMATS.keys())
    })
    group_genes: str = dataclasses.field(default='cluster', metadata={
        'help': 'Grouping for the detected genes in the output',
        'choices': ['cluster', 'gene', 'allele']
    })


class MainPipelineCombine(BaseScript[PipelineCombineInput, PipelineCombineOutput, PipelineCombineOpts]):
    """
    Main script for the pipeline combine tool.
    """

    def __init__(
            self,
            script_in: PipelineCombineInput,
            script_out: PipelineCombineOutput,
            opts: PipelineCombineOpts
    ) -> None:
        """
        Initializes the main script.
        :param script_in: Script input
        :param script_out: Script output
        :param opts: Script options
        :return: None
        """
        super().__init__(
            name='Pipeline Combine',
            version=f'CAMEL_{__VERSION__}',
            script_in=script_in,
            script_out=script_out,
            script_opts=opts
        )
        self._gene_detection_cols = []

    def _execute(self) -> None:
        """
        Runs this tool.
        """
        # Parse input data
        logger.debug(f'Parsing {len(self._script_in.tsv)} input files')
        records_out = []
        for path_input in self._script_in.tsv:
            isolate_data = self._parse_pipe_out(path_input)
            records_out.append(isolate_data)
        data_out = pd.DataFrame(records_out)
        logger.info(f'Summaries parsed for {len(data_out)} datasets')

        # Check if the pipelines are the same
        if len(data_out['pipeline_name'].unique()) > 1:
            logger.warning('Combining output of multiple pipelines (not recommended!)')

        # Filter columns
        cols_to_keep = data_out.apply(lambda x: MainPipelineCombine._filter_columns(
            x, self._script_opts.exclude, self._script_opts.include, self._script_opts.gene_format is not None), axis=0)
        data_out_filt = data_out.loc[:, cols_to_keep]
        logger.info(f'Keeping {len(data_out_filt.columns)}/{len(data_out.columns)} columns')

        # Re-order columns
        data_out_filt = data_out_filt.reindex(sorted(
            data_out_filt.columns, key=lambda x: MainPipelineCombine._get_sorting_key(x, self._gene_detection_cols)),
            axis=1)

        # Save the output file
        path_out = self._script_out.output
        data_out_filt.to_csv(path_out, sep='\t', index=False)
        logger.info(
            f'Output file created: {path_out} ({humanize.naturalsize(path_out.stat().st_size)})')

    def _parse_pipe_out(self, path_in: Path) -> dict[str, Any]:
        """
        Parses a pipeline summary output file.
        :param path_in: Input path
        :return: Parsed info as a dictionary
        """
        isolate_data = {}
        with path_in.open() as handle:
            for line in handle.readlines():
                parts = line.strip().split('\t')
                key = parts[0]
                value = '\t'.join(parts[1:])

                # Gene detection hits
                m = re.match(r'hits_(.*)', key)
                if m and (self._script_opts.gene_format is not None):
                    # Keep track of gene detection columns for sorting
                    self._gene_detection_cols.append(m.group(1))

                    # Add detected genes
                    for key, value in MainPipelineCombine._format_gene_detection_hits(
                            m.group(1), value, GENE_FORMATS[self._script_opts.gene_format], self._script_opts.group_genes):
                        isolate_data[key] = value

                # LREFinder genes
                if key == 'lrefinder_genes':
                    hits = ast.literal_eval(value)
                    for h in hits:
                        isolate_data[f"lrefinder_genes-{h['Gene']}"] = f"{h['Gene']} ({h['Template identity']})"
                else:
                    isolate_data[key] = value
        return isolate_data

    @staticmethod
    def _filter_columns(column: pd.Index, keys_exclude: Union[str, None], keys_include: Union[str, None],
                        genes_formatted: bool) -> bool:
        """
        Function to check if the target row should be filtered.
        :param column: Columns to check
        :param keys_exclude: Keys that are filtered out (with wildcard matching)
        :param keys_include: Keys that are included (no wildcard matching)
        :param genes_formatted: Boolean to indicate if genes are formatted
        :return: True if column should be kept, False otherwise
        """
        # Remove excluded columns
        if keys_exclude is not None:
            for key_excluded in keys_exclude.split(','):
                # Expand gene detection hits
                if genes_formatted and key_excluded.startswith('hits_'):
                    key = f"{key_excluded.replace('hits_', '')}*"
                else:
                    key = key_excluded
                if MainPipelineCombine.key_matches(key, str(column.name)):
                    return False
            return True

        # Retain included columns
        elif keys_include is not None:
            for key_included in keys_include.split(','):
                # Expand gene detection hits
                if genes_formatted and key_included.startswith('hits_'):
                    key = f"{key_included.replace('hits_', '')}*"
                else:
                    key = key_included
                if MainPipelineCombine.key_matches(key, str(column.name)):
                    return True
            return False

        # If both options are not specified, retain all columns
        return True

    @staticmethod
    def key_matches(pattern: str, key: str) -> bool:
        """
        Checks if the pattern matches the key.
        :param pattern: Pattern
        :param key: Key
        :return: True if matches, False otherwise
        """
        # Regular full length match
        if '*' not in pattern:
            if key == pattern:
                return True

        # Regex match
        if '*' in pattern:
            pattern = pattern.replace('*', '.*')
            m = re.match(pattern, key)
            if m:
                return True

        # Gene detection match
        if 'hits_' in pattern:
            if key == pattern.replace('hits_', ''):
                return True

        return False

    @staticmethod
    def _get_sorting_key(col_name: str, gene_detection_cols: list[str]) -> str:
        """
        Returns the sorting key for the given column.
        :param col_name: Column name
        :param gene_detection_cols: Gene detection columns
        :return: sorting key
        """
        try:
            next(c for c in gene_detection_cols if col_name.startswith(c))
            return col_name
        except StopIteration:
            # Not a gene detection column -> don't change order
            return ''

    @staticmethod
    def _format_gene_detection_hits(db_key: str, value: str, format_str: str, grouping: str) -> list[tuple[str, str]]:
        """
        Formats a gene entry.
        :param db_key: Database key
        :param value: Value
        :param format_str: Format string
        :param grouping: Option for determining the key for grouping detected genes
        :return: Key, formatted hit
        """
        hits = ast.literal_eval(value)
        if len(hits) == 0:
            return []
        output_rows = []
        for hit in hits:
            if grouping == 'cluster':
                key = f'{db_key}_{hit[0]}'
            elif grouping == 'gene':
                gene_name = re.sub(r'_\d+$', '', hit[1])
                key = f'{db_key}_{gene_name}'
            else:
                key = f'{db_key}_{hit[1]}'
            output_rows.append((key, format_str.format(hit=hit)))
        return output_rows

@click.command(name='pipeline_combine', short_help='Combines the output of multiple pipeline runs')
@click.argument("input", type=click.Path(path_type=Path), nargs=-1, required=True)
@cliutils.add_click_options_from_dataclass(PipelineCombineOutput)
@cliutils.add_click_options_from_dataclass(PipelineCombineOpts)
def main(**kwargs) -> None:
    """
    Combines the output of multiple pipeline runs.
    """
    script = MainPipelineCombine(
        script_in=PipelineCombineInput(tsv=[Path(x) for x in kwargs['input']]),
        script_out=PipelineCombineOutput(**cliutils.from_kwargs(PipelineCombineOutput, kwargs)),
        opts=PipelineCombineOpts(**cliutils.from_kwargs(PipelineCombineOpts, kwargs))
    )
    script.run()


if __name__ == '__main__':
    initialize_logging()
    main()
