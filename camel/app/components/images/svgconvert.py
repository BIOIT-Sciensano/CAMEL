from pathlib import Path

from camel.app.command.command import Command


class SVGConvert:
    """
    Converts SVG images to other formats.
    """

    @staticmethod
    def convert_svg(svg_file, output, size=None):
        """
        Converts a SVG image. The output format is determined based on the extension of the output file.
        :param svg_file: SVG input file
        :param output: Output file
        :param size: Image size
        :return: None
        """
        convert_command = Command(f'convert {svg_file} {output}')
        if size:
            convert_command.command += f' -resize {size}'
        convert_command.run(Path('.'))
