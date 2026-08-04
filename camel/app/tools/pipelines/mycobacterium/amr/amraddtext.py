from camelcore.app.io.tooliofile import ToolIOFile
from PIL import Image, ImageDraw, ImageFont

from camel.app.core.errors import InvalidToolInputError
from camel.app.core.tool import Tool
from camel.app.loggers import logger
from camel.resources import FONT_SANS, FONT_SANS_BOLD


class AMRAddText(Tool):
    """
    Adds text to the visualization of the resistance characterization.
    """

    def __init__(self) -> None:
        """"
        Initializes this tool.
        """
        super().__init__('Mycobacterium: visualization add text', '0.1')

    def _check_input(self) -> None:
        """
        Checks if the provided input is valid.
        :return: None
        """
        if 'PNG' not in self._tool_inputs:
            raise InvalidToolInputError("Input image (PNG) is required")
        if 'VAL_sample' not in self._tool_inputs:
            raise InvalidToolInputError("Sample name input is required")
        if 'coverage' not in self._input_informs:
            raise InvalidToolInputError("Coverage informs are required")
        if 'lineage' not in self._input_informs:
            raise InvalidToolInputError("Lineage informs are required")
        super()._check_input()

    def _execute_tool(self) -> None:
        """
        Executes this tool.
        :return: None
        """
        # Setup image
        img = Image.open(self._tool_inputs['PNG'][0].path)
        img_width, img_height = img.size
        logger.info(f"Image width: {img_width}, image height: {img_height}")
        draw = ImageDraw.Draw(img)

        # Get fonts
        font_large_bold = ImageFont.truetype(str(FONT_SANS_BOLD), 64)
        font_normal = ImageFont.truetype(str(FONT_SANS), 48)

        # Create messages
        messages = [
            (self._tool_inputs['VAL_sample'][0].value, font_large_bold, 24),
            (self._input_informs['lineage']['detected_lineage'].name, font_normal, 8),
            (f"Median coverage: {int(self._input_informs['coverage']['median_depth'])}X", font_normal, 8)
        ]
        x_by_msg, y_by_msg = self.calculate_positions(messages, img_width, img_height)

        # Add to text to image
        for msg, font, _ in messages:
            draw.text((x_by_msg[msg], y_by_msg[msg]), msg, (0, 0, 0), font=font)

        # Save output file
        output_path = self._folder / 'img_edited.png'
        logger.info(f"Saving image to: {output_path}")
        img.save(str(output_path))
        self._tool_outputs['PNG'] = [ToolIOFile(output_path)]

    def calculate_positions(self, messages: list[tuple[str, 'ImageFont', int]], img_width: int, img_height: int) \
            -> tuple[dict[str, int], dict[str, int]]:
        """
        This function calculates the position to align the added text horizontally and vertically.
        :param messages: Messages to add (Message str, font, spacing below)
        :param img_width: Original img width
        :param img_height: Original img height
        :return: x, y coordinates by message
        """
        height_by_msg = {msg: font.getsize(msg)[1] for msg, font, _ in messages}
        spacing_by_msg = {msg: spacing for msg, _, spacing in messages}
        logger.info(f"Height per message: {height_by_msg} - spacing: {spacing_by_msg}")
        total_height = sum(height_by_msg.values()) + sum([spacing for _, _, spacing in messages])

        # Align vertically
        height_y0 = (img_height - total_height) / 2
        y_by_msg = {}
        for message, _, _ in messages:
            height_added_components = sum([height_by_msg[m] for m in y_by_msg.keys()])
            spacing_added_components = sum([spacing_by_msg[m] for m in y_by_msg.keys()])
            y_by_msg[message] = int(height_y0 + height_added_components + spacing_added_components)
        logger.info(f"Calculated y-coordinates: {y_by_msg}")

        # Align horizontally
        x_by_msg = {msg: 0.5 * (img_width - font.getsize(msg)[0]) for msg, font, _ in messages}
        return x_by_msg, y_by_msg
