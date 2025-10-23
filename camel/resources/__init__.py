from importlib.resources import files
from pathlib import Path

DIR_RESOURCES = Path(__file__).parent.absolute()

# CSS
CSS_STYLE = files('camel').joinpath('resources/reports/style.css')

# FONTS
FONT_SANS = DIR_RESOURCES / 'fonts' / 'FreeSans.ttf'
FONT_SANS_BOLD = DIR_RESOURCES / 'fonts' / 'FreeSansBold.ttf'

# Static
LOGO_SCIENSANO = files('camel').joinpath('resources/reports/static/logo-sciensano.png')

# Citations
DIR_CITATIONS = DIR_RESOURCES / 'citations'
