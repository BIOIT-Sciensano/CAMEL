from importlib.resources import files
from pathlib import Path

PATH_JQUERY = Path(str(files('camel').joinpath('resources/javascript/jquery-3.2.1.min.js')))
PATH_CSS = Path(str(files('camel').joinpath('resources/css/style.css')))
