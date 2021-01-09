import os
from pathlib import Path

DIR = Path(os.path.abspath(__file__)).parent
DEFAULT_XML = DIR.joinpath('xml_templates/default_dock.xml')
DEFAULT_BATCH = DIR.joinpath('batch_templates/default_batch.sbatch')
DEFAULT_OPTIONS = DIR.joinpath('options_templates/default_options.txt')
