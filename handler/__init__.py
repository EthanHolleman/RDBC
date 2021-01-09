import os
from pathlib import Path
import logging

DIR = Path(os.path.abspath(__file__)).parent
DEFAULT_XML = DIR.joinpath('xml_templates/default_dock.xml')
DEFAULT_BATCH = DIR.joinpath('batch_templates/default_batch.sbatch')
DEFAULT_OPTIONS = DIR.joinpath('options_templates/default_options.txt')


def make_logger(name):
    '''Creates a logger object specific to the job.
    Logs are then saved in the meta directory of the job.

    Returns:
        Logger: Python logger object for this job.
    '''
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '%(levelname)s\t%(asctime)s\t%(name)s\t%(message)s')
    file_handler = logging.FileHandler('RDBC.log')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger
