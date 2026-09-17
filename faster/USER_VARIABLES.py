import os

CODE_DIRECTORY = os.path.dirname(__file__)

PROJECT_DIRECTORY, _ = os.path.split(CODE_DIRECTORY)
RESULTS_DIRECTORY = os.path.join(PROJECT_DIRECTORY, 'results')
DATA_DIRECTORY = CODE_DIRECTORY
