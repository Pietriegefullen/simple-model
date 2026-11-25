import os

CODE_DIRECTORY = os.path.dirname(__file__)
PROJECT_DIRECTORY, _ = os.path.split(CODE_DIRECTORY)
RESULTS_DIRECTORY = os.path.join(PROJECT_DIRECTORY, 'results')
CONVERTED_CHECKPOINTS_DIRECTORY = os.path.join(PROJECT_DIRECTORY, 'converted_checkpoints')
DATA_DIRECTORY = CODE_DIRECTORY
PLOT_DIRECTORY = os.path.join(PROJECT_DIRECTORY, 'plots')
