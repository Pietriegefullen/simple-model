import os

ROOT_DIRECTORY = os.path.dirname(__file__)
#ROOT_DIRECTORY = os.path.join('D:', os.sep, 'OneDrive - Universität Hamburg','simple model')

simple_model_dir, _ = os.path.split(ROOT_DIRECTORY)
LOG_DIRECTORY = os.path.join(simple_model_dir, 'optimization_log')