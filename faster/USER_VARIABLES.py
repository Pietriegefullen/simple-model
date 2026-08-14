import os

ROOT_DIRECTORY = os.path.dirname(__file__)

simple_model_dir, _ = os.path.split(ROOT_DIRECTORY)
LOG_DIRECTORY = os.path.join(simple_model_dir, '20260814_fit')
DATA_DIRECTORY = os.path.join(ROOT_DIRECTORY, 'data')
