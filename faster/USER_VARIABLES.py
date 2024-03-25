import os

ROOT_DIRECTORY = os.path.dirname(__file__)

simple_model_dir, _ = os.path.split(ROOT_DIRECTORY)
LOG_DIRECTORY = os.path.join(simple_model_dir, 'model_comparison_results')