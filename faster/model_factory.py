# model_versions/registry.py
from importlib import import_module

BUNDLES = {
    "0.3": "model_v0_3.model",
    "0.4": "model_v0_4.model",
}

def current_version():
    return max(BUNDLES.keys())

def implementation(version):
    try:
        return import_module(BUNDLES[str(version)])
    except KeyError as error:
        raise ValueError(
            f"No installed implementation for model version {version!r}"
        ) from error

def build_model(version, pathway_names, normalized):
    return implementation(version).build_model(pathway_names, normalized)

def pathway_by_name(version, name):
    return implementation(version).pathway_by_name(name)

if __name__ == '__main__':
    print(current_version())
