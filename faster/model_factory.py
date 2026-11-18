# model_versions/registry.py
from importlib import import_module

BUNDLES = {
    "0.3": "model_v0_3.model",
    "0.4": "model_v0_4.model",
    "0.4.1": "model_v0_4.model",
    "0.4.2": "model_v0_4.model",
    "0.4.3": "model_v0_4_3.model", # using C_ATOMS_BIOMASS
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

def build_model(version, pathway_names=None, normalized=True):
    """Build a versioned model.

    The two-argument historical form, ``build_model(pathways, normalized)``,
    continues to select the current version for interactive tools.
    """
    if not isinstance(version, str):
        if pathway_names is not None and not isinstance(pathway_names, bool):
            raise TypeError('Legacy build_model accepts only pathways and normalized.')
        normalized = normalized if pathway_names is None else pathway_names
        pathway_names = version
        version = current_version()
    return implementation(version).build_model(pathway_names, normalized)


def envelope_converted_model_parameters(version, normalized=True):
    """Return the versioned parameter catalogue for converted envelopes."""
    model_class = implementation(version).Model
    try:
        return model_class.envelope_converted_model_parameters(normalized)
    except AttributeError as error:
        raise ValueError(
            f'Model version {version!r} has no converted-envelope parameter catalogue.'
        ) from error


def hypothesis_model_parameters(version, normalized=True):
    """Return the versioned parameter catalogue for the hypothesis fit."""
    model_class = implementation(version).Model
    try:
        return model_class.hypothesis_model_parameters(normalized)
    except AttributeError as error:
        raise ValueError(
            f'Model version {version!r} has no hypothesis parameter catalogue.'
        ) from error

def pathway_by_name(version, name):
    return implementation(version).pathway_by_name(name)

if __name__ == '__main__':
    print(current_version())
