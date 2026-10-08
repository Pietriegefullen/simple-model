"""Human-readable labels for structural model IDs."""


# The four-character prefix is the stable first segment of a model hash.  Keep
# this mapping here so plots and result collection label model variants alike.
MODEL_ID_TO_VARIANT = {
    'PIRL': 'A',         # default
    '6FY6': 'B',         # no thermodynamics
    'KX3J': 'C',         # no Hydro
    'XIYF': 'old',
    'FA75': 'A',
    'EM33': 'B',
    '3NWE': 'B',         # equivalent to B
    'DXOO': 'C',
}


def model_variant_from_id(model_id):
    """Return the named variant for a full model ID, or ``None`` if unknown."""
    for model_prefix, variant in MODEL_ID_TO_VARIANT.items():
        if model_id.startswith('model-' + model_prefix):
            return variant
    return None
