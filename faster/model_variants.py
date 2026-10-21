"""Human-readable labels for structural model IDs."""


# The four-character prefix is the stable first segment of a model hash.  Keep
# this mapping here so plots and result collection label model variants alike.
MODEL_ID_TO_VARIANT = {
    'PIRL': 'A',         # default
    '6FY6': 'B',         # no thermodynamics (fitted)
    'KX3J': 'C',         # no Hydro (mit thermodynamics)
    'XIYF': 'old',
    'FA75': 'A',
    'EM33': 'B',
    '3NWE': 'B',         # equivalent to B
    'DXOO': 'C',

    # model version 0.2
    'ZCJ4': 'A',
    'LDHK': 'B',
    '6V72': 'E', # carbon use efficiency constant at 0.5 for all pathways.
    'O26L': 'F', # like E CUE = 0.5 for all AND all _v_max = 1.0
}


def model_variant_from_id(model_id):
    """Return the named variant for a full model ID, or ``None`` if unknown."""
    for model_prefix, variant in MODEL_ID_TO_VARIANT.items():
        if model_id.startswith('model-' + model_prefix):
            return variant
    return '?'
