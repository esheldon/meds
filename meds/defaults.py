default_config = {
    # buffer bounds for first cull of objects not on chip
    # in arcsec
    'bounds_buffer_uv': 16.0,

    # allowed values in the bitmask image
    'bitmask_allowed': 0,

    # cutout types in addition to 'image'.  Allowed values are
    # ['weight','seg','bmask']
    'cutout_types': [],

    # default output data types for images
    'image_dtype':'f4',
    'weight_dtype':'f4',
    'seg_dtype':'i4',
    'bmask_dtype':'i4',
}

