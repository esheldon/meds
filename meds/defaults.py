default_config = {
    # buffer bounds for first cull of objects not on chip
    # in arcsec
    'bounds_buffer_uv': 16.0,

    # buffer for "coadd" (first entry) image in the meds in pixels
    # objects within the coadd image bounds plus this buffer are 
    # "in" the coadd image
    'coadd_bounds_buffer_rowcol': 1e-3,
    
    # objects within 'coadd_bounds_buffer_rowcol' buffer around 
    # the bounds of the coadd are forced to be at the edge of 
    # the chip if 'force_into_coadd_bounds' is True    
    'force_into_coadd_bounds': True,

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

