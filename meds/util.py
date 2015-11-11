from __future__ import print_function
import os
import numpy
import tempfile

DEFVAL = -9999
IMAGE_INFO_TYPES = ['image','weight','seg','bmask','bkg']

def get_meds_output_struct(nobj, ncutout_max, extra_fields=None):
    """
    get the object_data structure, putting in default
    values

    parameters
    ----------
    nobj: int
        Number of entries in the array
    ncutout_max: int
        Number of cutouts to reserve in fields used for each cutout
    extra_fields: numpy descriptor
        A numpy type descriptor to add to the default set
    """
    dtype=get_meds_output_dtype(ncutout_max, extra_fields=extra_fields)
    data = numpy.zeros(nobj, dtype=dtype)

    noset=['ncutout']
    for d in dtype:
        name=d[0]
        if name not in noset:
            data[name] = DEFVAL

    return data

def get_meds_input_struct(nobj, extra_fields=None):
    """
    get the minimal object_data structure for input to the
    MEDSMaker, putting in default values

    parameters
    ----------
    nobj: int
        number of objects (length of array0
    extra_fields: numpy descriptor, optional
        optional extra fields to add 
    """
    dtype=get_meds_input_dtype(extra_fields=extra_fields)
    data = numpy.zeros(nobj, dtype=dtype)

    noset=['ncutout']
    for d in dtype:
        name=d[0]
        if name not in noset:
            data[name] = DEFVAL

    return data


def get_meds_input_dtype(extra_fields=None):
    """
    get the minimal set of entries for input to the MEDSMaker

    parameters
    ----------
    extra_fields: numpy descriptor
        A numpy type descriptor to add to the default set
    """
    dtype = [
        ('id', 'i8'),
        ('box_size', 'i8'),
        ('ra','f8'),
        ('dec','f8'),
    ]

    if extra_fields is not None:
        dtype += extra_fields

    return dtype



def get_meds_output_dtype(ncutout_max, extra_fields=None):
    """
    Full dtype for the object_data structure and file extension

    parameters
    ----------
    ncutout_max: int
        Number of cutouts to reserve in fields used for each cutout
    extra_fields: numpy descriptor
        A numpy type descriptor to add to the default set
    """

    dtype = [
        ('id', 'i8'),
        ('box_size', 'i8'),
        ('ra','f8'),
        ('dec','f8'),
        ('ncutout', 'i8'),
        ('file_id', 'i8', (ncutout_max,)),
        ('start_row', 'i8', (ncutout_max,)),
        ('orig_row', 'f8', (ncutout_max,)),
        ('orig_col', 'f8', (ncutout_max,)),
        ('orig_start_row', 'i8', (ncutout_max,)),
        ('orig_start_col', 'i8', (ncutout_max,)),
        ('cutout_row', 'f8', (ncutout_max,)),
        ('cutout_col', 'f8', (ncutout_max,)),
        ('dudrow', 'f8', (ncutout_max,)),
        ('dudcol', 'f8', (ncutout_max,)),
        ('dvdrow', 'f8', (ncutout_max,)),
        ('dvdcol', 'f8', (ncutout_max,)),
    ]

    if extra_fields is not None:
        dtype += extra_fields

    return dtype


def get_image_info_struct(nimage, path_len, wcs_len=None):
    """
    get the image info structure

    Set default scale to 1.0. The other fields are 0 for
    numbers, or blank for strings

    parameters
    ----------
    nimage: int
        number of images in array
    path_len: int
        length of path strings
    wcs_len: int, optional
        length of wcs strings. If not sent, wcs will not
        be present in the array
    """
    dt = get_image_info_dtype(path_len, wcs_len=wcs_len)

    data = numpy.zeros(nimage, dtype=dt)

    data['scale'] = 1.0

    return data

def get_image_info_dtype(path_len, wcs_len=None):
    """
    get the image_info dtype for the specified path string
    length and wcs string length

    parameters
    ----------
    path_len: int
        length of path strings
    wcs_len: int, optional
        length of wcs strings. If not sent, wcs will not
        be present in data type
    """

    path_fmt = 'S%d' % path_len

    dt=[]
    for ctype in IMAGE_INFO_TYPES:
        path_name = '%s_path' % ctype
        ext_name  = '%s_ext' % ctype

        dt += [
            (path_name, path_fmt),
            (ext_name,'i2'),
        ]

    dt += [
        ('image_id', 'i8'),
        ('image_flags', 'i8'),
        ('magzp', 'f4'),
        ('scale', 'f4'),
        ('position_offset','f8'),
    ]
    if wcs_len is not None:
        wcs_fmt = 'S%d' % wcs_len
        dt += [
            ('wcs',wcs_fmt),
        ]

    return dt



def make_wcs_positions(row, col, offset, inverse=False):
    """
    make a structure holding both the original wcs and zero-offset positions.
    This is only meant for converting between 1-offset and 0-offset

    wcs positions are called 'wcs_row','wcs_col' and zero offset are called
    'zrow','zcol'

    parameters
    ----------
    row: array
        rows in the image, wcs coords if inverse=False
    col: array
        columns in the image, wcs coords if inverse=False
    offset: float
        offset to subtract from the input positions
    inverse: bool, optional
        Set to True if the input are zero based
    """
    n=row.size
    dt=[('wcs_row','f8'),
        ('wcs_col','f8'),
        ('zrow','f8'),
        ('zcol','f8')]

    data=numpy.zeros(row.size, dtype=dt)

    if inverse:
        data['zrow'] = row
        data['zcol'] = col

        data['wcs_row'] = row + offset
        data['wcs_col'] = col + offset
    else:
        data['wcs_row'] = row
        data['wcs_col'] = col

        data['zrow'] = row - offset
        data['zcol'] = col - offset

    return data


