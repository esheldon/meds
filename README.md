meds
====

A Python library to create and read Multi Epoch Data Structures. A minimal C
library is also provided for reading MEDS files

## Documentation

Full instructions for installing and using the python and C libraries are
here: https://github.com/esheldon/meds/wiki

### A few examples for reading files
```python

import meds

# create a MEDS object for the given MEDS file
m=meds.MEDS(filename)

# read a cutout for object 35, cutout index 5
index=35
cutout_index=5
image=m.get_cutout(index, cutout_index)

# read the second cutout for this object
cutout_index=1
im = m.get_cutout(index, cutout_index)

# get other image types
seg  = m.get_cutout(index, cutout_index, type=’seg’)
wt   = m.get_cutout(index, cutout_index, type=’weight’)
mask = m.get_cutout(index, cutout_index, type=’bmask’)

# get a python list of all cutouts for this object
imlist  = m.get_cutout_list(index)
seglist = m.get_cutout_list(index,type=’seg’)

# The contents of the object data table is loaded when the MEDS object is
# created, and are accessible by name.

# number of cutouts
ncutout=m[’ncutout’][index]
for i in xrange(ncutout):
    imlist=m.get_cutout_list(i)
    # process the images

# get the jacobian of the WCS transformation
# as a dict
j=m.get_jacobian(index, cutout_index)

# as a numpy matrix
j=m.get_jacobian_matrix(index, cutout_index)

# list for all cutouts
jlist=m.get_jacobian_list(index)

# get the "ubserseg" weight map
wt=m.get_cweight_cutout_nearest(index, cutout_index)
```
