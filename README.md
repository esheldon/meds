meds
====

Python and C libraries to work with Multi Epoch Data Structures

Description of the file format here
    https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/Multi_Epoch_Data_Structure

python library
-----------------------------

This is a pure python library.  For reading the FITS files, the fitsio python
package is used.  fitsio https://github.com/esheldon/fitsio

MEDS python Docs are here
    https://github.com/esheldon/meds/blob/master/meds/meds.py

installing the python library
#############################

    in the usual place
        python setup.py install

    at a different prefix
        python setup.py install --prefix=/some/path

    make sure it is on your PYTHONPATH


C library
------------------------

This is a pure C library for working with MEDS.  Docs here
    https://github.com/esheldon/meds/blob/master/src/meds.h

installing the C library
#############################

    in the usual place
        make install

    in a different prefix
        make install prefix=/some/path

linking to the library
#######################

INclude "meds.h" and use the following to link against the library.  Make sure
to get the order correct

    CC  ... -lmeds -lcfitsio -lm ...


