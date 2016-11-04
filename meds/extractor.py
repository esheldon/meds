"""
MEDSExtractor
    A class to extract a subset of the objects in a MEDS file
    and write to a new file

MEDSCatalogExtractor
    A class to extract the catalog and metadata from a MEDS
    and write to a new file.
"""
from __future__ import print_function
import os
import fitsio
import numpy

def extract_range(meds_file, start, end, sub_file):
    """
    Extract a subset of objects and write a new meds file.

    If you want this as a temporary file, which will be cleaned 
    when you are done with it, use a MEDSExtractor object with
    cleanup=True
    """

    extractor=MEDSExtractor(meds_file, start, end, sub_file)

def extract_catalog(meds_file, sub_file):
    """
    Extract the catalog and meta data from a MEDS file and write
    to a new file

    If you want this as a temporary file, which will be cleaned when you are
    done with it, use a MEDSCatalogExtractor object with cleanup=True
    """

    extractor=MEDSCatalogExtractor(meds_file, sub_file)


class MEDSExtractor(object):
    """
    Class to extract a subset of objects and write a new meds file.

    Optionally clean up the new file when the object is destroyed.
    """
    def __init__(self, meds_file, start, end, sub_file,
                 cleanup=False, copy_all=False):
        self.meds_file=meds_file
        self.start=start
        self.end=end
        self.sub_file=sub_file
        self.cleanup=cleanup
        self.copy_all=copy_all

        # keep track of these in case copy_all is sent, so we
        # know what extra extensions to copy
        self.default_ext=[
            'object_data',
            'image_info',
            'metadata',
            'image_cutouts',
            'weight_cutouts',
            'seg_cutouts',
            'bmask_cutouts',
            'psf',
        ]
        self._check_inputs()
        self._extract()

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        if self.cleanup:
            if os.path.exists(self.sub_file):
                print('removing sub file:',self.sub_file)
                os.remove(self.sub_file)

    def _extract(self):
        
        with fitsio.FITS(self.meds_file) as infits:
            print('opening sub file:',self.sub_file)
            with fitsio.FITS(self.sub_file,'rw',clobber=True) as outfits:

                #
                # subset of object data table
                #
                obj_data = infits['object_data'][self.start:self.end+1]

                cstart, cend = self._get_row_range(obj_data)
                if cstart != -1:
                    # adjust to new start.  If cstart==-1 will all be -9999
                    obj_data['start_row'] -= cstart

                # psfs can have different dimensions
                if 'psf' in infits:
                    psf_cstart, psf_cend = self._get_psf_row_range(obj_data)
                    if psf_cstart != -1:
                        # adjust to new start.  If cstart==-1 will all be -9999
                        obj_data['psf_start_row'] -= psf_cstart

                outfits.write(obj_data, extname='object_data')

                #
                # copy all metadata and image info
                #
                iinfo=infits['image_info'][:]
                outfits.write(iinfo, extname='image_info')

                meta=infits['metadata'][:]
                outfits.write(meta, extname='metadata')

                #
                # extract cutouts for the requested objects
                #
                if cstart == -1:
                    self._write_dummy(outfits)
                else:
                    image_cutouts=infits['image_cutouts'][cstart:cend]
                    outfits.write(image_cutouts, extname='image_cutouts')
                    del image_cutouts

                    weight_cutouts=infits['weight_cutouts'][cstart:cend]
                    outfits.write(weight_cutouts, extname='weight_cutouts')
                    del weight_cutouts

                    seg_cutouts=infits['seg_cutouts'][cstart:cend]
                    outfits.write(seg_cutouts, extname='seg_cutouts')
                    del seg_cutouts

                    if 'bmask_cutouts' in infits:
                        bmask_cutouts=infits['bmask_cutouts'][cstart:cend]
                        outfits.write(bmask_cutouts, extname='bmask_cutouts')
                        del bmask_cutouts

                    if 'psf' in infits:
                        psfs=infits['psf'][psf_cstart:psf_cend]
                        outfits.write(psfs, extname='psf')
                        del psfs

                if self.copy_all:
                    for hdu in infits:
                        if hdu.has_data():
                            extname=hdu.get_extname()
                            if extname not in self.default_ext:
                                h=hdu.read_header()
                                data=hdu.read(header=True)
                                outfits.write(data, header=h, extname=extname)

    def _write_dummy(self, outfits):
        print('no objects with cutouts, writing dummy data')
        dummy=numpy.zeros(2, dtype='f4') + -9999
        outfits.write(dummy, extname='image_cutouts')
        dummy=numpy.zeros(2, dtype='f4')
        outfits.write(dummy, extname='weight_cutouts')
        dummy=numpy.zeros(2, dtype='i4') + -9999
        outfits.write(dummy, extname='seg_cutouts')
        dummy=numpy.zeros(2, dtype='i4') + -9999
        outfits.write(dummy, extname='bmask_cutouts')

    def _get_row_range(self, data):
        """
        get pixel range for this subset
        """
        w,=numpy.where( data['ncutout'] > 0)
        if w.size==0:
            return -1, -1


        ifirst = w[0]
        ilast  = w[-1]

        cstart    = data['start_row'][ifirst,0]

        ncutout_last = data['ncutout'][ilast]
        npix_last = data['box_size'][ilast]**2
        cend     = data['start_row'][ilast,ncutout_last-1] + npix_last

        return cstart, cend

    def _get_psf_row_range(self, data):
        """
        get pixel range for this subset
        """
        w,=numpy.where( data['ncutout'] > 0)
        if w.size==0:
            return -1, -1


        ifirst = w[0]
        ilast  = w[-1]

        cstart    = data['psf_start_row'][ifirst,0]

        ncutout_last = data['ncutout'][ilast]
        npix_last = data['psf_box_size'][ilast]**2
        cend     = data['psf_start_row'][ilast,ncutout_last-1] + npix_last

        return cstart, cend



    def _check_inputs(self):
        if self.meds_file==self.sub_file:
            raise ValueError("output file name equals input")

        if self.start > self.end:
            raise ValueError("found start > end: %d %d" % (self.start,self.end) )


class MEDSCatalogExtractor(object):
    """
    Class to extract the catalog and meta data and write a new file
    """
    def __init__(self, meds_file, new_file, cleanup=False):
        self.meds_file=meds_file
        self.new_file=new_file
        self.cleanup=cleanup

        self._check_inputs()
        self._extract()

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    def __del__(self):
        self.close()

    def close(self):
        if self.cleanup:
            if os.path.exists(self.new_file):
                print('removing cat only file:',self.new_file)
                os.remove(self.new_file)

    def _extract(self):
        
        with fitsio.FITS(self.meds_file) as infits:
            print('opening cat only file:',self.new_file)
            with fitsio.FITS(self.new_file,'rw',clobber=True) as outfits:

                #
                # object data table
                #
                obj_data = infits['object_data'][:]

                outfits.write(obj_data, extname='object_data')

                #
                # copy all metadata and image info
                #
                iinfo=infits['image_info'][:]
                outfits.write(iinfo, extname='image_info')

                meta=infits['metadata'][:]
                outfits.write(meta, extname='metadata')

    def _check_inputs(self):
        if self.meds_file==self.new_file:
            raise ValueError("output file name equals input")

