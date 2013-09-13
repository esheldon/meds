import os
import fitsio

def extract_range(meds_file, start, end, sub_file):
    """
    Extract a subset of objects and write a new meds file.

    If you want this as a temporary file, which will be cleaned 
    when you are done with it, use a MEDSExtractor object with
    cleanup=True
    """

    extractor=MEDSExtractor(meds_file, start, end, sub_file)

class MEDSExtractor(object):
    """
    Class to extract a subset of objects and write a new meds file.

    Optionally clean up the new file when the object is destroyed.
    """
    def __init__(self, meds_file, start, end, sub_file, cleanup=False):
        self.meds_file=meds_file
        self.start=start
        self.end=end
        self.sub_file=sub_file
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
            if os.path.exists(self.sub_file):
                print 'removing sub file:',self.sub_file
                os.remove(self.sub_file)

    def _extract(self):
        
        with fitsio.FITS(self.meds_file) as infits:
            print 'opening sub file:',self.sub_file
            with fitsio.FITS(self.sub_file,'rw',clobber=True) as outfits:

                #
                # subset of object data table
                #
                obj_data = infits['object_data'][self.start:self.end+1]

                cstart, cend = self._get_row_range(obj_data)

                # adjust to new start
                obj_data['start_row'] -= cstart

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
                image_cutouts=infits['image_cutouts'][cstart:cend]
                outfits.write(image_cutouts, extname='image_cutouts')
                del image_cutouts

                weight_cutouts=infits['weight_cutouts'][cstart:cend]
                outfits.write(weight_cutouts, extname='weight_cutouts')
                del weight_cutouts

                seg_cutouts=infits['seg_cutouts'][cstart:cend]
                outfits.write(seg_cutouts, extname='seg_cutouts')
                del seg_cutouts

    def _get_row_range(self, data):
        """
        get pixel range for this subset
        """
        cstart   = data['start_row'][0,0]

        ncutout  = data['ncutout'][-1]
        npix     = data['box_size'][-1]**2 * ncutout
        cend     = data['start_row'][-1,ncutout-1] + npix

        return cstart, cend



    def _check_inputs(self):
        if self.meds_file==self.sub_file:
            raise ValueError("output file name equals input")

        if self.start > self.end:
            raise ValueError("found start > end: %d %d" % (start,end) )

