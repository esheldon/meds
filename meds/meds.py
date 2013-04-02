"""
Defines the MEDS class to work with MEDS (Multi Epoch Data Structures)

See docs for the MEDS class for more info

Author: Erin Sheldon, BNL
"""

class MEDS(object):
    """
    Class to work with MEDS (Multi Epoch Data Structures)

    For details of the structure, see
    https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/Multi_Epoch_Data_Structure

    One can extract cutouts using get_cutout() and get_mosaic() and
    get_cutout_list()

    One can access all fields from the catalog using [field_name] notation. The
    number of entries is in the .size attribute. Note the actual fields in the
    catalog may change over time.  You can use get_cat() to get the full
    catalog as a recarray.

    The first cutout for an object is always from the coadd.

    public methods
    --------------
    get_cutout(iobj, icutout)
        Get a single cutout for the indicated entry
    get_mosaic(iobj)
        Get a mosaic of all cutouts associated with this coadd object
    get_cutout_list(iobj)
        Get an image list with all cutouts associated with this coadd object
    get_source_path(iobj, icutout)
        Get the source filename associated with the indicated cutout
    get_sky_path(iobj, icutout)
        Get the source sky filename associated with the indicated cutout
    get_source_info(iobj, icutout)
        Get all info about the source images
    get_cat()
        Get the catalog; extension 1
    get_image_info()
        Get the entire image info structure
    get_meta()
        Get all the metadata

    examples
    --------
    from deswl import meds
    m=meds.MEDS(filename)

    # number of coadd objects
    num=m.size

    # number of cutouts for object 35
    m['ncutout'][35]

    # get cutout 3 for object 35
    im=m.get_cutout(35,3)

    # get all the cutouts for object 35 as a single image
    mosaic=m.get_mosaic(35)

    # get all the cutouts for object 35 as a list of images
    im=m.get_cutout_list(35)

    # get a cutout for the weight map
    wt=m.get_cutout(35,3,type='weight')

    # get a cutout for the segmentation map
    seg=m.get_cutout(35,3,type='seg')

    # get the source filename for cutout 3 for object 35
    fname=m.get_source_path(35,3)

    # you can access any of the columns in the
    # catalog (stored as a recarray) directly

    # e.g. get the center in the cutout for use in image processing
    row = m['row_cutout'][35]
    col = m['col_cutout'][35]

    # source filename
    fname = m.get_source_path(35,3)

    # or you can just get the catalog to work with
    cat=m.get_cat()
    info=m.get_image_info()


    Fields in main catalog
    -----------------------

     id                 i4       id from coadd catalog
     ncutout            i4       number of cutouts for this object
     box_size           i4       box size for each cutout
     file_id            i4[NMAX] zero-offset id into the file names in the 
                                 second extension
     start_row          i4[NMAX] zero-offset, points to start of each cutout.
     orig_row           f8[NMAX] zero-offset position in original image
     orig_col           f8[NMAX] zero-offset position in original image
     orig_start_row     i4[NMAX] zero-offset start corner in original image
     orig_start_col     i4[NMAX] zero-offset start corner in original image
     cutout_row         f8[NMAX] zero-offset position in cutout imag
     cutout_col         f8[NMAX] zero-offset position in cutout image
     dudrow             f8[NMAX] jacobian of transformation 
                                 row,col->ra,dec tangent plane (u,v)
     dudcol             f8[NMAX]
     dvdrow             f8[NMAX]
     dvdcol             f8[NMAX]


    requirements
    ------------
    numpy
    fitsio https://github.com/esheldon/fitsio
    """
    def __init__(self, filename):
        import fitsio
        self._filename=filename
        
        self._fits=fitsio.FITS(filename)

        self._cat=self._fits["object_data"][:]
        self._image_info=self._fits["image_info"][:]
        self._meta=self._fits["metadata"][:]

    def get_cutout(self, iobj, icutout, type='image'):
        """
        Get a single cutout for the indicated entry

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of the cutout for this object.
        type: string, optional
            Cutout type. Default is 'image'.  Allowed
            values are 'image','weight'

        returns
        -------
        The cutout image
        """
        self._check_indices(iobj, icutout=icutout)

        box_size=self._cat['box_size'][iobj]
        start_row = self._cat['start_row'][iobj,icutout]
        row_end = start_row + box_size*box_size

        extname=self._get_extension_name(type)

        imflat = self._fits[extname][start_row:row_end]
        im = imflat.reshape(box_size,box_size)
        return im

    def get_mosaic(self, iobj, type='image'):
        """
        Get a mosaic of all cutouts associated with this coadd object

        parameters
        ----------
        iobj:
            Index of the object
        type: string, optional
            Cutout type. Default is 'image'.  Allowed
            values are 'image','weight'

        returns
        -------
        An image holding all cutouts
        """

        self._check_indices(iobj)

        ncutout=self._cat['ncutout'][iobj]
        box_size=self._cat['box_size'][iobj]

        start_row = self._cat['start_row'][iobj,0]
        row_end = start_row + box_size*box_size*ncutout

        extname=self._get_extension_name(type)

        mflat=self._fits[extname][start_row:row_end]
        mosaic=mflat.reshape(ncutout*box_size, box_size)

        return mosaic

    def get_cutout_list(self, iobj, type='image'):
        """
        Get an image list with all cutouts associated with this coadd object

        Note each individual cutout is actually a view into a larger
        mosaic of all images.

        parameters
        ----------
        iobj:
            Index of the object
        type: string, optional
            Cutout type. Default is 'image'.  Allowed
            values are 'image','weight'

        returns
        -------
        A list of images hold all cutouts.
        """

        mosaic=self.get_mosaic(iobj,type=type)
        ncutout=self._cat['ncutout'][iobj]
        box_size=self._cat['box_size'][iobj]
        return self._split_mosaic(mosaic, box_size, ncutout)

    def get_source_info(self, iobj, icutout):
        """
        Get the full source file information for the indicated cutout.

        Includes SE image and sky image

        parameters
        ----------
        iobj: 
            Index of the object
        """
        self._check_indices(iobj, icutout=icutout)
        ifile=self._cat['file_id'][iobj,icutout]
        return self._image_info[ifile]

    def get_source_path(self, iobj, icutout):
        """
        Get the source filename associated with the indicated cutout

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of the cutout for this object.

        returns
        -------
        The filename
        """

        info=self.get_source_info(iobj, icutout)
        return info['image_path']

    def get_sky_path(self, iobj, icutout):
        """
        Get the source filename associated with the indicated cutout

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of the cutout for this object.

        returns
        -------
        The filename
        """

        info=self.get_source_info(iobj, icutout)
        return info['sky_path']


    def get_cat(self):
        """
        Get the catalog
        """
        return self._cat

    def get_image_info(self):
        """
        Get all image information
        """
        return self._image_info
    
    def get_meta(self):
        """
        Get all the metadata
        """
        return self._meta

    def _get_extension_name(self, type):
        if type=='image':
            return "image_cutouts"
        elif type=="weight":
            return "weight_cutouts"
        elif type=="seg":
            return "seg_cutouts"
        else:
            raise ValueError("bad cutout type '%s'" % type)

    def _split_mosaic(self, mosaic, box_size, ncutout):
        imlist=[]
        for i in xrange(ncutout):
            r1=i*box_size
            r2=(i+1)*box_size
            imlist.append( mosaic[r1:r2, :] )

        return imlist


    def _check_indices(self, iobj, icutout=None):
        if iobj >= self._cat.size:
            raise ValueError("object index should be within "
                             "[0,%s)" % self._cat.size)

        ncutout=self._cat['ncutout'][iobj]
        if ncutout==0:
            raise ValueError("object %s has no cutouts" % iobj)

        if icutout is not None:
            if icutout >= ncutout:
                raise ValueError("requested cutout index %s for "
                                 "object %s should be in bounds "
                                 "[0,%s)" % (icutout,iobj,ncutout))

    def __repr__(self):
        return repr(self._fits[1])
    def __getitem__(self, item):
        return self._cat[item]

    @property
    def size(self):
        return self._cat.size
