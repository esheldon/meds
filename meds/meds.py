from __future__ import print_function
import numpy
import fitsio

try:
    xrange  # noqa: F821
except Exception:
    xrange = range

try:
    from . import _uberseg
    _have_c_ubserseg = True
except ImportError:
    print("could not load fast ubserseg")
    _have_c_ubserseg = False


class MEDS(object):
    """
    Class to work with MEDS (Multi Epoch Data Structures)

    For details of the data structure, see
    https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/Multi_Epoch_Data_Structure

    One can extract cutouts using `get_cutout()` and `get_mosaic()` and
    `get_cutout_list()`.

    One can access all fields from the catalog using `[field_name]` notation.
    The number of entries is in the `.size` attribute. Note the actual fields
    in the catalog may change over time.  You can use `get_cat()` to get the
    full catalog as a `recarray`.

    The first cutout for an object is always from the coadd.

    Parameters
    ----------
    filename : str
        The path to the MEDS file.

    Attributes
    ----------
    size : int
        The number of objects in the file.

    Methods
    -------
    close()
        Close the underlying FITS file.
    get_cutout(iobj, icutout, type='image')
        Get a single cutout for the indicated entry and image type.
    get_mosaic(iobj, type='image')
        Get a mosaic of all cutouts associated with this coadd object.
    get_cutout_list(iobj, type='image')
        Get an image list with all cutouts associated with this coadd object.
    get_psf(iobj, icutout)
        Get a single psf image for the indicated entry.
    get_psf_list(iobj)
        Get a list of psf images.
    get_cweight_cutout(iobj, icutout, restrict_to_seg=False)
        Composite the weight and seg maps, interpolating seg map from the
        coadd.
    get_cweight_mosaic(iobj)
        Composite the weight and seg maps, interpolating seg map from the
        coadd.
    get_cweight_cutout_list(iobj)
        Composite the weight and seg maps, interpolating seg map from the
        coadd.
    get_uberseg(iobj, icutout, fast=True)
        Get the cweight map and zero out pixels not nearest to central object.
    get_cweight_cutout_nearest
        Alias for get_uberseg.
    get_uberseg_list(iobj, fast=True)
        Composite the weight and seg maps, interpolating seg map from the
        coadd.
    get_cweight_cutout_nearest_list
        Alias for get_uberseg_list.
    get_cseg_cutout(iobj, icutout)
        Interpolate the coadd seg onto the plane of the cutout.
    get_cseg_mosaic(iobj)
        Interpolate the coadd seg onto the planes of the cutouts. Get
        a big mosaic of all.
    get_cseg_cutout_list(iobj)
        Interpolate the coadd seg onto the planes of the cutouts.
        Get a list of all seg cutouts.
    get_cseg_weight(iobj, icutout, use_canonical_cen=False)
        Get the largest circularly masked weight map that does not
        interesect any other objects seg map.
    interpolate_coadd_seg_mosaic(iobj)
        Get a mosaic of interpolated seg maps.
    interpolate_coadd_seg(iobj, icutout)
        Interpolate the coadd segmentation map onto the SE image frame.
    get_source_info(iobj, icutout)
        Get the full source file information for the indicated cutout.
    get_source_path(iobj, icutout)
        Get the source filename associated with the indicated cutout.
    get_cat()
        Get the catalog.
    get_image_info()
        Get the image information.
    get_meta()
        Get the metadata.
    get_jacobian(iobj, icutout)
        Get the jacobian as a dictionary.
    get_jacobian_matrix(iobj, icutout)
        Get the jacobian as a numpy array.
    get_jacobian_list(iobj)
        Get the list of jacobians for all cutouts
        for this object.
    get_number(iobj)
        Get the segmentation map number.
    get_cutout_rowcol(iobj, icutout)
        Get cutout_row, cutout_col for the specified object
        and epoch.

    Notes
    -----
    The underlying catalog (e.g., that obtained with `get_cat`) has atleast the
    following columns:

        id                 i8       id from coadd catalog
        ncutout            i8       number of cutouts for this object
        box_size           i8       box size for each cutout
        file_id            i8[NMAX] zero-offset id into the file names in the
                                    second extension
        start_row          i8[NMAX] zero-offset, points to start of each cutout
        orig_row           f8[NMAX] zero-offset position in original image
        orig_col           f8[NMAX] zero-offset position in original image
        orig_start_row     i8[NMAX] zero-offset start corner in original image
        orig_start_col     i8[NMAX] zero-offset start corner in original image
        cutout_row         f8[NMAX] zero-offset position in cutout imag
        cutout_col         f8[NMAX] zero-offset position in cutout image
        dudrow             f8[NMAX] jacobian of transformation
                                    row,col->ra,dec tangent plane (u,v)
        dudcol             f8[NMAX] ...
        dvdrow             f8[NMAX] ...
        dvdcol             f8[NMAX] ...

    Examples
    --------
    >>> from deswl import meds
    >>> m = meds.MEDS(filename)

    >>> # number of coadd objects
    >>> num = m.size

    >>> # number of cutouts for object 35
    >>> m['ncutout'][35]

    >>> # get cutout 3 for object 35
    >>> im = m.get_cutout(35, 3)

    >>> # get all the cutouts for object 35 as a single image
    >>> mosaic = m.get_mosaic(35)

    >>> # get all the cutouts for object 35 as a list of images
    >>> im = m.get_cutout_list(35)

    >>> # get a cutout for the weight map
    >>> wt = m.get_cutout(35, 3, type='weight')

    >>> # get a cutout for the segmentation map
    >>> seg = m.get_cutout(35, 3, type='seg')

    >>> # get a cutout for the bitmask map
    >>> seg = m.get_cutout(35, 3, type='bmask')

    >>> # get the source filename for cutout 3 for object 35
    >>> fname = m.get_source_path(35, 3)

    >>> # you can access any of the columns in the
    >>> # catalog (stored as a recarray) directly
    >>> # e.g. get the center in the cutout for use in image processing
    >>> row = m['row_cutout'][35]
    >>> col = m['col_cutout'][35]

    >>> # source filename
    >>> fname = m.get_source_path(35, 3)

    >>> # or you can just get the catalog to work with
    >>> cat = m.get_cat()
    >>> info = m.get_image_info()
    >>> meta = m.get_meta()
    """
    def __init__(self, filename):
        self._filename = filename

        self._fits = fitsio.FITS(filename)

        self._cat = self._fits["object_data"][:]
        self._image_info = self._fits["image_info"][:]
        self._meta = self._fits["metadata"][:]

    def close(self):
        self._fits.close()

    def get_cutout(self, iobj, icutout, type='image'):
        """Get a single cutout for the indicated entry and image type.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int
            Index of the cutout for this object.
        type: string, optional
            Cutout type. Default is 'image'. Allowed values are 'image',
            'weight', 'seg', 'bmask', 'ormask', 'noise' or 'psf'.

        Returns
        -------
        cutout : np.array
            The cutout image.
        """

        if type == 'psf':
            return self.get_psf(iobj, icutout)

        self._check_indices(iobj, icutout=icutout)

        box_size = self._cat['box_size'][iobj]
        start_row = self._cat['start_row'][iobj, icutout]
        row_end = start_row + box_size*box_size

        extname = self._get_extension_name(type)

        imflat = self._fits[extname][start_row:row_end]
        im = imflat.reshape(box_size, box_size)
        return im

    def get_mosaic(self, iobj, type='image'):
        """Get a mosaic of all cutouts associated with this coadd object.

        Parameters
        ----------
        iobj : int
            Index of the object.
        type: string, optional
            Cutout type. Default is 'image'. Allowed values are 'image',
            'weight', 'seg', 'bmask', 'ormask', 'noise' or 'psf'.

        Returns
        -------
        mosaic : np.array
            An image holding all of the cutouts.
        """

        self._check_indices(iobj)

        ncutout = self._cat['ncutout'][iobj]
        box_size = self._cat['box_size'][iobj]

        start_row = self._cat['start_row'][iobj, 0]
        row_end = start_row + box_size*box_size*ncutout

        extname = self._get_extension_name(type)

        mflat = self._fits[extname][start_row:row_end]
        mosaic = mflat.reshape(ncutout*box_size, box_size)

        return mosaic

    def get_cutout_list(self, iobj, type='image'):
        """Get an image list with all cutouts associated with this coadd object.

        Note each individual cutout is actually a view into a larger
        mosaic of all images.

        Parameters
        ----------
        iobj : int
            Index of the object.
        type: string, optional
            Cutout type. Default is 'image'. Allowed values are 'image',
            'weight', 'seg', 'bmask', 'ormask', 'noise' or 'psf'.

        Returns
        -------
        image_list : list of np.arrays
            A list of images of each cutout.
        """

        if type == 'psf':
            return self.get_psf_list(iobj)

        mosaic = self.get_mosaic(iobj, type=type)
        return split_mosaic(mosaic)

    def has_psf(self):
        """
        returns True if psfs are in the file
        """
        return 'psf' in self._fits

    def get_psf(self, iobj, icutout):
        """Get a single psf image for the indicated entry.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int
            Index of the cutout for this object.

        Returns
        -------
        psf : np.array
            The PSF as a numpy array.
        """

        if not self.has_psf():
            raise RuntimeError("this MEDS file has no 'psf' extension")

        self._check_indices(iobj, icutout=icutout)

        cat = self._cat
        shape = self._get_psf_shape(iobj, icutout)
        npix = shape[0]*shape[1]

        start_row = self._cat['psf_start_row'][iobj, icutout]
        row_end = start_row + npix

        imflat = self._fits['psf'][start_row:row_end]
        im = imflat.reshape(shape)
        return im

    def _get_psf_shape(self, iobj, icutout):
        cat = self._cat
        if 'psf_row_size' in cat.dtype.names:
            if len(cat['psf_row_size'].shape) > 1:
                shape = (
                    cat['psf_row_size'][iobj,icutout],
                    cat['psf_col_size'][iobj,icutout],
                )
            else:
                shape = (
                    cat['psf_row_size'][iobj],
                    cat['psf_col_size'][iobj],
                )
        else:
            if len(cat['psf_box_size'].shape) > 1:
                box_size = self._cat['psf_box_size'][iobj, icutout]
            else:
                box_size = self._cat['psf_box_size'][iobj]

            shape = (box_size, box_size)

        return shape

    def get_psf_list(self, iobj):
        """Get a list of psf images.

        Parameters
        ----------
        iobj : int
            Index of the object.

        Returns
        -------
        list : list of np.arrays
            A list of the PSFs as numpy arrays.
        """
        ncut = self['ncutout'][iobj]
        return [self.get_psf(iobj, icut) for icut in xrange(ncut)]

    def get_cweight_cutout(self, iobj, icutout, restrict_to_seg=False):
        """Composite the weight and seg maps, interpolating seg map from the
        coadd.

        The weight is set to zero outside the region as defined in the coadd.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int
            Index of cutout.
        restrict_to_seg : bool, optional
            Set weights to zero for any pixel in the image not associated
            with the central object.

        Returns
        -------
        weight : np.array
            The weight map as a numpy array.
        """
        wt = self.get_cutout(iobj, icutout, type='weight')
        coadd_seg = self.get_cutout(iobj, 0, type='seg')
        cwt = self._make_composite_image(iobj, icutout, wt, coadd_seg,
                                         restrict_to_seg=restrict_to_seg)
        return cwt

    def get_cweight_mosaic(self, iobj, restrict_to_seg=False):
        """Composite the weight and seg maps, interpolating seg map from the
        coadd.

        The weight is set to zero outside the region as defined in the coadd.

        Parameters
        ----------
        iobj : int
            Index of the object.
        restrict_to_seg : bool, optional
            Set weights to zero for any pixel in the image not associated
            with the central object.

        Returns
        -------
        mosaic : np.array
            A mosaic of the composite weight maps.
        """
        wtmosaic = self.get_mosaic(iobj, type='weight')
        coadd_seg = self.get_cutout(iobj, 0, type='seg')

        # shares underlying storage
        wlist = split_mosaic(wtmosaic)

        for icutout, wt in enumerate(wlist):
            cwt = self._make_composite_image(iobj, icutout, wt, coadd_seg,
                                             restrict_to_seg=restrict_to_seg)
            wt[:, :] = cwt[:, :]

        return wtmosaic

    def get_cweight_cutout_list(self, iobj, restrict_to_seg=False):
        """Composite the weight and seg maps, interpolating seg map from the
        coadd.

        The weight is set to zero outside the region as defined in the coadd.

        Parameters
        ----------
        iobj : int
            Index of the object.
        restrict_to_seg : bool, optional
            Set weights to zero for any pixel in the image not associated
            with the central object.

        Returns
        -------
        list : list of np.arrays
            A list of the weight maps.
        """
        wtmosaic = self.get_cweight_mosaic(
            iobj, restrict_to_seg=restrict_to_seg)

        # shares underlying storage
        wlist = split_mosaic(wtmosaic)
        return wlist

    def get_uberseg(self, iobj, icutout, fast=True):
        """Get the cweight map and zero out pixels not nearest to central
        object.

        Adapted from Niall Maccrann and Joe Zuntz.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int
            Index of cutout.
        fast : bool, optional
            Use the fast C code.

        Returns
        -------
        weight : np.array
            The weight map as a numpy array.
        """

        weight = self.get_cweight_cutout(iobj, icutout)
        seg = self.interpolate_coadd_seg(iobj, icutout)

        # if only have sky and object, then just return
        if len(numpy.unique(seg)) == 2:
            return weight

        # First get all indices of all seg map pixels which contain an object
        # i.e. are not equal to zero

        obj_inds = numpy.where(seg != 0)

        # the seg map holds the sextractor number, 1 offset
        object_number = self['number'][iobj]

        if fast and _have_c_ubserseg:
            # call fast c code with tree
            Nx, Ny = seg.shape
            Ninds = len(obj_inds[0])
            seg = seg.astype(numpy.int32)
            weight = weight.astype(numpy.float32, copy=False)
            obj_inds_x = obj_inds[0].astype(numpy.int32, copy=False)
            obj_inds_y = obj_inds[1].astype(numpy.int32, copy=False)
            _uberseg.uberseg_tree(
                seg, weight, Nx, Ny,
                object_number, obj_inds_x, obj_inds_y, Ninds)
        else:
            # Then loop through pixels in seg map, check which obj ind it is
            # closest to.  If the closest obj ind does not correspond to the
            # target, set this pixel in the weight map to zero.

            for i, row in enumerate(seg):
                for j, element in enumerate(row):
                    obj_dists = (i-obj_inds[0])**2 + (j-obj_inds[1])**2
                    ind_min = numpy.argmin(obj_dists)

                    segval = seg[obj_inds[0][ind_min], obj_inds[1][ind_min]]
                    if segval != object_number:
                        weight[i, j] = 0.

        return weight

    get_cweight_cutout_nearest = get_uberseg

    def get_uberseg_list(self, iobj, fast=True):
        """Get the cweight map and zero out pixels not nearest to central
        object.

        Adapted from Niall Maccrann and Joe Zuntz.

        Parameters
        ----------
        iobj : int
            Index of the object.
        fast : bool, optional
            Use the fast C code.

        Returns
        -------
        list : list of np.arrays
            A list of the weight maps.
        """

        useg_list = []
        for i in xrange(self['ncutout'][iobj]):
            uberseg = self.get_uberseg(
                iobj,
                i,
                fast=fast,
            )

            useg_list.append(uberseg)

        return useg_list

    get_cweight_cutout_nearest_list = get_uberseg_list

    def get_cseg_cutout(self, iobj, icutout):
        """Interpolate the coadd seg onto the plane of the cutout.

        The seg is set to zero outside the region as defined in the coadd,
        and to the "number" field from sextractor inside the region.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int
            Index of cutout.

        Returns
        -------
        seg : np.array
            The segmentation map.
        """
        seg = self.get_cutout(iobj, icutout, type='seg')
        seg[:, :] = self.get_number(iobj)

        coadd_seg = self.get_cutout(iobj, 0, type='seg')
        cseg = self._make_composite_image(iobj, icutout, seg, coadd_seg)
        return cseg

    def get_cseg_mosaic(self, iobj):
        """Interpolate the coadd seg onto the planes of the cutouts. Get
        a big mosaic of all of the epochs.

        The seg is set to zero outside the region as defined in the coadd,
        and to the "number" field from sextractor inside the region.

        Parameters
        ----------
        iobj : int
            Index of the object.

        Returns
        -------
        mosaic : np.array
            A mosaic of the segmentation maps.
        """
        segmosaic = self.get_mosaic(iobj, type='seg')
        segmosaic[:, :] = self.get_number(iobj)

        coadd_seg = self.get_cutout(iobj, 0, type='seg')

        # shares underlying storage
        wlist = split_mosaic(segmosaic)

        for icutout, seg in enumerate(wlist):
            cseg = self._make_composite_image(iobj, icutout, seg, coadd_seg)
            seg[:, :] = cseg[:, :]

        return segmosaic

    def get_cseg_cutout_list(self, iobj):
        """Interpolate the coadd seg onto the planes of the cutouts.
        Get a list of all seg cutouts.

        The seg is set to zero outside the region as defined in the coadd,
        and to the "number" field from sextractor inside the region.

        Parameters
        ----------
        iobj : int
            Index of the object.

        Returns
        -------
        list : list of np.arrays
            A list of the segmentation maps.
        """
        segmosaic = self.get_cseg_mosaic(iobj)

        # shares underlying storage
        seglist = split_mosaic(segmosaic)
        return seglist

    def get_cseg_weight(self, iobj, icutout, use_canonical_cen=False):
        """Get the largest circularly masked weight map that does not
        interesect any other objects seg map.

        If there are no other objects in the scene, the regular weight
        map is returned.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int
            Index of the cutout for this object.
        use_canonical_cen : bool, optional
            If `True`, use the canonical center of the image,
            (# of pixels - 1)/2, to compute the circular mask. The default of
            `False` uses the center recorded in the MEDS object data.

        Returns
        -------
        weight : np.array
            The masked weight map.
        """
        seg = self.get_cutout(iobj, icutout, type='seg')
        weight = self.get_cutout(iobj, icutout, type='weight')
        number = self['number'][iobj]

        wother = numpy.where((seg != 0) & (seg != number))
        if wother[0].size == 0:
            # no other objects in the stamp
            return weight

        if use_canonical_cen:
            row, col = (numpy.array(weight.shape)-1.0) / 2.0
        else:
            row = self['cutout_row'][iobj, icutout]
            col = self['cutout_col'][iobj, icutout]

        rows, cols = numpy.mgrid[
            0:weight.shape[0],
            0:weight.shape[1],
        ]

        rows = rows.astype('f8')
        cols = cols.astype('f8')

        rows -= row
        cols -= col

        r2 = rows**2 + cols**2

        minr2 = r2[wother].min()

        # now set the weight to zero for radii larger than that
        wkeep = numpy.where(r2 < minr2)
        new_weight = numpy.zeros(weight.shape)
        if wkeep[0].size > 0:
            new_weight[wkeep] = weight[wkeep]

        return new_weight

    def interpolate_coadd_seg_mosaic(self, iobj):
        """Get a mosaic of interpolated seg maps.

        Parameters
        ----------
        iobj : int
            Index of the object.

        Returns
        -------
        mosaic : np.array
            A mosaic of the segmentation maps.
        """
        seg_mosaic = self.get_mosaic(iobj, type='seg')

        segs = split_mosaic(seg_mosaic)

        for icutout, seg in enumerate(segs):
            if icutout == 0:
                continue

            seg[:, :] = self.interpolate_coadd_seg(iobj, icutout)
        return seg_mosaic

    def interpolate_coadd_seg(self, iobj, icutout):
        """Interpolate the coadd segmentation map onto the SE image frame.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int
            Index of cutout.

        Returns
        -------
        seg : np.array
            Interpolate segmentation map.
        """

        coadd_seg = self.get_cutout(iobj, 0, type='seg')
        if icutout == 0:
            return coadd_seg

        jmatrix = self.get_jacobian_matrix(iobj, icutout)
        cen = self.get_cutout_rowcol(iobj, icutout)

        return self._interpolate_coadd_seg_image(
            iobj,
            jmatrix,
            cen=cen,
        )

    def _interpolate_coadd_seg_image(self, iobj, jmatrix, cen=None):
        """Interpolate the coadd segmentation map onto the SE image frame.

        Parameters
        ----------
        iobj : int
            Index of the object.
        jmatrix: matrix
            The jacobian matrix for the image.
        """

        coadd_seg = self.get_cutout(iobj, 0, type='seg')

        seg = 0*coadd_seg.copy()

        if cen is None:
            cen = (numpy.array(seg.shape)-1.0)/2.0

        rowcen, colcen = cen

        coadd_jacob = self.get_jacobian_matrix(iobj, 0)
        coadd_rowcen, coadd_colcen = self.get_cutout_rowcol(iobj, 0)

        # rows in SE seg mape
        rows, cols = numpy.mgrid[0:seg.shape[0], 0:seg.shape[1]]
        rowsrel = rows-rowcen
        colsrel = cols-colcen

        # this will raise a numpy.linalg.linalg.LinAlgError exception
        cjinv = coadd_jacob.getI()

        # convert pixel coords in SE cutout to u,v
        u = rowsrel*jmatrix[0, 0] + colsrel*jmatrix[0, 1]
        v = rowsrel*jmatrix[1, 0] + colsrel*jmatrix[1, 1]

        # now convert into pixels for coadd
        crow = coadd_rowcen + u*cjinv[0, 0] + v*cjinv[0, 1]
        ccol = coadd_colcen + u*cjinv[1, 0] + v*cjinv[1, 1]

        crow = crow.round().astype('i8')
        ccol = ccol.round().astype('i8')

        # clipping makes the notation easier
        crow = crow.clip(0, coadd_seg.shape[0]-1)
        ccol = ccol.clip(0, coadd_seg.shape[1]-1)

        seg[rows, cols] = coadd_seg[crow, ccol]

        return seg

    def get_source_info(self, iobj, icutout):
        """Get the full source file information for the indicated cutout.

        Includes SE image and sky image.

        Parameters
        ----------
        iobj : int
            Index of the object.

        Returns
        -------
        image_info : tuple
            A tuple of the image info array entries.
        """
        self._check_indices(iobj, icutout=icutout)
        ifile = self._cat['file_id'][iobj, icutout]
        return self._image_info[ifile]

    def get_source_path(self, iobj, icutout):
        """Get the source filename associated with the indicated cutout.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int
            Index of the cutout for this object.

        Returns
        -------
        path : str
            The filename of the source image.
        """

        info = self.get_source_info(iobj, icutout)
        return info['image_path']

    def get_cat(self):
        """Get the catalog.
        """
        return self._cat

    def get_image_info(self):
        """Get the image information.
        """
        return self._image_info

    def get_meta(self):
        """Get the metadata.
        """
        return self._meta

    def get_jacobian(self, iobj, icutout):
        """Get the jacobian as a dictionary.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int
            Index of the cutout for this object.

        Returns
        -------
        jacobian_dict : dict
            A dictionary of the jacobian entries

                row0
                col0
                dudrow
                dudcol
                dvdcol
                dvdrow
        """
        self._check_indices(iobj, icutout=icutout)

        row0, col0 = self.get_cutout_rowcol(iobj, icutout)
        dudrow = self['dudrow'][iobj, icutout]
        dudcol = self['dudcol'][iobj, icutout]
        dvdrow = self['dvdrow'][iobj, icutout]
        dvdcol = self['dvdcol'][iobj, icutout]

        return {'row0': row0,
                'col0': col0,
                'dudrow': dudrow,
                'dudcol': dudcol,
                'dvdrow': dvdrow,
                'dvdcol': dvdcol}

    def get_jacobian_matrix(self, iobj, icutout):
        """Get the jacobian as a numpy array.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int integer
            Index of the cutout for this object.

        Returns
        -------
        jacobian : np.array
            A 2x2 matrix of the jacobian
                dudrow dudcol
                dvdrow dvdcol
        """
        jacob = numpy.matrix(numpy.zeros((2, 2)), copy=False)

        jacob[0, 0] = self['dudrow'][iobj, icutout]
        jacob[0, 1] = self['dudcol'][iobj, icutout]
        jacob[1, 0] = self['dvdrow'][iobj, icutout]
        jacob[1, 1] = self['dvdcol'][iobj, icutout]

        return jacob

    def get_jacobian_list(self, iobj):
        """Get the list of jacobians for all cutouts
        for this object.

        Parameters
        ----------
        iobj : int
            Index of the object.

        Returns
        -------
        jacobian_list : list of dicts
            A list of the jacobian entries as a dictionary.
        """
        self._check_indices(iobj)
        jlist = []
        for icutout in xrange(self['ncutout'][iobj]):
            j = self.get_jacobian(iobj, icutout)
            jlist.append(j)

        return jlist

    def get_number(self, iobj):
        """Get the segmentation map number.

        Parameters
        ----------
        iobj : int
            Index of the object.

        Returns
        -------
        number : int
            The segmentation map number.
        """
        if 'number' not in self._cat.dtype.names:
            return iobj+1
        else:
            return self._cat['number'][iobj]

    def get_cutout_rowcol(self, iobj, icutout):
        """Get cutout_row, cutout_col for the specified object
        and epoch.

        Parameters
        ----------
        iobj : int
            Index of the object.
        icutout : int integer
            Index of the cutout for this object.

        Returns
        -------
        row, col : tuple of floats
            The location in the cutout image.
        """
        row = self['cutout_row'][iobj, icutout]
        col = self['cutout_col'][iobj, icutout]

        return row, col

    def _make_composite_image(
            self, iobj, icutout, im, coadd_seg, restrict_to_seg=False):
        """
        Internal routine to composite the coadd seg onto another image,
        meaning set zero outside the region

        for the coadd this is easy, but for SE cutouts we need to use the
        jacobian to transform between SE and coadd coordinate systems
        """

        cim = im.copy()

        coadd_rowcen, coadd_colcen = self.get_cutout_rowcol(iobj, 0)
        rowcen, colcen = self.get_cutout_rowcol(iobj, icutout)

        segid = coadd_seg[int(coadd_rowcen), int(coadd_colcen)]

        if icutout == 0:
            # this cutout is the coadd
            logic = (coadd_seg != segid)
            if not restrict_to_seg:
                logic = logic & (coadd_seg != 0)

            # w = numpy.where( (coadd_seg != segid) & (coadd_seg != 0) )
            w = numpy.where(logic)
            if w[0].size != 0:
                cim[w] = 0.0
        else:
            rows, cols = numpy.mgrid[0:cim.shape[0], 0:cim.shape[1]]
            rows = rows-rowcen
            cols = cols-colcen

            se_jacob = self.get_jacobian_matrix(iobj, icutout)
            coadd_jacob = self.get_jacobian_matrix(iobj, 0)

            try:
                cjinv = coadd_jacob.getI()
            except numpy.linalg.linalg.LinAlgError:
                print('coadd jacobian is singular, setting weight to zero')
                cim[:, :] = 0.0
                return cim

            # convert pixel coords in SE cutout to u,v
            u = rows*se_jacob[0, 0] + cols*se_jacob[0, 1]
            v = rows*se_jacob[1, 0] + cols*se_jacob[1, 1]

            # now convert into pixels for coadd
            crow = coadd_rowcen + u*cjinv[0, 0] + v*cjinv[0, 1]
            ccol = coadd_colcen + u*cjinv[1, 0] + v*cjinv[1, 1]

            crow = crow.round().astype('i8')
            ccol = ccol.round().astype('i8')

            '''
            wbad=numpy.where(   (crow < 0) | (crow >= coadd_seg.shape[0])
                              & (ccol < 0) | (ccol >= coadd_seg.shape[1]) )
            if wbad[0].size != 0:
                cim[wbad] = 0
            '''
            # clipping makes the notation easier
            crow = crow.clip(0, coadd_seg.shape[0]-1)
            ccol = ccol.clip(0, coadd_seg.shape[1]-1)

            logic = (coadd_seg[crow, ccol] != segid)
            if not restrict_to_seg:
                logic = logic & (coadd_seg[crow, ccol] != 0)

            # wbad = numpy.where( (coadd_seg[crow,ccol] != segid ) &
            #                     (coadd_seg[crow,ccol] != 0) )
            wbad = numpy.where(logic)

            if wbad[0].size != 0:
                cim[wbad] = 0

        return cim

    def _get_extension_name(self, type):
        if type == 'image':
            return "image_cutouts"
        elif type == "weight":
            return "weight_cutouts"
        elif type == "seg":
            return "seg_cutouts"
        elif type == "bmask":
            return "bmask_cutouts"
        elif type == "ormask":
            return "ormask_cutouts"
        elif type == "noise":
            return "noise_cutouts"
        else:
            ext = "%s_cutouts" % type
            if ext not in self._fits:
                raise ValueError("bad cutout type '%s'" % type)
            else:
                return ext

    def _check_indices(self, iobj, icutout=None):
        if iobj >= self._cat.size:
            raise ValueError("object index should be within "
                             "[0,%s)" % self._cat.size)

        ncutout = self._cat['ncutout'][iobj]
        if ncutout == 0:
            raise ValueError("object %s has no cutouts" % iobj)

        if icutout is not None:
            if icutout >= ncutout:
                raise ValueError("requested cutout index %s for "
                                 "object %s should be in bounds "
                                 "[0,%s)" % (icutout, iobj, ncutout))

    def __repr__(self):
        return repr(self._fits['object_data'])

    def __getitem__(self, item):
        return self._cat[item]

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.close()

    @property
    def size(self):
        return self._cat.size


def split_mosaic(mosaic):
    """Split the mosaic into a list of images.

    The images in the list share memory with the original.

    Parameters
    ----------
    mosaic : np.array
        The mosaic image.

    Returns
    -------
    image_list : list of np.arrays
        List of the images.
    """
    box_size = mosaic.shape[1]
    ncutout = mosaic.shape[0]//box_size

    imlist = []
    for i in xrange(ncutout):
        r1 = i*box_size
        r2 = (i+1)*box_size
        imlist.append(mosaic[r1:r2, :])

    return imlist


def reject_outliers(imlist, wtlist, nsigma=5.0, A=0.3):
    """Set the weight for outlier pixels to zero.

    Algorithm
    ---------

    Reject pixels for which

     | im - med | > n*sigma_i + A*|med|

    where mu is the median

    I actually do

        wt*(im-med)**2 > (n + A*|med|*sqrt(wt))**2

    We wrongly assume the images all align, but this is ok as long as nsigma is
    high enough

    If the number of images is < 3 then the weight maps are not modified

    Credit
    ------
    Algorithm based on discussions with Daniel Gruen.

    Parameters
    ----------
    imlist : list of np.arrays
        List of images.
    wtlist : list of np.arrays
        List of weight maps for each image.
    nsigma : float, optional
        The minimum number of sigmas a pixel needs to deviate from the median
        to be flagged as an outlier.
    A : float
        The minimum fraction of the median a pixel needs to deviate from the
        median to be flagged as an outlier.

    Returns
    -------
    nreject : int
        The number of rejected pixels.
    """

    nreject = 0

    nim = len(imlist)
    if nim < 3:
        return nreject

    dims = imlist[0].shape
    imstack = numpy.zeros((nim, dims[0], dims[1]))

    for i, im in enumerate(imlist):
        imstack[i, :, :] = im

    med = numpy.median(imstack, axis=0)

    for i in xrange(nim):
        im = imlist[i]
        wt = wtlist[i]

        wt.clip(0.0, out=wt)

        ierr = numpy.sqrt(wt)

        # wt*(im-med)**2
        chi2_image = im.copy()
        chi2_image -= med
        chi2_image *= chi2_image
        chi2_image *= wt

        # ( n + A*|med|*sqrt(wt) )**2
        maxvals = numpy.abs(med)
        maxvals *= A
        maxvals *= ierr
        maxvals += nsigma
        maxvals *= maxvals

        w = numpy.where(chi2_image > maxvals)

        if w[0].size > 0:
            wt[w] = 0.0
            nreject += w[0].size

    return nreject
