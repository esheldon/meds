"""
Defines the MEDS class to work with MEDS (Multi Epoch Data Structures)

See docs for the MEDS class for more info

    Copyright (C) 2013, Erin Sheldon, BNL

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function
try:
    xrange
except:
    xrange=range

import numpy
import fitsio

try:
    from . import _uberseg
    _have_c_ubserseg=True
except ImportError:
    print("could not load fast ubserseg")
    _have_c_ubserseg=False

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

    get_cweight_cutout(iobj, icutout)
        Composite the weight and seg maps, interpolating seg map from the coadd
    get_cweight_mosaic(iobj)
        Composite the weight and seg maps, interpolating seg map from the coadd
        get all maps in a mosaic
    get_cweight_cutout_list(iobj)
        Composite the weight and seg maps, interpolating seg map from the coadd
        get all maps in a list

    get_cseg_cutout(iobj, icutout)
        Interpolate the coadd seg onto the plane of the cutout.
    get_cseg_mosaic(self, iobj)
        Interpolate the coadd seg onto the planes of the cutouts. Get
        a big mosaic of all.
    get_cseg_cutout_list(self, iobj)
        Interpolate the coadd seg onto the planes of the cutouts.
        Get a list of all seg cutouts.

    get_source_path(iobj, icutout)
        Get the source filename associated with the indicated cutout
    get_source_info(iobj, icutout)
        Get all info about the source images
    get_cat()
        Get the catalog; extension 1
    get_image_info()
        Get the entire image info structure
    get_meta()
        Get all the metadata
    get_jacobian(iobj, icutout)
        Get the jacobian as a dict
    get_jacobian_matrix(iobj, icutout)
        Get the jacobian as a numpy matrix
    get_jacobian_list(iobj)
        Get the list of jacobians for all cutouts for this object.

    stand alone function

    split_mosaic(mosaic)
        Split the mosaic into a list of images.
 


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

    # get a cutout for the bitmask map
    seg=m.get_cutout(35,3,type='bmask')

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
    meta=m.get_meta()


    Required Fields in main catalog
    --------------------------------

     id                 i8       id from coadd catalog
     ncutout            i8       number of cutouts for this object
     box_size           i8       box size for each cutout
     file_id            i8[NMAX] zero-offset id into the file names in the 
                                 second extension
     start_row          i8[NMAX] zero-offset, points to start of each cutout.
     orig_row           f8[NMAX] zero-offset position in original image
     orig_col           f8[NMAX] zero-offset position in original image
     orig_start_row     i8[NMAX] zero-offset start corner in original image
     orig_start_col     i8[NMAX] zero-offset start corner in original image
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
        self._filename=filename
        
        self._fits=fitsio.FITS(filename)

        self._cat=self._fits["object_data"][:]
        self._image_info=self._fits["image_info"][:]
        self._meta=self._fits["metadata"][:]

    def close(self):
        self._fits.close()
    
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
            values are 'image','weight','seg','bmask'

        returns
        -------
        The cutout image
        """

        if type=='psf':
            return self.get_psf(iobj,icutout)

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
            values are 'image','weight','seg'

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
            values are 'image','weight','seg'

        returns
        -------
        A list of images hold all cutouts.
        """

        if type=='psf':
            return self.get_psf_list(iobj)

        mosaic=self.get_mosaic(iobj,type=type)
        ncutout=self._cat['ncutout'][iobj]
        box_size=self._cat['box_size'][iobj]
        return split_mosaic(mosaic)

    def get_psf(self, iobj, icutout):
        """
        Get a single psf image for the indicated entry

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of the cutout for this object.

        returns
        -------
        The psf image
        """

        if 'psf' not in self._fits:
            raise RuntimeError("this MEDS file has no 'psf' extension")

        self._check_indices(iobj, icutout=icutout)

        cat=self._cat
        if len(cat['psf_box_size'].shape) > 1:
            box_size=self._cat['psf_box_size'][iobj,icutout]
        else:
            box_size=self._cat['psf_box_size'][iobj]

        start_row = self._cat['psf_start_row'][iobj,icutout]
        row_end = start_row + box_size*box_size

        imflat = self._fits['psf'][start_row:row_end]
        im = imflat.reshape(box_size,box_size)
        return im

    def get_psf_list(self, iobj):
        """
        get a list of psf images
        """
        ncut=self['ncutout'][iobj]
        return [self.get_psf(iobj, icut) for icut in xrange(ncut)]

    def get_cweight_cutout(self, iobj, icutout, restrict_to_seg=False):
        """
        Composite the weight and seg maps, interpolating seg map from the coadd

        The weight is set to zero outside the region as defined in the coadd

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of cutout

        returns
        -------
        The weight map
        """
        wt=self.get_cutout(iobj, icutout, type='weight')
        coadd_seg=self.get_cutout(iobj, 0, type='seg')
        cwt=self._make_composite_image(iobj, icutout, wt, coadd_seg,
                                       restrict_to_seg=restrict_to_seg)
        return cwt

    def get_cweight_mosaic(self, iobj, restrict_to_seg=False):
        """
        Composite the weight and seg maps, interpolating seg map from the coadd

        The weight is set to zero outside the region as defined in the coadd

        parameters
        ----------
        iobj:
            Index of the object

        returns
        -------
        A composite of all weight maps
        """
        wtmosaic=self.get_mosaic(iobj, type='weight')
        coadd_seg=self.get_cutout(iobj, 0, type='seg')

        ncutout=self._cat['ncutout'][iobj]
        box_size=self._cat['box_size'][iobj]

        # shares underlying storage
        wlist = split_mosaic(wtmosaic)

        for icutout,wt in enumerate(wlist):
            cwt=self._make_composite_image(iobj, icutout, wt, coadd_seg,
                                           restrict_to_seg=restrict_to_seg)
            wt[:,:] = cwt[:,:]

        return wtmosaic

    def get_cweight_cutout_list(self, iobj):
        """
        Composite the weight and seg maps, interpolating seg map from the coadd

        The weight is set to zero outside the region as defined in the coadd

        parameters
        ----------
        iobj:
            Index of the object

        returns
        -------
        A list containing all weight maps
        """
        wtmosaic=self.get_cweight_mosaic(iobj)

        ncutout=self._cat['ncutout'][iobj]
        box_size=self._cat['box_size'][iobj]

        # shares underlying storage
        wlist = split_mosaic(wtmosaic)
        return wlist



    def get_uberseg(self, iobj, icutout, fast=True):
        """
        get the cweight map and zero out pixels not nearest to central object

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of cutout

        returns
        -------
        The weight map

        Adapted from Niall Maccrann and Joe Zuntz
        """

        weight = self.get_cweight_cutout(iobj, icutout)
        seg    = self.interpolate_coadd_seg(iobj, icutout)

        #if only have sky and object, then just return
        if len(numpy.unique(seg)) == 2:
            return weight

        # First get all indices of all seg map pixels which contain an object
        # i.e. are not equal to zero

        obj_inds = numpy.where(seg != 0)

        # the seg map holds the sextractor number, 1 offset
        object_number = self['number'][iobj]

        if fast and _have_c_ubserseg:
            #call fast c code with tree
            Nx,Ny = seg.shape
            Ninds = len(obj_inds[0])
            seg = seg.astype(numpy.int32)
            weight = weight.astype(numpy.float32,copy=False)
            obj_inds_x = obj_inds[0].astype(numpy.int32,copy=False)
            obj_inds_y = obj_inds[1].astype(numpy.int32,copy=False)
            _uberseg.uberseg_tree(seg,weight,Nx,Ny,object_number,obj_inds_x,obj_inds_y,Ninds)
        else:
            # Then loop through pixels in seg map, check which obj ind it is closest
            # to.  If the closest obj ind does not correspond to the target, set this
            # pixel in the weight map to zero.

            for i,row in enumerate(seg):
                for j, element in enumerate(row):
                    obj_dists = (i-obj_inds[0])**2 + (j-obj_inds[1])**2
                    ind_min=numpy.argmin(obj_dists)

                    segval = seg[obj_inds[0][ind_min],obj_inds[1][ind_min]]
                    if segval != object_number:
                        weight[i,j] = 0.

        return weight

    get_cweight_cutout_nearest = get_uberseg

    def get_uberseg_list(self, iobj, fast=True):
        """
        Composite the weight and seg maps, interpolating seg map from the coadd

        The weight is set to zero outside the region as defined in the coadd

        parameters
        ----------
        iobj:
            Index of the object

        returns
        -------
        A list containing all weight maps
        """

        useg_list=[]
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
        """
        Interpolate the coadd seg onto the plane of the cutout.

        The seg is set to zero outside the region as defined in the coadd,
        and to the "number" field from sextractor inside the region.

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of cutout

        returns
        -------
        The seg map
        """
        seg=self.get_cutout(iobj, icutout, type='seg')
        seg[:,:] = self.get_number(iobj)

        coadd_seg=self.get_cutout(iobj, 0, type='seg')
        cseg=self._make_composite_image(iobj, icutout, seg, coadd_seg)
        return cseg


    def get_cseg_mosaic(self, iobj):
        """
        Interpolate the coadd seg onto the planes of the cutouts. Get
        a big mosaic of all.

        The seg is set to zero outside the region as defined in the coadd,
        and to the "number" field from sextractor inside the region.

        parameters
        ----------
        iobj:
            Index of the object

        returns
        -------
        A mosaic of all seg maps
        """
        segmosaic=self.get_mosaic(iobj, type='seg')
        segmosaic[:,:] = self.get_number(iobj)

        coadd_seg=self.get_cutout(iobj, 0, type='seg')

        ncutout=self._cat['ncutout'][iobj]
        box_size=self._cat['box_size'][iobj]

        # shares underlying storage
        wlist = split_mosaic(segmosaic)

        for icutout,seg in enumerate(wlist):
            cseg=self._make_composite_image(iobj, icutout, seg, coadd_seg)
            seg[:,:] = cseg[:,:]

        return segmosaic

    def get_cseg_cutout_list(self, iobj):
        """
        Interpolate the coadd seg onto the planes of the cutouts.
        Get a list of all seg cutouts.

        The seg is set to zero outside the region as defined in the coadd,
        and to the "number" field from sextractor inside the region.

        parameters
        ----------
        iobj:
            Index of the object

        returns
        -------
        A list containing all seg maps
        """
        segmosaic=self.get_cseg_mosaic(iobj)

        ncutout=self._cat['ncutout'][iobj]
        box_size=self._cat['box_size'][iobj]

        # shares underlying storage
        seglist = split_mosaic(segmosaic)
        return seglist

    def interpolate_coadd_seg_mosaic(self, iobj):
        """
        get a mosaic of interpolated seg maps
        """
        seg_mosaic=self.get_mosaic(iobj,type='seg')

        segs=split_mosaic(seg_mosaic)

        for icutout,seg in enumerate(segs):
            if icutout==0:
                continue

            seg[:,:] = self.interpolate_coadd_seg(iobj, icutout)
        return seg_mosaic

    def interpolate_coadd_seg(self, iobj, icutout):
        """
        interpolate the coadd segmentation map onto the SE image frame

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of cutout
        """

        coadd_seg = self.get_cutout(iobj, 0, type='seg')
        if icutout==0:
            return coadd_seg

        jmatrix = self.get_jacobian_matrix(iobj, icutout)
        cen = self.get_cutout_rowcol(iobj, icutout)

        return self._interpolate_coadd_seg_image(
            iobj,
            jmatrix,
            cen=cen,
        )

    def _interpolate_coadd_seg_image(self, iobj, jmatrix, cen=None):
        """
        interpolate the coadd segmentation map onto the SE image frame

        parameters
        ----------
        iobj:
            Index of the object
        jmatrix: matrix
            The jacobian matrix for the image
        """

        coadd_seg = self.get_cutout(iobj, 0, type='seg')

        seg = 0*coadd_seg.copy()

        if cen is None:
            cen = (numpy.array(seg.shape)-1.0)/2.0

        rowcen, colcen = cen

        coadd_jacob=self.get_jacobian_matrix(iobj, 0)
        coadd_rowcen, coadd_colcen = self.get_cutout_rowcol(iobj, 0)

        # rows in SE seg mape
        rows,cols=numpy.mgrid[0:seg.shape[0], 0:seg.shape[1]]
        rowsrel = rows-rowcen
        colsrel = cols-colcen

        # this will raise a numpy.linalg.linalg.LinAlgError exception
        cjinv = coadd_jacob.getI()

        # convert pixel coords in SE cutout to u,v
        u = rowsrel*jmatrix[0,0] + colsrel*jmatrix[0,1]
        v = rowsrel*jmatrix[1,0] + colsrel*jmatrix[1,1]

        # now convert into pixels for coadd
        crow = coadd_rowcen + u*cjinv[0,0] + v*cjinv[0,1]
        ccol = coadd_colcen + u*cjinv[1,0] + v*cjinv[1,1]

        crow = crow.round().astype('i8')
        ccol = ccol.round().astype('i8')

        # clipping makes the notation easier
        crow = crow.clip(0,coadd_seg.shape[0]-1)
        ccol = ccol.clip(0,coadd_seg.shape[1]-1)

        seg[rows, cols] = coadd_seg[crow, ccol]

        return seg


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

    def get_jacobian(self, iobj, icutout):
        """
        Get the jacobian as a dict keyed by

            row0
            col0
            dudrow
            dudcol
            dvdcol
            dvdrow

        parameters
        ----------
        iobj:
            Index of the object
        icutout:
            Index of the cutout for this object.
        """
        self._check_indices(iobj,icutout=icutout)

        row0, col0 = self.get_cutout_rowcol(iobj, icutout)
        dudrow = self['dudrow'][iobj,icutout]
        dudcol = self['dudcol'][iobj,icutout]
        dvdrow = self['dvdrow'][iobj,icutout]
        dvdcol = self['dvdcol'][iobj,icutout]

        return {'row0':row0,
                'col0':col0,
                'dudrow':dudrow,
                'dudcol':dudcol,
                'dvdrow':dvdrow,
                'dvdcol':dvdcol}

    def get_jacobian_matrix(self, iobj, icutout):
        """
        Get the jacobian as a numpy matrix

        parameters
        ----------
        iobj:
            Index of the object
        icutout: integer
            Index of the cutout for this object.

        returns
        -------
        A 2x2 matrix of the jacobian
            dudrow dudcol
            dvdrow dvdcol
        """
        jacob=numpy.matrix( numpy.zeros( (2,2) ), copy=False)

        jacob[0,0] = self['dudrow'][iobj,icutout]
        jacob[0,1] = self['dudcol'][iobj,icutout]
        jacob[1,0] = self['dvdrow'][iobj,icutout]
        jacob[1,1] = self['dvdcol'][iobj,icutout]

        return jacob


    def get_jacobian_list(self, iobj):
        """
        Get the list of jacobians for all cutouts
        for this object.

        parameters
        ----------
        iobj:
            Index of the object
        """
        self._check_indices(iobj)
        jlist=[]
        for icutout in xrange(self['ncutout'][iobj]):
            j=self.get_jacobian(iobj, icutout)
            jlist.append(j)

        return jlist

    def get_number(self, iobj):
        """
        Old versions of the meds files did not have number
        from the sextractor catalog
        """
        if 'number' not in self._cat.dtype.names:
            return iobj+1
        else:
            return self._cat['number'][iobj]

    def get_cutout_rowcol(self, iobj, icutout):
        """
        get cutout_row, cutout_col for the specified object
        and epoch

        parameters
        ----------
        iobj:
            Index of the object
        icutout: integer
            Index of the cutout for this object.

        returns
        -------
        row,col the location in the cutout image
        """
        row=self['cutout_row'][iobj,icutout]
        col=self['cutout_col'][iobj,icutout]

        return row, col

    def _make_composite_image(self, iobj, icutout, im, coadd_seg, restrict_to_seg=False):
        """
        Internal routine to composite the coadd seg onto another image,
        meaning set zero outside the region

        for the coadd this is easy, but for SE cutouts we need to use the
        jacobian to transform between SE and coadd coordinate systems
        """
        
        cim=im.copy()

        coadd_rowcen, coadd_colcen = self.get_cutout_rowcol(iobj, 0)
        rowcen, colcen = self.get_cutout_rowcol(iobj, icutout)

        segid=coadd_seg[int(coadd_rowcen),int(coadd_colcen)]

        if icutout==0:
            # this cutout is the coadd
            logic=(coadd_seg != segid)
            if not restrict_to_seg:
                logic = logic & (coadd_seg != 0)

            #w=numpy.where( (coadd_seg != segid) & (coadd_seg != 0) )
            w=numpy.where(logic)
            if w[0].size != 0:
                cim[w] = 0.0
        else:
            rows,cols=numpy.mgrid[0:cim.shape[0], 0:cim.shape[1]]
            rows = rows-rowcen
            cols = cols-colcen

            se_jacob=self.get_jacobian_matrix(iobj, icutout)
            coadd_jacob=self.get_jacobian_matrix(iobj, 0)

            try:
                cjinv = coadd_jacob.getI()
            except numpy.linalg.linalg.LinAlgError:
                print('coadd jacobian is singular, setting weight to zero')
                cim[:,:] = 0.0
                return cim

            # convert pixel coords in SE cutout to u,v
            u = rows*se_jacob[0,0] + cols*se_jacob[0,1]
            v = rows*se_jacob[1,0] + cols*se_jacob[1,1]

            # now convert into pixels for coadd
            crow = coadd_rowcen + u*cjinv[0,0] + v*cjinv[0,1]
            ccol = coadd_colcen + u*cjinv[1,0] + v*cjinv[1,1]

            crow = crow.round().astype('i8')
            ccol = ccol.round().astype('i8')

            '''
            wbad=numpy.where(   (crow < 0) | (crow >= coadd_seg.shape[0])
                              & (ccol < 0) | (ccol >= coadd_seg.shape[1]) )
            if wbad[0].size != 0:
                cim[wbad] = 0
            '''
            # clipping makes the notation easier
            crow = crow.clip(0,coadd_seg.shape[0]-1)
            ccol = ccol.clip(0,coadd_seg.shape[1]-1)

            logic = (coadd_seg[crow,ccol] != segid )
            if not restrict_to_seg:
                logic = logic & (coadd_seg[crow,ccol] != 0)

            #wbad=numpy.where( (coadd_seg[crow,ccol] != segid ) & (coadd_seg[crow,ccol] != 0) )
            wbad=numpy.where(logic)

            if wbad[0].size != 0:
                cim[wbad] = 0

        return cim


    def _get_extension_name(self, type):
        if type=='image':
            return "image_cutouts"
        elif type=="weight":
            return "weight_cutouts"
        elif type=="seg":
            return "seg_cutouts"
        elif type=="bmask":
            return "bmask_cutouts"
        elif type=="ormask":
            return "ormask_cutouts"
        elif type=="noise":
            return "noise_cutouts"
        else:
            raise ValueError("bad cutout type '%s'" % type)


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
    """
    Split the mosaic into a list of images.

    The images in the list share memory with the original.
    """
    box_size=mosaic.shape[1]
    ncutout = mosaic.shape[0]//box_size

    imlist=[]
    for i in xrange(ncutout):
        r1=i*box_size
        r2=(i+1)*box_size
        imlist.append( mosaic[r1:r2, :] )

    return imlist

def reject_outliers(imlist, wtlist, nsigma=5.0, A=0.3):
    """
    Set the weight for outlier pixels to zero

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

    credit
    ------
    Algorithm based on discussions with Daniel Gruen
    """

    nreject=0

    nim=len(imlist)
    if nim < 3:
        return nreject

    dims=imlist[0].shape
    imstack = numpy.zeros( (nim, dims[0], dims[1]) )

    for i,im in enumerate(imlist):
        imstack[i,:,:] = im

    med=numpy.median(imstack, axis=0)

    for i in xrange(nim):
        im=imlist[i]
        wt=wtlist[i]

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

        w=numpy.where(chi2_image > maxvals)

        if w[0].size > 0:
            wt[w] = 0.0
            nreject += w[0].size

    return nreject


