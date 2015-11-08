"""
code to build MEDS files
"""
from __future__ import print_function
import json
import copy
import numpy
from numpy import where, zeros

# external requirements
import esutil as eu
import fitsio

from .util import \
        make_wcs_positions, \
        get_meds_output_struct, \
        get_meds_input_struct, \
        get_image_info_struct

from .bounds import Bounds
from .defaults import default_config


SUPPORTED_CUTOUT_TYPES = ['image','weight','seg','bmask']

class MEDSMaker(dict):
    """
    Write MEDS files.  See the docs at https://github.com/esheldon/meds
    for a description of the format

    parameters
    -----------
    obj_data: numpy array with fields
        Required fields are are 'id','box_size','ra','dec'.  For
        data types, see meds.util.get_meds_input_dtype
    image_info: numpy array with fields
        Information for each image.  For the required data type, see
        the meds.util.get_image_info_dtype() function.
    config: dict, optional
        Optional configuration parameters.  The available options
        are given XXX
    meta_data: numpy array with fields, optional
        Optional meta data to write.  This is typically a length
        one array, but can be anything in principle.
    """
    def __init__(self,
                 obj_data,
                 image_info,
                 config=None,
                 meta_data=None):

        self._load_config(config)
        self._set_extra_config()

        # make copies since we may alter some things
        self.image_info = image_info.copy()
        self.meta_data=copy.deepcopy(meta_data)
        self._set_obj_data(obj_data)

        self._force_box_sizes_even()


    def write(self, filename):
        """
        build the meds layout and write images
        """
        self._build_meds_layout()
        self._write_data(filename)

    def _write_data(self, filename):
        """
        run through and write cutouts from each SE file
        for each image

        Matt's notes
        1) figure out which objs overlap
        2) grab cutouts
        3) set weight maps properly (zero out bad pixels, areas off the chip)
        4) set bitmasks
        5) grab seg maps
        6) write to proper spot in each 1d image on disk

        """

        print("opening output MEDS file: '%s'" % filename)
        with fitsio.FITS(filename,'rw',clobber=True) as fits:
            self.fits=fits

            self._write_object_data()
            self._write_image_info()
            self._write_metadata()

            self._reserve_mosaic_images()

            self._write_cutouts('image')
            self._write_cutouts('weight')
            self._write_cutouts('seg')
            self._write_cutouts('bmask')

        print('output is in:',filename)

    def _write_object_data(self):
        """
        write the object data
        """

        print('writing object_data')
        self.fits.write(self.obj_data,
                        extname=self['object_data_extname'])

    def _write_image_info(self):
        """
        write the object data
        """

        print('writing image_info')
        self.fits.write(self.image_info,
                        extname=self['image_info_extname'])

    def _write_metadata(self):
        """
        write the object data
        """

        print('writing metadata')
        if self.meta_data is not None:
            self.fits.write(self.meta_data,
                            extname=self['metadata_extname'])


    def _reserve_mosaic_images(self):
        """
        reserve space on disk for each mosaic image
        """

        fits=self.fits

        dims=[self.total_pixels]

        for type in self['cutout_types']:
            print('    reserving %s mosaic' % type)
            extname=self['%s_cutout_extname' % type]
            dtype=self['%s_dtype' % type]
            fits.create_image_hdu(img=None,
                                  dtype=dtype,
                                  dims=dims,
                                  extname=extname)

    def _write_cutouts(self, cutout_type):
        """
        write the cutouts for the specified type
        """

        print('writing %s cutouts' % cutout_type)

        obj_data=self.obj_data

        nfile=self.image_info.size
        nobj=obj_data.size

        cutout_hdu = self._get_cutout_hdu(cutout_type)

        for file_id in xrange(nfile):

            pkey   = '%s_path' % cutout_type
            impath = self.image_info[pkey][file_id].strip()

            print('    %d/%d %s %s' % (file_id+1,nfile,cutout_type,impath))

            im_data = self._read_image(file_id, cutout_type)

            if im_data is None:
                print('    no %s specified for file' % cutout_type)
                continue

            for iobj in xrange(nobj):
                ncut=obj_data['ncutout'][iobj]

                for icut in xrange(ncut):
                    if obj_data['file_id'][iobj, icut] == file_id:
                        self._write_cutout(iobj,icut,cutout_hdu,im_data)

    def _write_cutout(self, iobj, icut, cutout_hdu, im_data):
        """
        extract a cutout and write it to the mosaic image
        """
        dims = im_data.shape

        d=self.obj_data
        orow = d['orig_start_row'][iobj,icut]
        ocol = d['orig_start_col'][iobj,icut]
        bsize = d['box_size'][iobj]
        start_row = d['start_row'][iobj,icut]

        orow_box, row_box = self._get_clipped_boxes(dims[0],orow,bsize)
        ocol_box, col_box = self._get_clipped_boxes(dims[1],ocol,bsize)

        read_im = im_data[orow_box[0]:orow_box[1],
                          ocol_box[0]:ocol_box[1]]

        subim = zeros( (bsize,bsize), dtype=read_im.dtype)

        subim[row_box[0]:row_box[1],
              col_box[0]:col_box[1]] = read_im

        cutout_hdu.write(subim, start=start_row)

    def _get_clipped_boxes(self, dim, start, bsize):
        """
        get clipped boxes for slicing

        If the box size goes outside the dimensions,
        trim them back

        parameters
        ----------
        dim: int
            Dimension of this axis
        start: int
            Starting position in the image for this axis
        bsize: int
            Size of box

        returns
        -------
        obox, box

        obox: [start,end]
            Start and end slice ranges in the original image
        box: [start,end]
            Start and end slice ranges in the output image
        """
        # slice range in the original image
        obox = [start, start+bsize]

        # slice range in the sub image into which we will copy
        box = [0, bsize]

        # rows
        if obox[0] < 0:
            obox[0] = 0
            box[0] = 0 - start

        im_max = dim
        diff= im_max - obox[1]
        if diff < 0:
            obox[1] = im_max
            box[1] = box[1] + diff

        return obox, box

    def _get_cutout_hdu(self, cutout_type):
        """
        get the cutout hdu object for the specified cutout type
        """
        tkey = '%s_cutouts' % cutout_type
        cutout_hdu = self.fits[tkey]
        return cutout_hdu

    def _read_image(self, file_id, cutout_type):
        """
        read an image, performing manipulations for
        some types

        images are background subtracted if a background
        file is specified.  The image is zerod where the
        bitmask is nonzero, if a bitmask file is specified.

        Similarly, weights are zerod where the bitmask is set.

        parameters
        ----------
        file_id: int
            The id into the image_info structure
        cutout_type: string
            'image','bkg','seg','bmask'
        """

        im = self._read_one_image(file_id, cutout_type)

        if cutout_type=='image':
            bkg = self._read_one_image(file_id, 'bkg')

            if bkg is not None:
                im -= bkg
            else:
                print('    no background for image')

            bmask = self._read_one_image(file_id, 'bmask')
            if bmask is not None:
                w=self._check_bad_bmask(bmask)
                im[w] = 0.0
            else:
                print('    no bmask for image')

            scale = self._get_scale(file_id)
            im *= scale

        elif cutout_type=='weight':

            bmask = self._read_one_image(file_id, 'bmask')

            if bmask is not None:
                w=self._check_bad_bmask(bmask)
                im[w] = 0.0
            else:
                print('    no bmask for image')

            scale = self._get_scale(file_id)
            im *= (1.0/scale**2)


        return im

    def _check_bad_bmask(self, bmask):
        """
        return indices with not-allowed bits set
        """

        binv = self['bitmask_allowed_inv']
        wbad = where( (bmask & binv) != 0)
        if wbad[0].size != 0:
            print('        found %d masked pixels' % wbad[0].size)
        return wbad


    def _read_one_image(self, file_id, cutout_type):
        """
        read a single image, no manipulations done here
        """
        info=self.image_info

        pkey = '%s_path' % cutout_type
        extkey = '%s_ext' % cutout_type
        impath=info[pkey][file_id].strip()
        ext = info[extkey][file_id]

        if impath.lower() == 'none' or impath=='':
            im=None
        else:
            if isinstance(ext, str):
                ext = ext.strip()

            with fitsio.FITS(impath) as fits:
                im = fits[ext].read()

        return im

    def _get_scale(self, file_id):
        """
        get the scale for the image if specified, else
        return 1.0
        """
        if 'scale' in self.image_info.dtype.names:
            return self.image_info['scale'][file_id]
        else:
            return 1.0


    def _build_meds_layout(self):
        """
        build the object data, filling in the stub we read

        note position offsets appear nowhere in this function
        """

        # box sizes are even
        half_box_size = self.obj_data['box_size']//2
        obj_data=self.obj_data

        nim  = self.image_info.size
        nobj = obj_data.size

        for file_id in xrange(nim):

            self._get_wcs(file_id)
            impath=self.image_info['image_path'][file_id].strip()
            position_offset=self.image_info['position_offset'][file_id]


            print("file %4d of %4d: '%s'" % (file_id+1,nim,impath))

            wcs = self._get_wcs(file_id)

            # monkey patching in the position_offset into wcs
            wcs.position_offset=position_offset

            q = self._do_rough_sky_cut(wcs, obj_data['ra'], obj_data['dec'])
            print('    first cut:  %6d of %6d objects' % (len(q),nobj))

            # this is the bottleneck
            pos = self._do_sky2image(wcs,
                                     obj_data['ra'][q],
                                     obj_data['dec'][q])

            # now test if in the actual image space.  Bounds are created
            # in the offset coords
            bnds = self._get_image_bounds(wcs)
            in_bnds = bnds.contains_points(pos['zrow'], pos['zcol'])
            q_rc, = numpy.where(in_bnds == True)
            print('    second cut: %6d of %6d objects' % (len(q_rc),len(q)))

            # compose them
            q = q[q_rc]

            # fill in the object_data structure

            # note q_rc since pos was created using obj_data[q]
            qrow = pos['zrow'][q_rc]
            qcol = pos['zcol'][q_rc]

            icut = obj_data['ncutout'][q]
            obj_data['file_id'][q,icut] = file_id
            obj_data['orig_row'][q,icut] = qrow
            obj_data['orig_col'][q,icut] = qcol

            #ostart_row = numpy.floor(qrow) - half_box_size[q]
            #ostart_col = numpy.floor(qcol) - half_box_size[q]
            # to be consistent with the C++ code
            ostart_row = qrow.astype('i4') - half_box_size[q]
            ostart_col = qcol.astype('i4') - half_box_size[q]
            crow       = qrow - ostart_row
            ccol       = qcol - ostart_col

            obj_data['orig_start_row'][q,icut] = ostart_row
            obj_data['orig_start_col'][q,icut] = ostart_col
            obj_data['cutout_row'][q,icut]     = crow
            obj_data['cutout_col'][q,icut]     = ccol

            # do jacobian, in original, not-offset coords
            # note q_rc since pos was created using obj_data[q]
            jacob = wcs.get_jacobian(pos['wcs_row'][q_rc],
                                     pos['wcs_col'][q_rc])

            # jacob is a tuple of arrays
            obj_data['dudcol'][q,icut] = jacob[0]
            obj_data['dudrow'][q,icut] = jacob[1]
            obj_data['dvdcol'][q,icut] = jacob[2]
            obj_data['dvdrow'][q,icut] = jacob[3]

            # increment
            obj_data['ncutout'][q] += 1

        self.obj_data = self._make_resized_data(obj_data)
        self._set_start_rows_and_pixel_count()

    def _set_start_rows_and_pixel_count(self):
        """
        set the total number of pixels in each mosaic
        """
        print('setting start rows and pixel count')
        data=self.obj_data
        nobj=data.size

        npix = (data['ncutout']*data['box_size']**2).sum()
        self.total_pixels = npix

        npix=0
        current_row = 0
        for iobj in xrange(nobj):
            ncut = data['ncutout'][iobj]
            if ncut > 0:
                bsize=data['box_size'][iobj]
                npix_per_cutout = bsize*bsize

                for icut in xrange(ncut):
                    data['start_row'][iobj,icut] = current_row
                    current_row += npix_per_cutout
                    npix += npix_per_cutout

        if self.total_pixels != npix:
            raise ValueError("total_pixels %d != "
                             "npix %d" % (self.total_pixels, npix))

        print('total pixels:',self.total_pixels)

    def _get_wcs(self, file_id):
        """
        either load the wcs from the image_info, or from
        the image header
        """
        if 'wcs' in self.image_info.dtype.names:
            wcs_string = self.image_info['wcs'][file_id]
            wcs_data = json.loads(wcs_string)
        else:
            impath=self.image_info['image_path'][file_id].strip()
            ext=self.image_info['image_ext'][file_id]
            wcs_data = fitsio.read_header(impath, ext=ext)

        wcs = eu.wcsutil.WCS(wcs_data)
        return wcs

    def _make_resized_data(self, odata):
        """
        make a new struct with ncutout-sized-arrays based on
        the actual maximum ncutout
        """

        nmax = odata['file_id'].shape[1]
        new_nmax = odata['ncutout'].max()
        if new_nmax < 2:
            new_nmax = 2
        temp_obj_data = odata

        nobj = temp_obj_data.size
        new_data = get_meds_output_struct(nobj, new_nmax,
                                          extra_fields=self['extra_fields'])

        tmpst = get_meds_output_struct(1, new_nmax)
        required_fields = tmpst.dtype.names

        for name in new_data.dtype.names:
            if name in temp_obj_data.dtype.names:

                lshape = len(new_data[name].shape)
                if lshape > 1 and name in required_fields:
                    new_data[name][:,:] = temp_obj_data[name][:,0:new_nmax]
                else:
                    new_data[name][:] = temp_obj_data[name][:]

        del temp_obj_data

        return new_data


    def _do_sky2image(self, wcs, ra, dec):
        """
        get image positions for the input radec. returns a structure
        with both wcs positions and zero offset positions
        """
        col,row = wcs.sky2image(ra,dec)

        positions = make_wcs_positions(row, col, wcs.position_offset)
        return positions


    def _do_rough_sky_cut(self, wcs, ra, dec):
        """
        rough sky bounds cut
        """

        sky_bnds = self._get_rough_sky_bounds(wcs)

        in_sky_bnds = sky_bnds.contains_points(ra, dec)
        q, = numpy.where(in_sky_bnds == True)

        return q

    def _get_rough_sky_bounds(self, wcs):
        """
        rough sky bounds for precut
        """
        ncol, nrow = wcs.get_naxis()

        # corners in default coord. system
        rows = numpy.array([0.0,    0.0, nrow-1, nrow-1])
        cols = numpy.array([0.0, ncol-1,    0.0, ncol-1])

        pos=make_wcs_positions(rows, cols, wcs.position_offset, inverse=True)
        ra,dec = wcs.image2sky(pos['wcs_col'], pos['wcs_row'])

        decrad = numpy.deg2rad(dec)
        rafac  = numpy.cos(decrad).min()

        rabuff  = self['bounds_buffer_radec']/rafac
        decbuff = self['bounds_buffer_radec']
        sky_bnds = Bounds(ra.min()  - rabuff,
                          ra.max()  + rabuff,
                          dec.min() - decbuff,
                          dec.max() + decbuff)
        return sky_bnds

    def _get_image_bounds(self, wcs):
        """
        separate out so we can make changes to offset code without
        altering calling function
        """

        ncol, nrow = wcs.get_naxis()

        rvals = numpy.array([1.0, nrow])
        cvals = numpy.array([1.0, ncol])

        pos = make_wcs_positions(rvals, cvals, wcs.position_offset)

        bnds = Bounds(pos['zrow'][0],
                      pos['zrow'][1],
                      pos['zcol'][0],
                      pos['zcol'][1])

        return bnds

    def _force_box_sizes_even(self):
        """
        box sizes are required to be even for MEDS files

        The DES maker will only make even box sizes, but eventually
        we will move this into the more general MEDSMaker that will
        take the catalogs as input
        """
        w,=numpy.where( (self.obj_data['box_size'] % 2) != 0)
        if w.size > 0:
            self.obj_data['box_size'][w] += 1


    def _set_cutout_types(self):

        cutout_types = copy.deepcopy(self['cutout_types'])

        # make sure 'image' is at the front
        if 'image' in cutout_types:
            cutout_types.remove('image')
        cutout_types = ['image'] + cutout_types

        bad_types=[]
        for ctype in cutout_types:
            if ctype not in SUPPORTED_CUTOUT_TYPES:
                bad_types.append(ctype)

        if len(bad_types) != 0:
            st=', '.join(bad_types)
            raise ValueError("unsupported cutout types: '%s'" % st)

        self['cutout_types'] = cutout_types
        print('writing cutouts for:',cutout_types)

    def _set_obj_data(self, obj_data):
        """
        copy the input data into a full object_data structure.

        check for required fields
        """
        self._check_required_obj_data_fields(obj_data)
        self.obj_data = self._get_full_obj_data(obj_data)

    def _check_required_obj_data_fields(self, obj_data):
        """
        make sure the input structure has the required fields
        """
        min_st = get_meds_input_struct(1)

        missing=[]
        for name in min_st.dtype.names:
            if name not in obj_data.dtype.names:
                missing.append(name)

        if len(missing) > 0:
            missing=', '.join(missing)
            raise ValueError("missing fields from obj_data: '%s'" % missing)

    def _get_full_obj_data(self, obj_data):
        """
        make a full object structure, adding in any extra fields from the
        input structure.  Copy over the common fields
        """
        nmax = self.image_info.size
        if nmax < 2:
            nmax = 2

        self._set_extra_fields(obj_data, nmax)

        nobj = obj_data.size
        new_obj_data = \
            get_meds_output_struct(nobj, nmax,
                                   extra_fields=self['extra_fields'])
        eu.numpy_util.copy_fields(obj_data, new_obj_data)

        return new_obj_data

    def _set_extra_fields(self, obj_data, nmax):
        """
        determine the tags in obj_data but not in the required
        fields for the output object_data
        """
        full_st = get_meds_output_struct(1, nmax)
        extra_fields = []

        for dt in obj_data.dtype.descr:
            name=dt[0]

            if name not in full_st.dtype.names:
                extra_fields.append(dt)

        self['extra_fields'] = extra_fields

    def _set_image_info(self, image_info):
        """
        set the image info and check for required fields
        """
        self._check_image_info(image_info)
        self.image_info = image_info.copy()


    def _check_image_info(self, image_info):
        """
        check required fields

        currently just make sure the structure is exactly
        like that in get_image_info_dtype
        """

        plen=2
        dt = numpy.dtype(get_image_info_dtype(plen))

        missing=[]
        for name in dt.names:
            if name not in image_info.dtype.names:
                missing.append(name)

        if len(missing) > 0:
            s=', '.join(missing)
            raise ValueError("missing image_info entries: '%s'" % s)

    def _set_extra_config(self):
        """
        set extra configuration parameters that are not user-controlled
        """
        self['object_data_extname']  = 'object_data'
        self['image_info_extname'] = 'image_info'
        self['metadata_extname'] = 'metadata'

        self['image_cutout_extname']  = 'image_cutouts'
        self['weight_cutout_extname'] = 'weight_cutouts'
        self['seg_cutout_extname']    = 'seg_cutouts'
        self['bmask_cutout_extname']  = 'bmask_cutouts'

    def _load_config(self, config):
        """
        load the default config, then load the input config
        """
        self.update(default_config)

        if config is not None:
            if not isinstance(config, dict):
                raise RuntimeError("config must be a dict, "
                                   "got %s" % type(config))
            self.update(config)

        self._set_cutout_types()

        # need this to be unsigned
        allowed = self['bitmask_allowed']
        allowed = numpy.array([allowed],dtype='u4')
        self['bitmask_allowed'] = allowed[0]
        self['bitmask_allowed_inv'] = ~allowed[0]
