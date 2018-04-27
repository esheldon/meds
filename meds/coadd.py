"""
TODO
    - Add EDGEBLEED to edge flags, totally rejecting
    - adapt desmeds stuff to use finalcut only
    - consider not including edge cutouts in meds, or at least
      marking them somehow? Those pixels currently get zero weight.
    - tag versions of meds and psc
    - copy this code to some repo
        - only this FLAG_MAP is des specific, and is only
          used for printing.  So maybe the meds library
    - psf check
    - interp edge of ps pixels?
        - cutting when full edge masked
    - what to do for coadded bmask? currently oring
    - make psf a valid type for meds base maker
        - it is currently dealt with specially
"""
from __future__ import print_function
import os
import numpy as np
import fitsio
from . import util
import yaml
from . import maker
from .meds import reject_outliers
from .defaults import __version__

try:
    xrange=xrange
except:
    xrange=range


# flags for why a coadd was not made
NO_GOOD_CUTOUTS=2**0

# flags for why a cutout wasn't included
CUTOUT_HITS_EDGE=2**1
CUTOUT_MASK_FRAC=2**2
CUTOUT_CENTRAL_MASKED=2**3
CUTOUT_FULL_EDGE_MASKED=2**4
CUTOUT_ALL_WEIGHT_ZERO=2**5

"""
FLAG_MAP={
    "BPM":          1,  #/* set in bpm (hot/dead pixel/column)        */
    "SATURATE":     2,  #/* saturated pixel                           */
    "INTERP":       4,  #/* interpolated pixel                        */
    "BADAMP":       8,  #/* Data from non-functional amplifier        */
    "CRAY":        16,  #/* cosmic ray pixel                          */
    "STAR":        32,  #/* bright star pixel                         */
    "TRAIL":       64,  #/* bleed trail pixel                         */
    "EDGEBLEED":  128,  #/* edge bleed pixel                          */
    "SSXTALK":    256,  #/* pixel potentially effected by xtalk from  */
                        #/*       a super-saturated source            */
    "EDGE":       512,  #/* pixel flag to exclude CCD glowing edges   */
    "STREAK":    1024,  #/* pixel associated with streak from a       */
                        #/*       satellite, meteor, ufo...           */
    "SUSPECT":   2048,  #/* nominally useful pixel but not perfect    */
    "FIXED":     4096,  #/* corrected by pixcorrect                   */
    "NEAREDGE":  8192,  #/* suspect due to edge proximity             */
    "TAPEBUMP": 16384,  #/* suspect due to known tape bump            */
}
FLAG_MAP_INV={}
for flagname in FLAG_MAP:
    FLAG_MAP_INV[FLAG_MAP[flagname]] = flagname
"""

class MEDSCoaddMaker(maker.MEDSMaker):
    """
    Convert an existing MEDS file to a postage-stamp coadd
    MEDS file using the input coadder
    """
    def __init__(self, config, coadder, tmpdir=None):


        self.tmpdir=tmpdir

        self._load_config(config)
        self._set_extra_config()

        self.coadder=coadder

        # just a reference
        self.m = coadder.get_original_meds()

        # we use 2, even though only the coadd, so that
        # in fits they are still arrays
        self._set_obj_data()

        self.total_pixels = (self.m['box_size']**2).sum()

        self._set_image_info()
        self._set_meta()
        self._set_psf_layout()

    def _load_config(self, config):

        super(MEDSCoaddMaker,self)._load_config(config)
        for t in ['ormask','noise']:
            if t not in self['cutout_types']:
                self['cutout_types'] += [t]


    def write(self, filename, obj_range=None):
        """
        write all the data
        """

        print("opening output MEDS file: '%s'" % filename)
        with fitsio.FITS(filename,'rw',clobber=True) as fits:
            self.fits=fits

            self._write_image_info()
            self._write_metadata()

            self._reserve_mosaic_images()

            self._write_cutouts(obj_range=obj_range)

            # we fill this in as we go
            self._write_object_data()

        print('output is in:',filename)

    def _write_cutouts(self, obj_range=None):
        """
        get the coadd observations and write the
        data
        """
        if obj_range is None:
            obj_range=[0,self.m.size]

        self.last_row_start=0
        self.last_size=0
        self.last_psf_row_start=0
        self.last_psf_size=0

        for iobj in xrange(obj_range[0], obj_range[1]):
            print("%d %d/%d" % (iobj,iobj+1,obj_range[1]-1))

            self._write_object_cutouts(iobj)

    def _write_object_cutouts(self, iobj):
        """
        make the coadd and write all the data
        """
        icut=0
        d=self.obj_data
        coadd_obs, flags = self.coadder.get_coadd_obs(iobj)
        d['coadd_flags'][iobj] = flags
        if coadd_obs is not None:
            assert d['box_size'][iobj]==coadd_obs.image.shape[0]

            d['ncoadd'][iobj] = coadd_obs.meta['ncoadd']

            imflags=coadd_obs.meta['imflags']
            norig=len(imflags)
            d['orig_flags'][iobj,0:norig] = imflags
            d['start_row'][iobj,icut]=\
                    self.last_row_start+self.last_size
            d['psf_start_row'][iobj,icut]=\
                    self.last_psf_row_start+self.last_psf_size

            j = coadd_obs.jacobian
            pj = coadd_obs.psf.jacobian
            d['cutout_row'][iobj,icut] = j.row0
            d['cutout_col'][iobj,icut] = j.col0
            d['dudrow'][iobj,icut] = j.dudrow
            d['dudcol'][iobj,icut] = j.dudcol
            d['dvdrow'][iobj,icut] = j.dvdrow
            d['dvdcol'][iobj,icut] = j.dvdcol

            d['psf_box_size'][iobj] = coadd_obs.psf.image.shape[0]
            d['psf_cutout_row'][iobj,icut] = pj.row0
            d['psf_cutout_col'][iobj,icut] = pj.col0

            for cutout_type in self['cutout_types']+['psf']:
                self._write_cutout(
                    iobj,
                    icut,
                    coadd_obs,
                    cutout_type,
                )

            self.last_row_start     = d['start_row'][iobj,icut]
            self.last_size          = d['box_size'][iobj]**2
            self.last_psf_row_start = d['psf_start_row'][iobj,icut]
            self.last_psf_size      = d['psf_box_size'][iobj]**2

            d['ncutout'][iobj] = 1

    def _write_cutout(self, iobj, icut, coadd_obs, cutout_type):
        cutout_hdu = self._get_cutout_hdu(cutout_type)

        if cutout_type=='psf':
            im_data = coadd_obs.psf.image
            sname='psf_start_row'

        else:

            sname='start_row'
            if cutout_type=='image':
                im_data = coadd_obs.image
            elif cutout_type=='weight':
                im_data = coadd_obs.weight
            elif cutout_type=='seg':
                im_data = coadd_obs.seg
            elif cutout_type=='bmask':
                im_data = coadd_obs.bmask
            elif cutout_type=='ormask':
                im_data = coadd_obs.ormask
            elif cutout_type=='noise':
                im_data = coadd_obs.noise
            else:
                raise ValueError("bad cutout type: '%s'" % cutout_type)

        start_row=self.obj_data[sname][iobj,icut]
        cutout_hdu.write(im_data.ravel(), start=start_row)

    def _set_obj_data(self):
        nmax=2

        # subtract 1 for coadd
        nmax_orig=max(self.m['ncutout'].max()-1,2)
        extra_fields=[
            ('number','i8'),
            ('psf_box_size','i4'),
            ('psf_cutout_row','f8',nmax),
            ('psf_cutout_col','f8',nmax),
            ('psf_start_row','i8',nmax),
            ('orig_ncutout','i4'),
            ('orig_flags','i4',nmax_orig),
            ('ncoadd','i4'), # how many images got coadded
            ('coadd_flags','i4'),
        ]
        d=util.get_meds_output_struct(
            self.m.size,
            nmax,
            extra_fields=extra_fields,
        )

        cat=self.m.get_cat()

        d['ncutout']=0
        d['orig_ncutout'] = self.m['ncutout']-1 # -1 for original coadd
        d['ncoadd']=0
        d['file_id']=0
        d['coadd_flags']=0

        copy_fields=['id','box_size','ra','dec','number']
        for f in copy_fields:
            d[f] = cat[f]

        self.obj_data=d
        self.ncutout_max=nmax

    def _set_psf_layout(self):
        """
        set the box sizes and start row for each psf image
        """
        raise NotImplementedError("implement psf specifics in a sub class")

        '''
        obj_data = self.obj_data

        max_psf_size = self['max_psf_size']
        max_npixels_per = int(max_psf_size)**2
        nobj=self.m.size
        self.total_psf_pixels = max_npixels_per*nobj
        # fake the data for later
        self.psf_data=1
        '''

    def _set_image_info(self):
        """
        fake image info
        """
        self.image_info=util.get_image_info_struct(1, 10)

    def _set_meta(self):
        """
        set the meta data
        """
        import esutil as eu

        ometa=self.m.get_meta()

        dt=[]
        nv=np.__version__
        ev=eu.__version__
        fv=fitsio.__version__
        mv=__version__
        mfv=maker.MEDS_FMT_VERSION
        dt += [
            ('numpy_version','S%d' % len(nv)),
            ('esutil_version','S%d' % len(ev)),
            ('fitsio_version','S%d' % len(fv)),
            ('meds_version','S%d' % len(mv)),
            ('meds_fmt_version','S%d' % len(mfv)),
        ]
        names=[ d[0] for d in dt]

        for d in ometa.dtype.descr:
            n=d[0]
            if n not in names:
                dt.append( d )

        meta = np.zeros(1, dtype=dt)
        eu.numpy_util.copy_fields(ometa, meta)

        meta['numpy_version'] = nv
        meta['esutil_version'] = ev
        meta['fitsio_version'] = fv
        meta['meds_version'] = mv
        meta['meds_fmt_version'] = mfv

        self.meta_data=meta
 
class MEDSCoadder(dict):
    """
    generate postage stamp coadds and PSF coadds
    from the input MEDS object
    """
    def __init__(self,
                 config,
                 meds_obj,
                 psfmap,
                 seed,
                 make_plots=False):
        self._set_config(config)
        self.make_plots=make_plots

        self._set_bad_flags()
        self._set_edge_flags()
        self.rng=np.random.RandomState(seed)
        self.m = meds_obj
        self.psfmap = psfmap

        self._set_target_jacobian()

    def get_original_meds(self):
        """
        get a ref to the original meds object
        """
        return self.m

    def get_coadd_obs(self, iobj):
        """
        coadd the requested set of objects
        """
        import psc

        flags=0
        obslist, imflags = self._get_obslist(iobj)
        
        ncoadd=len(obslist)
        if ncoadd==0:
            flags |= NO_GOOD_CUTOUTS
            print("    nothing to coadd")
            return None, flags
    
        bmask=obslist[0].bmask*0
        ormask=obslist[0].bmask*0
        for obs in obslist:
            ormask |= obs.bmask

        coadder=psc.Coadder(
            obslist,
            jacobian=self.target_jacobian,
            **self['coadd']
        )
        coadd_obs = coadder.get_coadd()
        coadd_obs.bmask = bmask
        coadd_obs.ormask = ormask

        coadd_obs.seg = self._get_interpolated_coadd_seg(iobj)
        if self.make_plots:
            original_coadd=self.m.get_cutout(iobj, 0)
            original_seg=self.m.get_cutout(iobj, 0,type='seg')
            original_wt=self.m.get_cutout(iobj, 0,type='weight')
            self._show_coadd(iobj, coadd_obs,
                             original_coadd, original_wt, original_seg)

        meta={
            'ncoadd':ncoadd,
            'imflags':imflags,
        }
        coadd_obs.update_meta_data(meta)
        return coadd_obs, flags

    def _get_interpolated_coadd_seg(self, iobj):
        jmatrix = self._get_jacobian_matrix()
        return self.m.interpolate_coadd_seg_image(iobj, jmatrix)

    def _get_jacobian_matrix(self):
        j=self.target_jacobian
        return np.matrix(
            [ [j.dudrow, j.dudcol],
              [j.dvdrow, j.dvdcol], ]
        )

    def _set_config(self, config):
        self.update(config)

        self['coadd'] = self.get('coadd',{})

        self['dither_psfs'] = self['coadd'].pop('dither_psfs',True)


    def _set_bad_flags(self):
        if 'bmask_flags' not in self:
            bad_flags=None
        else:
            bad_flags = sum(self['bmask_flags'])

        self.bad_flags=bad_flags

    def _set_edge_flags(self):
        if 'edge_flags' not in self:
            edge_flags=None
        else:
            edge_flags = sum(self['bmask_flags'])

        self.edge_flags=edge_flags

    def _set_target_jacobian(self):
        import ngmix
        # center doesn't matter

        w,= np.where(self.m['ncutout'] > 1)
        dudrow = np.median(self.m['dudrow'][w,1])
        dudcol = np.median(self.m['dudcol'][w,1])
        dvdrow = np.median(self.m['dvdrow'][w,1])
        dvdcol = np.median(self.m['dvdcol'][w,1])

        self.target_jacobian=ngmix.Jacobian(
            row=15, col=15,
            dudrow=dudrow,
            dudcol=dudcol,
            dvdrow=dvdrow,
            dvdcol=dvdcol,
        )


    def _get_obslist(self, iobj):
        import ngmix

        m=self.m

        make_plots=self.make_plots

        obslist=ngmix.ObsList()

        imlist = m.get_cutout_list(iobj)
        if len(imlist)==1:
            # there is just the coadd; return the
            # empty obslist
            imflags=None
            return obslist, imflags

        wtlist = m.get_cutout_list(iobj, type='weight')
        bmlist = m.get_cutout_list(iobj, type='bmask')


        imlist=imlist[1:]
        wtlist=wtlist[1:]
        bmlist=bmlist[1:]

        if self['reject_outliers']:
            nreject = reject_outliers(
                imlist,
                wtlist,
            )
            if nreject > 0:
                print("    rejected",nreject)
            obslist.meta['nreject'] = nreject
        else:
            obslist.meta['nreject']=0

        ncutout = len(imlist)
        imflags = np.zeros(ncutout, dtype='i4')

        seglist = [m.interpolate_coadd_seg(iobj,i)
                   for i in range(1,ncutout+1)]

        # note starting at 1 assuming coadd is at 0
        for i in range(ncutout):
            icut=i+1

            print("    cutout %d/%d" % (i+1,ncutout))
            im=imlist[i]
            wt=wtlist[i]
            seg=seglist[i]
            bmask=bmlist[i]

            file_id=m['file_id'][iobj, icut]
            orig_row=m['orig_row'][iobj, icut]
            orig_col=m['orig_col'][iobj, icut]

            if self['symmetrize_mask']:
                self._symmetrize_weight(wt)
                self._symmetrize_bmask(bmask)

            if make_plots:
                self._show_epoch(iobj,icut,im,seg,wt,bmask)

            if self._hits_the_edge(bmask):
                print("    cutout hits the edge")
                imflags[i] |= CUTOUT_HITS_EDGE
                continue

            mask_frac = self._get_mask_frac(bmask, wt)
            if mask_frac > self['max_masked_frac']:
                print("    mask fraction %g exceeds "
                      "maximum %g" % (mask_frac, self['max_masked_frac']))
                imflags[i] |= CUTOUT_MASK_FRAC
                continue

            if self._central_region_is_masked(bmask, wt):
                print("    central region is masked")
                imflags[i] |= CUTOUT_CENTRAL_MASKED
                continue

            if self._ps_edge_pixels_are_masked(bmask, wt):
                print("    full edge of PS is masked")
                imflags[i] |= CUTOUT_FULL_EDGE_MASKED
                continue

            wt_logic = wt > 0.0
            w=np.where(wt_logic)
            if w[0].size == 0:
                print("    all weight is zero")
                imflags[i] |= CUTOUT_ALL_WEIGHT_ZERO
                continue

            jdict = m.get_jacobian(iobj, icut)
            jacobian = ngmix.Jacobian(
                row=jdict['row0'],
                col=jdict['col0'],
                dudrow=jdict['dudrow'],
                dudcol=jdict['dudcol'],
                dvdrow=jdict['dvdrow'],
                dvdcol=jdict['dvdcol'],
            )

            ccen = (np.array(im.shape)-1.0)/2.0
            offset_pixels=dict(
                row_offset=jacobian.row0-ccen[0],
                col_offset=jacobian.col0-ccen[1],
            )
            meta={'offset_pixels':offset_pixels}

            """
            print("image ccen:",ccen)
            print("orig_row: %.3f orig_col: %.3f" % (orig_row, orig_col))
            print("stamp_row: %.3f stamp_col: %.3f" % (jdict['row0'],jdict['col0']))
            """

            obs = ngmix.Observation(
                im,
                weight=wt,
                bmask=bmask,
                jacobian=jacobian,
                meta=meta,
            )

            obs.psf=self._get_psf_obs(
                obs, file_id, meta, orig_row, orig_col,
            )

            obs.seg=seg

            # interpolate im, weight, and a noise map
            # .noise attribute is added
            obs.meta['ninterp']=self._interp_bad_pixels(obs)

            if self.make_plots and obs.meta['ninterp'] > 0:
                self._show_epoch(iobj,icut,obs.image,seg,
                                 obs.weight,obs.bmask,
                                 extra='interp')
            obslist.append(obs)

        return obslist, imflags

    def _hits_the_edge(self, bmask):
        if self.edge_flags is not None:
            w=np.where( ((bmask & self.edge_flags) != 0) )
            hits_edge = w[0].size > 0
        else:
            hits_edge = False

        return hits_edge

    def _get_mask_frac(self, bmask, wt):
        logic=(wt <= 0.0)
        if self.bad_flags is not None:
            logic = logic | ((bmask & self.bad_flags) != 0)

            """
            for flag in FLAG_MAP_INV:
            #for flag in self['bmask_flags']:
                w=np.where( (bmask & flag) != 0 )
                if w[0].size > 0:
                    frac=w[0].size/float(bmask.size)
                    name=FLAG_MAP_INV[flag]
                    tup=(name,w[0].size,bmask.size,frac)
                    print("        %s %d/%d %g" % tup)
            """

        w=np.where(logic)
        mask_frac = float(w[0].size)/bmask.size
        return mask_frac

    def _symmetrize_weight(self, wt):
        """
        symmetrize zero weight pixels
        """
        assert wt.shape[0] == wt.shape[1]

        wt_rot=np.rot90(wt)
        wzero = np.where(wt_rot == 0.0)

        if wzero[0].size > 0:
            wt[wzero] = 0.0

    def _symmetrize_bmask(self, bmask):
        """
        symmetrize masked pixels
        """
        if self.bad_flags is None:
            return

        assert bmask.shape[0] == bmask.shape[1]

        bm_rot = np.rot90(bmask)

        wbad = np.where( (bm_rot & self.bad_flags) != 0)
        if wbad[0].size > 0:
            bmask[wbad] = bm_rot[wbad]

    def _central_region_is_masked(self, bmask, wt):
        cconf=self['central_region']
        if cconf['check']:
            rad=cconf['size']/2.0

            cen = (np.array(wt.shape)-1.0)/2.0
            low=int(round(cen[0]-rad))
            high=int(round(cen[0]+rad+1))

            logic = (wt[low:high, low:high] > 0.0)
            if self.bad_flags is not None:
                logic = logic & \
                    ((bmask[low:high, low:high] & self.bad_flags) == 0 )
            ok = np.all(logic)
        else:
            ok=True

        return not ok

    def _ps_edge_pixels_are_masked(self, bmask, wt):
        """
        we can't interpolate near the edge, so cut extreme
        cases
        """
        conf=self['stamp_edge']
        if conf['check']:
            bad_flags=self.bad_flags

            sides_bad = (
                np.all(wt[0, :]==0.0)
                or
                np.all(wt[-1,:]==0.0)
                or
                np.all(wt[:, 0]==0.0)
                or
                np.all(wt[:,-1]==0.0)
            )
            if not sides_bad:

                if bad_flags is not None:
                    sides_bad = (
                        np.all( (bmask[0, :] & bad_flags) != 0 )
                        or
                        np.all( (bmask[-1,:] & bad_flags) != 0 )
                        or
                        np.all( (bmask[:, 0] & bad_flags) != 0 )
                        or
                        np.all( (bmask[:,-1] & bad_flags) != 0 )
                    )


        return sides_bad


    def _get_psf_obs(self, obs, file_id, meta, row, col):
        """
        psfex specific code here

        for psfex we need to add 0.5 to get an offset
        that is the same as used for the object
        """
        raise NotImplementedError("implement psf specifics in a sub class")

    def _interp_bad_pixels(self, obs):
        """
        interpolate flagged pixels

        add a noise image that is also interpolated, if needed
        """

        iconf=self['interp']
        assert iconf['type']=="cubic","only cubic interpolation for now"

        im=obs.image
        weight=obs.weight

        if not obs.has_bmask():
            obs.bmask = np.zeros(im.shape, dtype='i4')

        bmask = obs.bmask

        bmravel = bmask.ravel()
        wtravel = weight.ravel()

        bad_logic = (wtravel <= 0.0)

        if self.bad_flags is not None:
            bad_logic = bad_logic | ( (bmravel & self.bad_flags) != 0 )

        wbad,=np.where(bad_logic)

        if wbad.size > 0:
            print("        interpolating %d/%d masked or zero weight "
                  "pixels" % (wbad.size,im.size))

            yy, xx = np.mgrid[0:im.shape[0], 0:im.shape[1]]

            x = xx.ravel()
            y = yy.ravel()

            yx = np.zeros( (x.size, 2) )
            yx[:,0] = y
            yx[:,1] = x

            wgood, = np.where(bad_logic == False)

            # make this a config option?
            medwt=np.median(wtravel[wgood])
            wtravel[wbad] = medwt

            im_interp = self._do_interp(yx, im, wgood, wbad)

            # still reference to ravelled one, should be ok
            tmp_noise = self._make_noise_image(weight)

            noise_interp = self._do_interp(yx, tmp_noise, wgood, wbad)

            obs.image  = im_interp
            obs.noise  = noise_interp

        else:
            obs.noise = self._make_noise_image(weight)
        
        return wbad.size

    def _make_noise_image(self, weight):
        """
        create a noise image based on the input weight map
        """

        err_image = np.sqrt( 1.0/weight )
        noise_image = self.rng.normal(size=weight.shape)
        noise_image *= err_image

        return noise_image


    def _do_interp(self, yx, im, wgood, wbad):
        import scipy.interpolate

        im_ravel = im.ravel()

        ii = scipy.interpolate.CloughTocher2DInterpolator(
            yx[wgood,:],
            im_ravel[wgood],
            fill_value=0.0,
        )

        im_interp = im.copy()
        im_interp_ravel = im_interp.ravel()

        vals = ii(yx[wbad,:])
        im_interp_ravel[wbad] = vals

        return im_interp


    def _show_epoch(self, iobj, icut, im, seg, wt, bmask, extra=None):
        import biggles
        import images

        tab=biggles.Table(2,2,aspect_ratio=1.0)

        tab[0,0] = images.view(im, title='im', show=False)
        tab[0,1] = images.view(seg, title='seg', show=False)
        tab[1,0] = images.view(wt, title='weight', show=False)
        tab[1,1] = images.view(bmask, title='bmask', show=False)
        
        d='plots'
        if not os.path.exists(d):
            try:
                os.makedirs(d)
            except:
                pass

        fname=[
            'object',
            '%06d' % iobj,
            '%02d' % icut,
        ]
        if extra is not None:
            fname += [extra]
        fname = '-'.join(fname)
        fname += '.png'
        fname=os.path.join(d, fname)

        #tab.show()
        print("    writing:",fname)
        tab.write_img(1500,1500,fname)

    def _show_coadd(self, iobj, obs, original_coadd, original_wt, original_seg):
        import biggles
        import images

        tab=biggles.Table(3,3,aspect_ratio=1.0)

        nl=None
        tab[0,0] = images.view(
            np.rot90( obs.image.transpose(),k=2),
            #obs.image,
            nonlinear=nl,
            title='coadd (trans/rot)',
            show=False,
        )
        tab[0,1] = images.view(original_coadd, nonlinear=nl, title='orig coadd', show=False)

        diff = np.rot90( obs.image.transpose(),k=2) - original_coadd
        tab[0,2] = images.view(diff, title='diff', show=False)

        cmin,cmax = obs.weight.min(), obs.weight.max()
        tab[1,0] = images.view(
            np.rot90( obs.weight.transpose(),k=2),
            #obs.weight,
            nonlinear=nl,
            title='weight (trans/rot) minmax: %.3g/%.3g' % (cmin,cmax),
            show=False,
        )
        cmin,cmax = original_wt.min(), original_wt.max()
        tab[1,1] = images.view(original_wt,
                               nonlinear=nl,
                               title='orig weight minmax: %.3g/%.3g' % (cmin,cmax),
                               show=False)

        tab[1,2] = images.view(obs.noise, title='noise', show=False)

        #tab[0,2] = images.view(obs.seg, title='seg', show=False)
        tab[2,0] = images.view(
            np.rot90( obs.seg.transpose(),k=2),
            #obs.seg,
            title='seg (trans/rot)',
            show=False,
        )
        tab[2,1] = images.view(original_seg, title='orig seg', show=False)

        tab[2,2] = images.view(obs.psf.image, title='coadded psf', show=False)
        
        d='plots'
        if not os.path.exists(d):
            try:
                os.makedirs(d)
            except:
                pass

        fname=[
            'object',
            '%06d' % iobj,
            'coadd',
        ]
        fname = '-'.join(fname)
        fname += '.png'
        fname=os.path.join(d, fname)

        #tab.show()
        print("    writing:",fname)
        tab.write_img(1500,1500,fname)


