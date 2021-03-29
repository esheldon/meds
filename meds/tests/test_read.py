import os
import pytest
import tempfile
import numpy as np
import meds
from ._fakemeds import make_fake_meds


@pytest.mark.parametrize('with_psf', [False, True])
@pytest.mark.parametrize('with_bmask', [False, True])
@pytest.mark.parametrize('with_ormask', [False, True])
@pytest.mark.parametrize('with_noise', [False, True])
def test_medsreaders_smoke(with_psf, with_bmask, with_ormask, with_noise):

    rng = np.random.RandomState(542)

    # assert meds.meds._have_c_ubserseg

    cutout_types = ['image', 'weight', 'seg']
    if with_bmask:
        cutout_types += ['bmask']
    if with_ormask:
        cutout_types += ['ormask']
    if with_noise:
        cutout_types += ['noise']
    if with_psf:
        cutout_types += ['psf']

    with tempfile.TemporaryDirectory() as tdir:
        fname = os.path.join(tdir, 'test-meds.fits')
        make_fake_meds(
            fname=fname, rng=rng,
            cutout_types=cutout_types,
        )

        # test reading from single meds file
        m = meds.MEDS(fname)

        image_info = m.get_image_info()
        meta = m.get_meta()
        assert 'medsconf' in meta.dtype.names

        for iobj in range(m.size):
            box_size = m['box_size'][iobj]
            ncutout = m['ncutout'][iobj]

            imlist = m.get_cutout_list(iobj, type='image')
            mos = m.get_mosaic(iobj, type='image')
            assert mos.shape == (ncutout * box_size, box_size)

            cseg_mos = m.interpolate_coadd_seg_mosaic(iobj)
            assert cseg_mos.shape == (ncutout * box_size, box_size)

            cw_mos = m.get_cweight_mosaic(iobj)
            assert cw_mos.shape == (ncutout * box_size, box_size)

            assert len(imlist) == ncutout, ('iobj: %s' % iobj)

            if with_psf:
                assert m.has_psf()
                psf_list = m.get_cutout_list(iobj, type='psf')
                assert len(psf_list) == ncutout

                psf_list = m.get_psf_list(iobj)
                assert len(psf_list) == ncutout
            else:
                with pytest.raises(ValueError):
                    m.get_cutout_list(iobj, type='psf')

            if with_bmask:
                bmask_list = m.get_cutout_list(iobj, type='bmask')
                assert len(bmask_list) == ncutout
            else:
                with pytest.raises(ValueError):
                    m.get_cutout_list(iobj, type='bmask')

            if with_ormask:
                ormask_list = m.get_cutout_list(iobj, type='ormask')
                assert len(ormask_list) == ncutout
            else:
                with pytest.raises(ValueError):
                    m.get_cutout_list(iobj, type='ormask')

            if with_noise:
                noise_list = m.get_cutout_list(iobj, type='noise')
                assert len(noise_list) == ncutout
            else:
                with pytest.raises(ValueError):
                    m.get_cutout_list(iobj, type='noise')

            for ctype in cutout_types:
                for icut in range(ncutout):
                    cut = m.get_cutout(iobj, icut, type=ctype)
                    if ctype == 'psf':
                        psf_box_size = m['psf_box_size'][iobj]
                        assert cut.shape == (psf_box_size, ) * 2
                    else:
                        assert cut.shape == (box_size, ) * 2

                    if ctype == 'seg' and icut == 0:
                        # we put this in, not required
                        assert np.any(cut == m['number'][iobj])

            # uberseg

            for fast in [True, False]:
                # hmm... fast c code is inaccessible for some reason
                useglist = m.get_uberseg_list(iobj, fast=fast)
                assert len(useglist) == ncutout
                for icut in range(ncutout):
                    uberseg = m.get_uberseg(iobj, icut, fast=fast)
                    assert uberseg.shape == (box_size, ) * 2
                    assert np.any(uberseg != 0)

            cseg_list = m.get_cseg_cutout_list(iobj)
            assert len(cseg_list) == ncutout

            cwt_list = m.get_cweight_cutout_list(iobj)
            assert len(cwt_list) == ncutout

            for icut in range(ncutout):
                cseg_cutout = m.get_cseg_cutout(iobj, icut)
                assert cseg_cutout.shape == (box_size, ) * 2

                cs = m.interpolate_coadd_seg(iobj, icut)
                assert cs.shape == (box_size, ) * 2

                cw = m.get_cweight_cutout(iobj, icut)
                assert cw.shape == (box_size, ) * 2

                for use_canonical_cen in [True, False]:
                    cseg_weight = m.get_cseg_weight(
                        iobj, icut, use_canonical_cen=use_canonical_cen,
                    )
                    assert cseg_weight.shape == (box_size, ) * 2

            jlist = m.get_jacobian_list(iobj)
            assert len(jlist) == ncutout

            for icut in range(ncutout):
                source_info = m.get_source_info(iobj, icut)
                ifile = m['file_id'][iobj, icut]
                for f in source_info.dtype.names:
                    assert source_info[f] == image_info[f][ifile]
                source_path = m.get_source_path(iobj, icut)
                assert source_path == image_info['image_path'][ifile]

                j = m.get_jacobian(iobj, icut)
                row0, col0 = m.get_cutout_rowcol(iobj, icut)
                assert j['row0'] == row0
                assert j['col0'] == col0
                assert j['dudrow'] == m['dudrow'][iobj, icut]
                assert j['dudcol'] == m['dudcol'][iobj, icut]
                assert j['dvdrow'] == m['dvdrow'][iobj, icut]
                assert j['dvdcol'] == m['dvdcol'][iobj, icut]

                jm = m.get_jacobian_matrix(iobj, icut)
                assert jm[0, 0] == j['dudrow']
                assert jm[0, 1] == j['dudcol']
                assert jm[1, 0] == j['dvdrow']
                assert jm[1, 1] == j['dvdcol']

        # test an error that can occur
        with pytest.raises(ValueError):
            m.get_cutout_list(0, type='blah')
        with pytest.raises(ValueError):
            m.get_cutout_list(1000, 0)
        with pytest.raises(ValueError):
            m.get_cutout_list(0, 1000)


def test_outlier_rejection():
    rng = np.random.RandomState(11)
    with tempfile.TemporaryDirectory() as tdir:
        fname = os.path.join(tdir, 'test-meds.fits')
        make_fake_meds(fname=fname, rng=rng)

        m = meds.MEDS(fname)
        imax = m['ncutout'].argmax()

        imlist = m.get_cutout_list(imax)
        wtlist = m.get_cutout_list(imax, type='weight')

        dims = imlist[0].shape
        rowbad, colbad = ((np.array(dims)-1)/2).astype('i4')

        imlist[0][rowbad, colbad] = 1.e9

        meds.meds.reject_outliers(imlist, wtlist)
        assert wtlist[0][rowbad, colbad] == 0.0
