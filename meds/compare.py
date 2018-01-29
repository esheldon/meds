"""
code for comparing two meds files
"""

from __future__ import print_function

import numpy
from .meds import MEDS

class Comparator(object):
    def __init__(self, file1, file2, id_name='id'):
        self.m1 = MEDS(file1)
        self.m2 = MEDS(file2)
        self.id_name = id_name

        self.image_rms_tol = 0.8
        self.weight_tol = 1.0e-4

    def go(self, nrand=None):
        """
        compare the two files
        """

        self._match()

        self.compare_object_data()
        for type in ['image','weight','bmask','seg']:
            print("checking",type)
            self.compare_images(type, nrand=nrand)

    def compare_object_data(self):
        """
        compare all fields, requiring exact match
        """

        print("checking object data")
        c1 = self.m1.get_cat()
        c2 = self.m2.get_cat()

        if c1.dtype.names != c2.dtype.names:
            raise ValueError("cat names don't match")

        for n in c1.dtype.names:
            check = numpy.allclose(
                c1[n][self.ind1],
                c2[n][self.ind2],
            )
            if not check:
                import esutil as eu
                #raise ValueError("field %s does not match" % n)
                diff = (c1[n][self.ind1] - c2[n][self.ind2]).ravel()
                rms = diff.std()
                mn_clip, rms_clip = eu.stat.sigma_clip(diff,silent=True)
                print("field %s does not match.  rms %g "
                      "clipped: %g" % (n, rms,rms_clip))

    def compare_images(self, type, nrand=None):
        """
        compare images of the specified type

        parameters
        ----------
        type: string
            Type of images to compare, e.g. 'image','weight' etc.
        nrand: integer, optional
            Compare at most this many images.  A random subset will be
            drawn
        """

        if nrand is not None:
            nrand = min(nrand, self.ind1.size)
            ids = self._get_random_subset(nrand)
        else:
            ids = numpy.arange(self.ind1.size)

        ntot = ids.size
        for icount, index in enumerate(ids):
            print("    %d/%d %d" % (icount+1,ntot,index))

            i1 = self.ind1[index]
            i2 = self.ind2[index]
            ncut=self._check_ncutout(i1, i2)

            for icut in range(ncut):
                im2 = self.m2.get_cutout(i2, icut, type=type)
                im1 = self.m1.get_cutout(i1, icut, type=type)
                #if type=='image':
                #    wt=self.m2.get_cutout(i2, icut, type='weight')
                #    err=numpy.sqrt(1.0/wt.max())
                #    print("        compare to noise:",err)

                if type=='image':
                    self._compare_images(im1, im2)
                elif type == 'weight':
                    self._compare_weights(im1, im2)
                else:
                    self._compare_images_exact(im1, im2)

    def _check_ncutout(self, i1, i2):
        """
        check the number of cutouts agrees
        """
        n1 = self.m1['ncutout'][i1] 
        n2 = self.m2['ncutout'][i2] 
        if n1 != n2:
            raise ValueError("ncutout disagrees for objects "
                             "%d:%d %d:%d" % (i1,n1,i2,n2))

        return n1

    def _compare_images(self, im1, im2):
        """
        compare the images.

        We only check the rms
        """
        if im1.shape != im2.shape:
            raise ValueError("shapes do not "
                             "match: %s %s" % (im1.shape, im2.shape))

        diff = im1-im2
        std = diff.std()
        if std > self.image_rms_tol:
            print("        rms %g greater than "
                  "tolerance %g" % (std,self.image_rms_tol))

            if False:
                import images
                images.compare_images(im1, im2, width=1500, height=1500)
                if 'q'==raw_input("hit a key: (q to quit)"):
                    stop


    def _compare_weights(self, im1, im2):
        """
        compare the images
        """
        if im1.shape != im2.shape:
            raise ValueError("shapes do not "
                             "match: %s %s" % (im1.shape, im2.shape))

        diff = im1-im2
        check = numpy.abs(diff) < self.weight_tol
        if not numpy.all(check):
            import images
            wbad = numpy.where(check == False)
            print("        %d matched worse than %g.  rms is %g" % \
                  (wbad[0].size, self.weight_tol, diff.std()))

    def _compare_images_exact(self, im1, im2):
        """
        compare the images
        """
        if im1.shape != im2.shape:
            raise ValueError("shapes do not "
                             "match: %s %s" % (im1.shape, im2.shape))

        wbad=numpy.where(im1 != im2)

        if wbad[0].size > 0:
            raise ValueError("%d/%d pixels did "
                             "not match" % (wbad[0].size, im1.size))

    def _get_random_subset(self, nrand):
        """
        get random unique subset of range [0,n)
        """
        import esutil as eu
        return eu.random.random_indices(self.ind1.size, nrand)

    def _match(self):
        """
        match by id
        """
        import esutil as eu
        self.ind1, self.ind2 = eu.numpy_util.match(
            self.m1[self.id_name],
            self.m2[self.id_name],
        )

        if self.ind1.size == 0:
            raise RuntimeError("no ids matched")
        else:
            print("%d/%d matched" % (self.ind1.size, self.m1.size))
