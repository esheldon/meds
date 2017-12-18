"""
code for comparing two meds files
"""

import numpy
from .meds import MEDS


class Comparator(object):
    def __init__(self, file1, file2):
        self.m1 = MEDS(file1)
        self.m2 = MEDS(file2)
        self.tol=1.0e-4

    def go(self, nrand=None):
        """
        compare the two files
        """

        self._match()
        self.compare_images("image", nrand=nrand)

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
            print("%d/%d %d" % (icount+1,ntot,index))

            i1 = self.ind1[index]
            i2 = self.ind2[index]
            print("ncut:",self.m1['ncutout'][i1],self.m2['ncutout'][i2])
            print(i1,i2)
            ncut=self._check_ncutout(i1, i2)

            for icut in range(ncut):
                print("    reading im2",index,icut)
                im2 = self.m2.get_cutout(i2, icut, type=type)
                print("    reading im1",index,icut)
                im1 = self.m1.get_cutout(i1, icut, type=type)

                self._compare_images(im1, im2)

    def _check_ncutout(self, i1, i2):
        """
        check the number of cutouts agrees
        """
        n1 = self.m1['ncutout'][i1] 
        n2 = self.m2['ncutout'][i2] 
        if n1 != n2:
            raise ValueError("ncutout disagrees for objects %d:%d %d:%d" % (i1,n1,i2,n2))

        return n1

    def _compare_images(self, im1, im2):
        """
        compare the images
        """
        if im1.shape != im2.shape:
            raise ValueError("shapes do not match: %s %s" % (im1.shape, im2.shape))
        #print(  numpy.abs(im1 - im2).max() )

        diff = im1-im2
        check = numpy.abs(im1 - im2) < self.tol
        if not numpy.all(check):
            import images
            wbad = numpy.where(check == False)
            #print("%d matched worse than %g.  rms is %g" % (wbad[0].size, self.tol, diff.std()))

            images.compare_images(im1, im2, width=1500, height=1500)
            if 'q'==input("hit a key: (q to quit)"):
                stop

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
            self.m1['id'],
            self.m2['id'],
        )

        if self.ind1.size == 0:
            raise RuntimeError("no ids matched")
        else:
            print("%d/%d matched" % (self.ind1.size, self.m1.size))
