"""
MEDSNumberExtractor
    A class to extract a subset of the objects in a MEDS file
    and write to a new file using only the object numbers
"""
from __future__ import print_function

import os
import fitsio
import numpy
from .extractor import get_psf_row_range


def extract_numbers(meds_file, numbers, sub_file):
    """
    Extract a subset of objects and write a new meds file.

    If you want this as a temporary file, which will be cleaned
    when you are done with it, use a MEDSNumberExtractor object with
    cleanup=True
    """

    MEDSNumberExtractor(meds_file, numbers, sub_file)


class MEDSNumberExtractor(object):
    """
    Class to extract a subset of objects and write a new meds file.

    Optionally clean up the new file when the object is destroyed.
    """

    def __init__(self, meds_file, numbers, sub_file, cleanup=False):
        self.meds_file = meds_file
        nums = numpy.array(numbers, dtype=int)
        unums = numpy.unique(nums)
        assert len(nums) == len(unums), (
            "Input numbers must be unique! len(numbers) = %ld, len(unique(numbers)) = %ld"  # noqa
            % (len(nums), len(unums))
        )
        self.numbers = nums
        self.sub_file = sub_file
        self.cleanup = cleanup
        self._check_inputs()

        q = numpy.argsort(self.numbers)
        self.numbers = self.numbers[q]

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
                print("removing sub file:", self.sub_file)
                os.remove(self.sub_file)

    def _get_inds(self, data):
        inds = []
        for number in self.numbers:
            if (
                number <= len(data["number"])
                and number >= 1
                and number == data["number"][number - 1]
            ):
                q = [number - 1]
            else:
                (q,) = numpy.where(number == data["number"])
            assert (
                len(q) == 1
            ), "Could not find or found duplicate number: number = %ld, found %ld" % (  # noqa
                number,
                len(q),
            )
            inds.append(q[0])
        inds = numpy.array(inds, dtype=int)
        return inds

    def _get_row_ranges(self, data):
        """
        get pixel range for this subset
        """
        (w,) = numpy.where(data["ncutout"] > 0)
        if w.size == 0:
            return [(-1, -1)]

        ranges = []
        start = True
        for w in range(len(data)):
            if data["ncutout"][w] > 0:
                if start:
                    ranges.append(
                        [
                            data["start_row"][w, 0],
                            data["start_row"][w, 0]
                            + data["box_size"][w] ** 2 * data["ncutout"][w],
                        ]
                    )
                    start = False
                else:
                    if data["start_row"][w, 0] != ranges[-1][-1]:
                        # add new if not back to back
                        ranges.append(
                            [
                                data["start_row"][w, 0],
                                data["start_row"][w, 0]
                                + data["box_size"][w] ** 2 * data["ncutout"][w],  # noqa
                            ]
                        )
                    else:
                        # just expand range of pixels if back to back
                        ranges[-1][-1] = (
                            data["start_row"][w, 0]
                            + data["box_size"][w] ** 2 * data["ncutout"][w]
                        )

        if len(ranges) == 0:
            ranges.append([-1, -1])

        # send back as list of tuples
        tup_ranges = [tuple(rng) for rng in ranges]
        return tup_ranges

    def _get_data_from_ranges(self, ranges, data):
        out_data = []
        dtype = None
        for rng in ranges:
            if dtype is None:
                arr = data[rng[0]: rng[1]]
                dtype = arr.dtype
            out_data.extend(list(data[rng[0]: rng[1]]))
        out_data = numpy.array(out_data, dtype=dtype)
        return out_data

    def _extract(self):

        with fitsio.FITS(self.meds_file) as infits:
            print("opening sub file:", self.sub_file)
            with fitsio.FITS(self.sub_file, "rw", clobber=True) as outfits:

                #
                # subset of object data table
                #
                inds = self._get_inds(infits["object_data"][:])
                obj_data = infits["object_data"][inds]

                ranges = self._get_row_ranges(obj_data)
                if ranges[0][0] != -1:
                    # adjust to new start.
                    # If ranges[0][0]==-1 will all be -9999
                    loc = 0
                    for i in range(len(obj_data)):
                        for j in range(obj_data["ncutout"][i]):
                            obj_data["start_row"][i, j] = loc
                            loc += obj_data["box_size"][i] ** 2

                if "psf" in infits:
                    psf_ranges = get_psf_row_ranges(obj_data)
                    if psf_ranges[0][0] != -1:
                        # adjust to new start.
                        # If ranges[0][0]==-1 will all be -9999
                        loc = 0
                        for iobj in range(len(obj_data)):
                            for icut in range(obj_data["ncutout"][iobj]):

                                obj_data["psf_start_row"][iobj, icut] = loc

                                rs = obj_data["psf_row_size"][iobj, icut]
                                cs = obj_data["psf_col_size"][iobj, icut]
                                loc += rs * cs

                outfits.write(obj_data, extname="object_data")

                #
                # copy all metadata and image info
                #
                iinfo = infits["image_info"][:]
                outfits.write(iinfo, extname="image_info")

                meta = infits["metadata"][:]
                outfits.write(meta, extname="metadata")

                #
                # extract cutouts for the requested objects
                #
                if ranges[0][0] == -1:
                    self._write_dummy(outfits)
                else:
                    image_cutouts = self._get_data_from_ranges(
                        ranges, infits["image_cutouts"]
                    )
                    outfits.write(image_cutouts, extname="image_cutouts")
                    del image_cutouts

                    weight_cutouts = self._get_data_from_ranges(
                        ranges, infits["weight_cutouts"]
                    )
                    outfits.write(weight_cutouts, extname="weight_cutouts")
                    del weight_cutouts

                    seg_cutouts = self._get_data_from_ranges(
                        ranges, infits["seg_cutouts"]
                    )
                    outfits.write(seg_cutouts, extname="seg_cutouts")
                    del seg_cutouts

                    bmask_cutouts = self._get_data_from_ranges(
                        ranges, infits["bmask_cutouts"]
                    )
                    outfits.write(bmask_cutouts, extname="bmask_cutouts")
                    del bmask_cutouts

                    if "psf" in infits:
                        psf_cutouts = self._get_data_from_ranges(
                            psf_ranges, infits["psf"]
                        )
                        outfits.write(psf_cutouts, extname="psf")
                        del psf_cutouts

    def _write_dummy(self, outfits):
        print("no objects with cutouts, writing dummy data")
        dummy = numpy.zeros(2, dtype="f4") + -9999
        outfits.write(dummy, extname="image_cutouts")
        dummy = numpy.zeros(2, dtype="f4")
        outfits.write(dummy, extname="weight_cutouts")
        dummy = numpy.zeros(2, dtype="i4") + -9999
        outfits.write(dummy, extname="seg_cutouts")
        dummy = numpy.zeros(2, dtype="i4") + -9999
        outfits.write(dummy, extname="bmask_cutouts")

    def _check_inputs(self):
        if self.meds_file == self.sub_file:
            raise ValueError("output file name equals input")

        if len(self.numbers) == 0:
            raise ValueError("one must extract at least one object")


def get_psf_row_ranges(data):
    """
    get pixel range for this subset
    """
    (w,) = numpy.where(data["ncutout"] > 0)
    if w.size == 0:
        return [(-1, -1)]

    ranges = []
    for w in range(len(data)):
        if data["ncutout"][w] > 0:
            rr = get_psf_row_range(data[w: w + 1])
            ranges.append(rr)

    if len(ranges) == 0:
        ranges.append([-1, -1])

    # send back as list of tuples
    tup_ranges = [tuple(rng) for rng in ranges]
    return tup_ranges
