#include <stdlib.h>
#include <stdio.h>
#include "meds.h"

int main(int argc, char **argv)
{
    if (argc < 4) {
        printf("usage: meds-get-cutout meds_file outfile iobj [icutout]\n"
               "  extract a mosaic of all cutouts for the object of index\n"
               "  iobj.  If icutout is sent, only extract the cutout \n"
               "  at that index.\n"
               "  iobj are zero offset, equal to NUMBER-1 in the \n"
               "  catalog.  icutout are zero offset with the coadd cutout\n"
               "  at icutout==0\n");
        exit(1);
    }

    const char *meds_file=argv[1];
    const char *outfile=argv[2];
    int index=atoi(argv[3]);

    struct meds *meds=meds_open(meds_file);
    if (!meds) {
        fprintf(stderr,"error reading meds, exiting\n");
        exit(1);
    }

    struct meds_cutout *image = NULL;
    if (argc > 4) {
        int icutout=atoi(argv[4]);
        image = meds_get_cutout(meds, index, icutout);
    } else {
        image = meds_get_mosaic(meds, index);
    }

    if (image) {
        int status=0;
        int clobber=1;
        meds_cutout_write_fits(image, outfile, clobber, &status);
    } else {
        fprintf(stderr,"failed to extract cutouts\n");
    }

    image=meds_cutout_free(image);
    meds=meds_free(meds);
}
