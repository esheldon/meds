#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "meds.h"

int main(int argc, char **argv)
{
    if (argc < 3) {
        printf("usage: test meds_file type\n");
        printf("   type=1 just run through cutouts\n");
        printf("   type=2 run through cutouts and weights\n");
        printf("   type=3 run through cutouts, weights, and seg maps\n");
        exit(1);
    }

    const char *meds_file=argv[1];
    int type=atoi(argv[2]);

    printf("opening meds file: %s\n", meds_file);
    struct meds *meds=meds_open(meds_file);

    if (!meds) {
        fprintf(stderr,"error reading meds, exiting\n");
        exit(1);
    }

    struct meds_cutout *mosaic=NULL;
    struct meds_cutout *wmosaic=NULL;
    struct meds_icutout *smosaic=NULL;

    long nobj=meds_get_size(meds);
    printf("running through %ld objects\n", nobj);
    for (long iobj=0; iobj<nobj; iobj++) {
        if (meds_get_ncutout(meds, iobj) > 0) {

            mosaic=meds_get_mosaic(meds, iobj);
            if (type==2) {
                wmosaic=meds_get_weight_mosaic(meds, iobj);
            }
            if (type==3) {
                smosaic=meds_get_seg_mosaic(meds, iobj);
            }

            mosaic=meds_cutout_free(mosaic);
            if (type==2) {
                wmosaic=meds_cutout_free(wmosaic);
            }
            if (type==3) {
                smosaic=meds_icutout_free(smosaic);
            }
        }
    }

    meds=meds_free(meds);
}
