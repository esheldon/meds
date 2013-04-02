/*
   Basic library to work with a Multi Epoch Data Structure
    https://cdcvs.fnal.gov/redmine/projects/deswlwg/wiki/Multi_Epoch_Data_Structure

   If you have your own image library, you can use the meds_get_cutoutp and
   meds_get_mosaicp to get raw pointers that can be wrapped by your library.

   You can also use image structures returned by meds_get_cutout and
   meds_get_mosaic

   A few examples (for more info, see actual definitions in this header)
   ---------------------------------------------------------------------

   #include <meds.h>

   struct meds *meds=meds_open(filename);
   meds_print(meds, stdout);

   long nobj=meds_get_size(meds);

   long iobj=35;
   long ncutout=meds_get_ncutout(meds, iobj);
   if (ncutout == 0) {
       printf("no cutouts for object %ld\n", iobj);
       ....
   }

   //
   // get the 4th cutout for this object (at index=3)
   //

   long icutout=3;

   // using a cutout structure.  Will return NULL if not found
   struct meds_cutout *cutout=meds_get_cutout(meds, iobj, icutout);

   printf("nrow: %ld ncol: %ld\n", 
       CUTOUT_NROW(cutout), CUTOUT_NCOL(cutout));

   long row=5, col=8;
   printf("pixel [%ld,%ld]: %g\n", 
       row, col, CUTOUT_GET(cutout, row, col));

   // get the center pixel in the cutout
   double rowcen=0,colcen=0;
   meds_get_cutout_cen(meds, iobj, icutout, &rowcen, &colcen);
   printf("center pixel [%lf,%lf]: %g\n", 
       rowcen, colcen, CUTOUT_GET(cutout, (int)rowcen, (int)colcen));


   // get the weight image
   struct meds_cutout *wcutout=meds_get_weight_cutout(meds, iobj, icutout);

   // get the seg map: note type is meds_icutout
   struct meds_icutout *scutout=meds_get_seg_cutout(meds, iobj, icutout);
   printf("center pixel seg [%lf,%lf]: %g\n", 
       rowcen, colcen, CUTOUT_GET(scutout, (int)rowcen, (int)colcen));

   // get struct with input coadd file path, etc.
   const struct meds_meta *meta=meds_get_meta(meds);

   // using the pointer interface; good if you have your own image library
   long nrow=0, ncol=0;
   double *pix=meds_get_cutoutp(meds, iobj, icutout, &nrow, &ncol);
   double *wpix=meds_get_weight_cutoutp(meds, iobj, icutout, &nrow, &ncol);


   // what file did the cutout come from?
   // works for cutouts and weights, which come from the same file
   const char *name=meds_get_source_path(meds, iobj, icutout);
   // sky image path
   const char *skyname=meds_get_sky_path(meds, iobj, icutout);


   //
   // get a mosaic of all cutouts for this object
   //

   // using the cutout structure
   struct meds_cutout *mosaic=meds_get_mosaic(meds, iobj);
   struct meds_cutout *wmosaic=meds_get_weight_mosaic(meds, iobj);

   printf("ncutout: %ld\n", MOSAIC_NCUTOUT(mosaic));
   printf("per cutout, nrow: %ld ncol: %ld\n", 
       CUTOUT_NROW(mosaic), CUTOUT_NCOL(mosaic));


   // This should agree with printout for the single cutout above
   printf("value at [%ld,%ld]\n", row, col,
       MOSAIC_GET(mosaic, icutout, row, col));
   printf("value at center [%lf,%lf]\n", rowcen, colcen,
       MOSAIC_GET(mosaic, icutout, (int)rowcen, (int)colcen));


   // use the pointer interface; good if you have your own image library
   long ncutout=0, nrow=0, ncol=0;
   double *mpix=meds_get_mosaicp(meds, iobj, &ncutout, &nrow, &ncol);


   // get the distortion matrix, jacobian of the transformation between row,col
   // and ra,dec tangent plane.
   // see definition of struct meds_distort
   const struct meds_distort *dist=meds_get_distortion(meds,iobj,icutout);


   free(pix);pix=NULL;
   free(wpix);wpix=NULL;
   free(mpix);mpix=NULL;

   // free the cutout structures.  They are set to NULL.
   cutout  = meds_cutout_free(cutout);
   wcutout = meds_cutout_free(wcutout);
   scutout = meds_icutout_free(scutout);
   mosaic  = meds_cutout_free(mosaic);
   wmosaic = meds_cutout_free(wmosaic);

   
   // the meds_obj structure contains additional information such as where the
   // cutouts were located in the original source images.  see the struct
   // definition for details, but note that this structure may change in the
   // future.
   //
   // TODO add getters so one does not need to use a meds_obj struct directly

   // get a meds_obj structure for an object

   const struct meds_obj *obj=meds_get_obj(meds, iobj);
   meds_obj_print(obj, stdout);

*/
#ifndef _MEDS_INCLUDE_HEADER_GUARD
#define _MEDS_INCLUDE_HEADER_GUARD

#include <fitsio.h>

//#define MEDS_NCOLUMNS 14

#define MEDS_DEFVAL -9999

// jacobian taking (row,col) -> tangent plane ra,dec (u,v)
struct meds_distort {
    double dudrow;
    double dudcol;
    double dvdrow;
    double dvdcol;
};

// this structure might change; e.g. we might implement 
// not-square cutouts
struct meds_obj {
    long ncutout;    // number of cutouts for this object, including coadd
    long box_size;   // cutout size is box_sizeXbox_size.  Note we might make
                     // cutouts not square in the future so don't rely on this
    long *file_id;   // index into the image_info structure for this cutout
    long *start_row; // zero-offset row in the big mosaic image of all cutouts
                        // for all objects in the MEDS file
    double  *orig_row;  // zero-offset center row in original image
    double  *orig_col;  // zero-offset center col in original image
    long    *orig_start_row; // zero-offset start row in original image
    long    *orig_start_col; // zero-offset start col in original image
    double  *cutout_row;     // zero-offset center row in cutout image
    double  *cutout_col;     // zero-offset center col in cutout image

    // distortion jacobian from pixels (row,col) to tangent plane ra,dec (u,v)
    struct meds_distort *distortion;
};

struct meds_cat {
    long size;
    struct meds_obj *data;
};

struct meds_image_info {
    long size;
    char* image_path;
    char* sky_path;
    char* seg_path;
};
struct meds_info_cat {
    long size;
    struct meds_image_info *data;
};

struct meds_meta {
    char *cat_file;
    char *coadd_file;
    char *coadd_srclist;
    char *cutout_file;

    // optional
    char *coaddcat_file;
    char *medsconf;
    int min_boxsize;
    int max_boxsize;
};
void meds_meta_print(const struct meds_meta *self, FILE *stream);

struct meds {
    char *meds_path;
    fitsfile *fits;
    struct meds_cat *cat;
    struct meds_info_cat *image_info;
    struct meds_meta *meta;
};


// open a meds structure
struct meds *meds_open(const char *filename);

// free the structure.  Returns NULL.  Use like this:
//    m=meds_free(m);
struct meds *meds_free(struct meds *self);

// number of entries in the catalog
long meds_get_size(const struct meds *self);

// the meds cutout filename
const char *meds_get_path(const struct meds *self);

// get the metadata structure with input file paths etc.
const struct meds_meta *meds_get_meta(const struct meds *self);

// get an entry in the catalog
// if index does not exist, NULL is returned
// note you shouldn't generally rely on this structure as it
// may change
const struct meds_obj *meds_get_obj(const struct meds *self, long iobj);

// return obj->ncutout for the indicated object
// if index does not exist, zero is returned
long meds_get_ncutout(const struct meds *self, long iobj);

// get object location in the cutout
// returns false if the object or cutout does not exist
int meds_get_cutout_cen(const struct meds *self,
                        long iobj,
                        long icutout,
                        double *row,
                        double *col);

// you can also work directly with the catalog
// again, remember the meds_obj structure may change
const struct meds_cat *meds_get_cat(const struct meds *self);

// read a single cutout as a simple pointer
// 
// returns the pixel array.  The number of rows and columns are stored
// in nrow,ncol entered by the user.
double *meds_get_cutoutp(const struct meds *self,
                         long iobj,
                         long icutout,
                         long *nrow,
                         long *ncol);

// read a cutout mosaic as a simple pointer
//
// returns the pixel array.  The number of cutouts, rows and columns are stored
// in nrow,ncol entered by the user.
double *meds_get_mosaicp(const struct meds *self,
                         long iobj,
                         long *ncutout,
                         long *nrow,
                         long *ncol);

// same but for weight image cutouts
double *meds_get_weight_cutoutp(const struct meds *self,
                                long iobj,
                                long icutout,
                                long *nrow,
                                long *ncol);

double *meds_get_weight_mosaicp(const struct meds *self,
                                long iobj,
                                long *ncutout,
                                long *nrow,
                                long *ncol);

// same but for seg image cutouts
int *meds_get_seg_cutoutp(const struct meds *self,
                          long iobj,
                          long icutout,
                          long *nrow,
                          long *ncol);

int *meds_get_seg_mosaicp(const struct meds *self,
                          long iobj,
                          long *ncutout,
                          long *nrow,
                          long *ncol);


// get info for the source image of the indicated cutout
const struct meds_image_info *meds_get_source_info(const struct meds *self,
                                                   long iobj,
                                                   long icutout);

// get the file_id for the source image of the indicated cutout.  This
// is used to get source file information
long meds_get_source_file_id(const struct meds *self,
                             long iobj,
                             long icutout);

// get the filename for the source image of the indicated cutout
const char *meds_get_source_path(const struct meds *self,
                                 long iobj,
                                 long icutout);

// get a reference to the distortion for this obj and citout
const struct meds_distort
*meds_get_distortion(const struct meds *self,
                     long iobj,
                     long icutout);

// print tools
void meds_print(const struct meds *self, FILE* stream);
void meds_obj_print(const struct meds_obj *obj, FILE* stream);
void meds_image_info_print(const struct meds_image_info *self, FILE* stream);


//
// simple image class for cutouts/mosaics
//
struct meds_cutout {
    long ncutout;

    long mosaic_size;       // ncutout*nrow*ncol
    long mosaic_nrow;       // ncutout*nrow
    long mosaic_ncol;       // ncol

    long cutout_size;       // nrow*ncol
    long cutout_nrow;       // per cutout
    long cutout_ncol;       // per cutout
    double **rows;
};
struct meds_icutout {
    long ncutout;

    long mosaic_size;       // ncutout*nrow*ncol
    long mosaic_nrow;       // ncutout*nrow
    long mosaic_ncol;       // ncol

    long cutout_size;       // nrow*ncol
    long cutout_nrow;       // per cutout
    long cutout_ncol;       // per cutout
    int **rows;
};


#define CUTOUT_SIZE(im) ((im)->cutout_size)
#define CUTOUT_NROW(im) ((im)->cutout_nrow)
#define CUTOUT_NCOL(im) ((im)->cutout_ncol)

#define CUTOUT_GET(im, row, col)  ( *((im)->rows[(row)] + (col)) )

#define CUTOUT_GET_ROW(im, row) ( (im)->rows[(row) )


#define MOSAIC_NCUTOUT(im) ((im)->ncutout)
#define MOSAIC_SIZE(im) ((im)->mosaic_size)
#define MOSAIC_NROW(im) ((im)->mosaic_nrow)
#define MOSAIC_NCOL(im) ((im)->mosaic_ncol)

#define MOSAIC_GET(im, cutout, row, col)                       \
    ( *((im)->rows[(cutout)*(im)->cutout_nrow + (row)] + (col)) )

#define MOSAIC_GET_ROW(im, cutout, row)                       \
    ((im)->rows[(cutout)*(im)->cutout_nrow + (row)] )

// read a single cutout.  Use CUTOUT_GET or MOSAIC_GET(0, row,col) to
// access pixels
struct meds_cutout *meds_get_cutout(const struct meds *self,
                                    long iobj,
                                    long icutout);

// read a cutout mosaic.  Use MOSAIC_GET to access pixels
struct meds_cutout *meds_get_mosaic(const struct meds *self,
                                    long iobj);

struct meds_cutout *meds_get_weight_cutout(const struct meds *self,
                                           long iobj,
                                           long icutout);
struct meds_cutout *meds_get_weight_mosaic(const struct meds *self,
                                           long iobj);

struct meds_icutout *meds_get_seg_cutout(const struct meds *self,
                                         long iobj,
                                         long icutout);
struct meds_icutout *meds_get_seg_mosaic(const struct meds *self,
                                         long iobj);


// returns NULL, use like this
//   cutout=meds_cutout_free(cutout);
struct meds_cutout *meds_cutout_free(struct meds_cutout *self);
struct meds_icutout *meds_icutout_free(struct meds_icutout *self);

#endif
