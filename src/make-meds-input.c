/*
   Read a fits Source Extractor fits file and write out an ascii
   file for input to the make-cutouts program.

   Input minsize and maxsize will be rounded up to nearest 2^N or
   3*2^N

       ra dec row col box_size

   The box sizes will be size 2^N or 3*2^N

   The input fits catalog should have columns 
       alphamodel_j2000, deltamodel_j2000
       x_image y_image
       xmin_image xmax_image ymin_image ymax_image
       flux_radius
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fitsio.h>

#define CAT_HDU 2
#define SIGMA_FAC 5

#define FWHM_FAC 2.3548200450309493

static int SIZES[]={2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,
                    384,512,768,1024,1536,2048,3072,4096,6144};
static long NSIZES=sizeof(SIZES)/sizeof(int);

struct obj {
    double ra;
    double dec;
    double row;
    double col;
    int rowmin;
    int rowmax;
    int colmin;
    int colmax;

    double flux_radius;
    double ellipticity;

    // output
    int box_size;
};

struct colnums {
    int ra;
    int dec;
    int row;
    int col;
    int rowmin;
    int rowmax;
    int colmin;
    int colmax;
    int flux_radius;
    //int ellipticity;
    int a_world;
    int b_world;
};

struct cat {
    long size;
    struct obj *data;
};

struct cat *cat_new(long size)
{
    struct cat *self=calloc(1,sizeof(struct cat));
    if (!self) {
        fprintf(stderr,"could not allocate struct cat\n");
        exit(1);
    }
    self->size=size;

    self->data=calloc(size, sizeof(struct obj));
    if (!self->data) {
        fprintf(stderr,"could not %ld struct obj\n",size);
        exit(1);
    }
    return self;
}
struct cat *cat_free(struct cat *self)
{
    if (self) {
        if (self->data) {
            free(self->data);
        }
        free(self);
    }
    return NULL;
}

int get_colnum(fitsfile *fits, const char *colname)
{
    int status=0;
    int colnum=0;
    if (fits_get_colnum(fits, 0, (char*)colname, &colnum, &status)) {
        fits_report_error(stderr,status);
        exit(1);
    }
    return colnum;
}

void get_colnums(fitsfile *fits, struct colnums *colnums)
{
    colnums->ra=get_colnum(fits, "ALPHAMODEL_J2000");
    colnums->dec=get_colnum(fits, "DELTAMODEL_J2000");
    colnums->row=get_colnum(fits, "Y_IMAGE");
    colnums->col=get_colnum(fits, "X_IMAGE");

    colnums->rowmin=get_colnum(fits, "YMIN_IMAGE");
    colnums->rowmax=get_colnum(fits, "YMAX_IMAGE");
    colnums->colmin=get_colnum(fits, "XMIN_IMAGE");
    colnums->colmax=get_colnum(fits, "XMAX_IMAGE");

    colnums->flux_radius=get_colnum(fits, "FLUX_RADIUS");
    //colnums->ellipticity=get_colnum(fits, "ELLIPTICITY");
    colnums->a_world=get_colnum(fits, "A_WORLD");
    colnums->b_world=get_colnum(fits, "B_WORLD");
}


long get_nrows(fitsfile *fits)
{
    int fitserr=0;
    long nrows=0;
    fits_get_num_rows(fits, &nrows, &fitserr);
    if (fitserr!=0) {
        fits_report_error(stderr,fitserr);
        exit(1);
    }

    return nrows;
}

fitsfile *open_fits(const char *filename)
{
    fitsfile *fits=NULL;
    int fitserr=0;

    fits_open_file(&fits,filename,READWRITE,&fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        exit(1);
    }
    return fits;
}
void fits_goto_hdu(fitsfile *fits, int hdu)
{
    int fitserr=0;
    fits_movabs_hdu(fits,hdu,0,&fitserr);
    if (fitserr != 0) {
        fits_report_error(stderr,fitserr);
        exit(1);
    }
}

double fits_load_col_dbl(fitsfile *fits, int colnum, LONGLONG row)
{
    int status=0;
    int nullval=0;
    LONGLONG firstelem=1;
    LONGLONG nelem=1;
    double val=0;
    if (fits_read_col_dbl(fits, colnum, row, firstelem, nelem,
                          nullval, &val, NULL, &status)) {
        fits_report_error(stderr,status);
        exit(1);
    }
    return val;
}
double fits_load_col_int(fitsfile *fits, int colnum, LONGLONG row)
{
    int status=0;
    int nullval=0;
    LONGLONG firstelem=1;
    LONGLONG nelem=1;
    int val=0;
    if (fits_read_col_int(fits, colnum, row, firstelem, nelem,
                          nullval, &val, NULL, &status)) {
        fits_report_error(stderr,status);
        exit(1);
    }
    return val;
}


void load_rows(fitsfile *fits, struct cat *cat)
{
    long nrows=cat->size;
    fprintf(stderr,"loading %ld rows\n", nrows);

    struct colnums colnums={0};
    get_colnums(fits, &colnums);

    struct obj *obj=cat->data;
    for (long row=1; row<=nrows; row++) {

        obj->ra=fits_load_col_dbl(fits, colnums.ra, row);
        obj->dec=fits_load_col_dbl(fits, colnums.dec, row);
        obj->row=fits_load_col_dbl(fits, colnums.row, row);
        obj->col=fits_load_col_dbl(fits, colnums.col, row);

        obj->rowmin=fits_load_col_int(fits, colnums.rowmin, row);
        obj->rowmax=fits_load_col_int(fits, colnums.rowmax, row);
        obj->colmin=fits_load_col_int(fits, colnums.colmin, row);
        obj->colmax=fits_load_col_int(fits, colnums.colmax, row);

        obj->flux_radius=fits_load_col_dbl(fits, colnums.flux_radius, row);

        double a_world=fits_load_col_dbl(fits, colnums.a_world, row);
        double b_world=fits_load_col_dbl(fits, colnums.b_world, row);

        obj->ellipticity=1.0 - b_world/a_world;

        obj++;
    }
}
struct cat *read_cat(const char *fname)
{
    int fitserr=0;
    fitsfile *fits=open_fits(fname);
    fits_goto_hdu(fits, CAT_HDU);

    long nrows=get_nrows(fits);
    struct cat *cat=cat_new(nrows);

    load_rows(fits, cat);

    fits_close_file(fits, &fitserr);

    return cat;
}

// write the full object
void obj_write(struct obj *self, FILE *stream)
{
    fprintf(stream,
      "%.16g %.16g %.16g %.16g %d %d %d %d %.16g %.16g %d\n",
      self->ra, self->dec, self->row, self->col,
      self->rowmin, self->rowmax, self->colmin, self->colmax,
      self->flux_radius, self->ellipticity, self->box_size);
}


// write the full catalog
void cat_write(struct cat *self, FILE *stream)
{
    struct obj *obj=self->data;
    for (long i=0; i<self->size; i++) {
        obj_write(obj, stream);
        obj++;
    }
}

// this is the actual input to make-cutouts
void cat_write_meds(struct cat *self)
{
    struct obj *obj=self->data;
    for (long i=0; i<self->size; i++) {
        printf("%.16g %.16g %.16g %.16g %d\n",
               obj->ra, obj->dec, obj->row, obj->col,
               obj->box_size);
        obj++;
    }
}


// round a size up to the nearest 2^N or 3*2^N
int fft_round_size(int size)
{
    int newsize=0;
    int found=0;
    for (long i=0; i<NSIZES; i++) {
        newsize=SIZES[i];
        if (newsize >= size) {
            found=1;
            break;
        }
    }
    if (!found) {
        fprintf(stderr,"Size %d is out of bounds [0,%d], setting to %d\n",
                size,newsize,newsize);
    }
    return newsize;
}

// convert 
int get_sigma_size(double flux_radius, double ellipticity)
{
    double sigma=flux_radius*2./FWHM_FAC;
    double drad = sigma*SIGMA_FAC;

    drad *= (1. + ellipticity);

    drad = ceil(drad);

    // box size is twice the radius
    int box_size = 2*( (int)drad );

    return box_size;
}
int get_box_size(struct obj *self, int minsize, int maxsize)
{
    int rowsize=self->rowmax-self->rowmin+1;
    int colsize=self->colmax-self->colmin+1;

    int sigma_size=get_sigma_size(self->flux_radius, self->ellipticity);

    int box_size=rowsize > colsize ? rowsize : colsize;

    box_size = sigma_size > box_size ? sigma_size: box_size;

    box_size = box_size < minsize ? minsize : box_size;
    box_size = box_size > maxsize ? maxsize : box_size;

    box_size = fft_round_size(box_size);
    return box_size;
}

void set_box_sizes(struct cat *self, int minsize, int maxsize)
{
    struct obj *obj=self->data;
    for (long i=0; i<self->size; i++) {
        obj->box_size = get_box_size(obj, minsize, maxsize);
        obj++;
    }
}
int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr,"usage: make-meds-input fitsfile minsize maxsize\n");
        fprintf(stderr,"  results go to stdout\n");
        exit(1);
    }

    const char *infile=argv[1];
    int minsize=atoi(argv[2]);
    int maxsize=atoi(argv[3]);

    minsize=fft_round_size(minsize);
    maxsize=fft_round_size(maxsize);

    fprintf(stderr,"minsize: %d\n", minsize);
    fprintf(stderr,"maxsize: %d\n", maxsize);

    struct cat *cat=read_cat(infile);

    set_box_sizes(cat, minsize, maxsize);

    cat_write_meds(cat);
    cat=cat_free(cat);
}
