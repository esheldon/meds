/*


   Two error styles are used here.  cfitsio uses a true return value to mean
   error, my routines use false to mean an error.  The user should only see
   the "false is error" in the public api.

   Everywhere zero-offset is used except when calling the cfitsio routines,
   where they are converted to 1-offset.

*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fitsio.h>

#include "meds.h"

//
// some generic static methods and variables
//

static char *MEDS_COLNAMES[]={
    "ncutout", "box_size",
    "file_id", "start_row", 
    "orig_row", "orig_col","orig_start_row","orig_start_col",
    "cutout_row","cutout_col",
    "dudrow","dudcol","dvdrow","dvdcol"};
static int MEDS_NCOLUMNS = sizeof(MEDS_COLNAMES)/sizeof(char*);

static long *alloc_longs(long n) {
    long *ptr=malloc(n*sizeof(long));
    if (!ptr) {
        fprintf(stderr,"could not %ld longs\n", n);
        exit(1);
    }
    for (long i=0; i<n; i++) {
        ptr[i]=MEDS_DEFVAL;
    }
    return ptr;
}
static double *alloc_doubles(long n) {
    double *ptr=malloc(n*sizeof(double));
    if (!ptr) {
        fprintf(stderr,"could not %ld doubles\n", n);
        exit(1);
    }
    for (long i=0; i<n; i++) {
        ptr[i]=MEDS_DEFVAL;
    }
    return ptr;
}
static int *alloc_ints(long n) {
    int *ptr=malloc(n*sizeof(int));
    if (!ptr) {
        fprintf(stderr,"could not %ld ints\n", n);
        exit(1);
    }
    for (long i=0; i<n; i++) {
        ptr[i]=MEDS_DEFVAL;
    }
    return ptr;
}


static void print_doubles(const double* vals, long n, const char *name,
                          FILE* stream)
{
    fprintf(stream,"%-14s : ", name);
    for (long i=0; i<n; i++) {
        fprintf(stream," %lf", vals[i]);
    }
    fprintf(stream,"\n");
}
static void print_longs(const long *vals, long n, const char *name,
                        FILE* stream)
{
    fprintf(stream,"%-14s : ", name);
    for (long i=0; i<n; i++) {
        fprintf(stream," %ld", vals[i]);
    }
    fprintf(stream,"\n");
}

static fitsfile *fits_open(const char *filename)
{
    int status=0;
    fitsfile *fits=NULL;
    if (fits_open_file(&fits, filename, READONLY, &status)) {
        fits_report_error(stderr,status);
        return NULL;
    }
    return fits;
}
static fitsfile *fits_close(fitsfile *fits)
{
    int status=0;
    if (fits_close_file(fits, &status)) {
        fits_report_error(stderr,status);
    }
    return NULL;
}

static long get_nrow(fitsfile *fits)
{
    int status=0;
    long nrow=0;
    if (fits_get_num_rows(fits, &nrow, &status)) {
        fits_report_error(stderr,status);
        return 0;
    }
    return nrow;
}

static int get_colnum(fitsfile *fits, const char *colname)
{
    int status=0;
    int colnum=0;
    if (fits_get_colnum(fits, 0, (char*)colname, &colnum, &status)) {
        fits_report_error(stderr,status);
        return MEDS_DEFVAL;
    }
    return colnum;
}
static int get_colnums(fitsfile *fits, int *colnums)
{
    for (int i=0; i<MEDS_NCOLUMNS; i++) {
        int colnum=get_colnum(fits, MEDS_COLNAMES[i]);
        if (colnum < 0) {
            return 0;
        }
        colnums[i] = colnum;
    }
    return 1;
}

// input rows are zero offset, get converted to 1 offset
static int fits_load_col_dbl(fitsfile *fits, 
                             int colnum,
                             LONGLONG row,
                             LONGLONG nelem,
                             double *data)
{
    int status=0;
    int nullval=0;
    LONGLONG firstelem=1;
    if (fits_read_col_dbl(fits, colnum, 1+row, firstelem, nelem,
                          nullval, data, NULL, &status)) {
        fits_report_error(stderr,status);
        return 0;
    }
    return 1;
}
static int fits_load_col_lng(fitsfile *fits, 
                             int colnum,
                             LONGLONG row,
                             LONGLONG nelem,
                             long *data)
{
    int status=0;
    int nullval=0;
    LONGLONG firstelem=1;
    if (fits_read_col_lng(fits, colnum, 1+row, firstelem, nelem,
                          nullval, data, NULL, &status)) {
        fits_report_error(stderr,status);
        return 0;
    }
    return 1;
}

static int fits_load_sub_dbl(fitsfile *fits, 
                             LONGLONG row,
                             LONGLONG nelem,
                             double *data)
{
    int status=0;
    LONGLONG firstpixels[1];

    // note col comes first
    firstpixels[0] = 1+row;

    if (fits_read_pixll(fits, TDOUBLE, firstpixels, nelem,
                        0, data, NULL, &status)) {
        fits_report_error(stderr,status);
        return 0;
    }
    return 1;
}
static int fits_load_sub_int(fitsfile *fits, 
                             LONGLONG row,
                             LONGLONG nelem,
                             int *data)
{
    int status=0;
    LONGLONG firstpixels[1];

    // note col comes first
    firstpixels[0] = 1+row;

    if (fits_read_pixll(fits, TINT, firstpixels, nelem,
                        0, data, NULL, &status)) {
        fits_report_error(stderr,status);
        return 0;
    }
    return 1;
}


// size of a fixed-length array column (per row)
static long get_array_col_size(fitsfile *fits, const char *colname)
{
    int status=0;
    int maxdim=2;
    int naxis = 0;
    long naxes[2]={0};

    int colnum=0;
    if (fits_get_colnum(fits, 0, (char*)colname, &colnum, &status)) {
        fits_report_error(stderr,status);
        return 0;
    }
    if (fits_read_tdim(fits, colnum, maxdim, &naxis,
                       naxes, &status)) {
        fits_report_error(stderr,status);
        return 0;
    }
    return naxes[0];
}

//
// struct meds_distort
//

static struct meds_distort *alloc_distortions(long n)
{
    struct meds_distort *ds=malloc(n*sizeof(struct meds_distort));
    if (!ds) {
        fprintf(stderr,"Could not allocate %ld distortion structs\n", n);
        exit(1);
    }
    for (long j=0; j<n; j++) {
        ds[j].dudrow=MEDS_DEFVAL;
        ds[j].dudcol=MEDS_DEFVAL;
        ds[j].dvdrow=MEDS_DEFVAL;
        ds[j].dvdcol=MEDS_DEFVAL;
    }

    return ds;
}

void meds_distort_print(struct meds_distort *self, FILE *stream)
{
    fprintf(stream,"%g %g %g %g\n",
            self->dudrow,
            self->dudcol,
            self->dvdrow,
            self->dvdcol);
}
static void print_distortions(const struct meds_distort *self,
                              long n,
                              FILE *stream)
{

    fprintf(stream,"%-14s : ", "dudrow");
    for (long i=0; i<n; i++) {
        fprintf(stream," %lf", self[i].dudrow);
    }
    fprintf(stream,"\n");

    fprintf(stream,"%-14s : ", "dudcol");
    for (long i=0; i<n; i++) {
        fprintf(stream," %lf", self[i].dudcol);
    }
    fprintf(stream,"\n");

    fprintf(stream,"%-14s : ", "dvdrow");
    for (long i=0; i<n; i++) {
        fprintf(stream," %lf", self[i].dvdrow);
    }
    fprintf(stream,"\n");

    fprintf(stream,"%-14s : ", "dvdcol");
    for (long i=0; i<n; i++) {
        fprintf(stream," %lf", self[i].dvdcol);
    }
    fprintf(stream,"\n");
}


//
// struct meds_obj
//

static void alloc_fields(struct meds_obj *self, long ncutout)
{
    self->file_id        = alloc_longs(ncutout);
    self->start_row      = alloc_longs(ncutout);
    self->orig_row       = alloc_doubles(ncutout);
    self->orig_col       = alloc_doubles(ncutout);

    self->orig_start_row = alloc_longs(ncutout);
    self->orig_start_col = alloc_longs(ncutout);

    self->cutout_row     = alloc_doubles(ncutout);
    self->cutout_col     = alloc_doubles(ncutout);

    self->distortion     = alloc_distortions(ncutout);
}
static void free_fields(struct meds_obj *self)
{
    if (self) {
        free(self->file_id);
        free(self->start_row);
        free(self->orig_row);
        free(self->orig_col);

        free(self->orig_start_row);
        free(self->orig_start_col);

        free(self->cutout_row);
        free(self->cutout_col);

        free(self->distortion);
    }
}

void meds_obj_print(const struct meds_obj *self, FILE* stream)
{
    fprintf(stream,"%-14s : %ld\n", "ncutout", self->ncutout);
    fprintf(stream,"%-14s : %ld\n", "box_size",   self->box_size);
    print_longs(self->file_id, self->ncutout, "file_id", stream);
    print_longs(self->start_row, self->ncutout, "start_row", stream);
    print_doubles(self->orig_row, self->ncutout, "orig_row", stream);
    print_doubles(self->orig_col, self->ncutout, "orig_col", stream);
    print_longs(self->orig_start_row, self->ncutout, "orig_start_row", stream);
    print_longs(self->orig_start_col, self->ncutout, "orig_start_col", stream);
    print_doubles(self->cutout_row, self->ncutout, "cutout_row", stream);
    print_doubles(self->cutout_col, self->ncutout, "cutout_col", stream);

    print_distortions(self->distortion, self->ncutout, stream);
}

//
// struct meds_cat
//

static struct meds_cat *meds_cat_new(long size, long ncutout_max)
{
    struct meds_cat *self=calloc(1, sizeof(struct meds_cat));
    if (!self) {
        fprintf(stderr,"could not allocate struct meds_cat\n");
        exit(1);
    }

    self->size=size;
    self->data=calloc(size, sizeof(struct meds_obj));
    if (!self->data) {
        fprintf(stderr,"could not allocate %ld meds_obj structures\n", size);
        exit(1);
    }

    struct meds_obj *obj=self->data;
    for (long i=0; i<size; i++) {
        alloc_fields(obj, ncutout_max);
        obj->ncutout=ncutout_max;
        obj++;
    }

    return self;
}

static struct meds_cat *meds_cat_free(struct meds_cat *self)
{
    if (self) {
        struct meds_obj *obj=self->data;
        for (long i=0; i<self->size; i++) {
            free_fields(obj);
            obj++;
        }
        free(self->data);
        free(self);
    }
    return NULL;
}

// this is kind of a monster
static int load_table(struct meds_cat *self, fitsfile *fits)
{
    int colnums[MEDS_NCOLUMNS];
    if (!get_colnums(fits, colnums)) {
        return 0;
    }

    int ndist=20;
    double *dudrow=alloc_doubles(ndist);
    double *dudcol=alloc_doubles(ndist);
    double *dvdrow=alloc_doubles(ndist);
    double *dvdcol=alloc_doubles(ndist);

    struct meds_obj *obj=self->data;
    for (long row=0; row < self->size; row++) {

    
        if (!fits_load_col_lng(fits, colnums[0], row, 1, &obj->ncutout))
            return 0;
        if (!fits_load_col_lng(fits, colnums[1], row, 1, &obj->box_size))
            return 0;
        if (!fits_load_col_lng(fits, colnums[2], row, obj->ncutout, obj->file_id))
            return 0;
        if (!fits_load_col_lng(fits, colnums[3], row, obj->ncutout, obj->start_row))
            return 0;
        if (!fits_load_col_dbl(fits, colnums[4], row, obj->ncutout, obj->orig_row))
            return 0;
        if (!fits_load_col_dbl(fits, colnums[5], row, obj->ncutout, obj->orig_col))
            return 0;
        if (!fits_load_col_lng(fits, colnums[6], row, obj->ncutout, obj->orig_start_row))
            return 0;
        if (!fits_load_col_lng(fits, colnums[7], row, obj->ncutout, obj->orig_start_col))
            return 0;
        if (!fits_load_col_dbl(fits, colnums[8], row, obj->ncutout, obj->cutout_row))
            return 0;
        if (!fits_load_col_dbl(fits, colnums[9], row, obj->ncutout, obj->cutout_col))
            return 0;

        if (obj->ncutout > ndist) {
            ndist=obj->ncutout;
            dudrow=realloc(dudrow,ndist*sizeof(double));
            dudcol=realloc(dudcol,ndist*sizeof(double));
            dvdrow=realloc(dvdrow,ndist*sizeof(double));
            dvdcol=realloc(dvdcol,ndist*sizeof(double));
        }
        if (!fits_load_col_dbl(fits, colnums[10], row, obj->ncutout, dudrow))
            return 0;
        if (!fits_load_col_dbl(fits, colnums[11], row, obj->ncutout, dudcol))
            return 0;
        if (!fits_load_col_dbl(fits, colnums[12], row, obj->ncutout, dvdrow))
            return 0;
        if (!fits_load_col_dbl(fits, colnums[13], row, obj->ncutout, dvdcol))
            return 0;

        for (long j=0; j<obj->ncutout; j++) {
            obj->distortion[j].dudrow=dudrow[j];
            obj->distortion[j].dudcol=dudcol[j];
            obj->distortion[j].dvdrow=dvdrow[j];
            obj->distortion[j].dvdcol=dvdcol[j];
        }
        obj++;
    }

    free(dudrow);
    free(dudcol);
    free(dvdrow);
    free(dvdcol);
    return 1;
}

static struct meds_cat *read_object_table(fitsfile *fits)
{
    int status=0;
    if (fits_movnam_hdu(fits, BINARY_TBL, "object_data", 0, &status)) {
        fits_report_error(stderr,status);
        return NULL;
    }

    long nrow=get_nrow(fits);
    if (!nrow) {
        return NULL;
    }
    //printf("found %ld objects in main table\n", nrow);

    long ncutout_max = get_array_col_size(fits, "file_id");
    if (!nrow) {
        return NULL;
    }
    //printf("ncutout max: %ld\n", ncutout_max);

    struct meds_cat *self=meds_cat_new(nrow, ncutout_max);

    if (!load_table(self, fits)) {
        self=meds_cat_free(self);
        return NULL;
    }
    return self;
}


//
// struct meds_image_info
//

void meds_image_info_print(const struct meds_image_info *self, FILE* stream)
{
    fprintf(stream,"image_path: %s\n", self->image_path);
    fprintf(stream,"sky_path:   %s\n", self->sky_path);
    fprintf(stream,"seg_path:   %s\n", self->seg_path);
}


//
// struct meds_info_cat
//

static struct meds_info_cat *meds_info_cat_new(long size, long namelen_max)
{
    struct meds_info_cat *self=calloc(1, sizeof(struct meds_info_cat));
    if (!self) {
        fprintf(stderr,"could not allocate struct meds_info_cat\n");
        exit(1);
    }

    self->size=size;
    self->data=calloc(size, sizeof(struct meds_image_info));
    if (!self->data) {
        fprintf(stderr,"could not allocate %ld meds_image_info structures\n", size);
        exit(1);
    }

    struct meds_image_info *info=self->data;
    for (long i=0; i<size; i++) {
        info->image_path = calloc(namelen_max, sizeof(char));
        info->sky_path = calloc(namelen_max, sizeof(char));
        info->seg_path = calloc(namelen_max, sizeof(char));

        if (!info->image_path || !info->sky_path) {
            fprintf(stderr,"failed to allocate filename of size %ld\n", namelen_max);
            exit(1);
        }
        info->size=namelen_max;
        info++;
    }

    return self;
}

static struct meds_info_cat *meds_info_cat_free(struct meds_info_cat *self)
{
    if (self) {
        struct meds_image_info *info=self->data;
        for (long i=0; i<self->size; i++) {
            free(info->image_path);
            free(info->sky_path);
            free(info->seg_path);
            info++;
        }
        free(self->data);
        free(self);
    }
    return NULL;
}


static int load_image_info(struct meds_info_cat *self, fitsfile *fits)
{
    int status=0;
    int colnum=1;
    char* nulstr=" ";
    LONGLONG firstelem=1;

    struct meds_image_info *info=self->data;
    for (long i=0; i<self->size; i++) {
        long row = i+1;

        colnum=1;
        if (fits_read_col_str(fits, colnum, row, firstelem, 1,
                              nulstr, &info->image_path, NULL, &status)) {
            fits_report_error(stderr,status);
            return 0;
        }
        colnum=2;
        if (fits_read_col_str(fits, colnum, row, firstelem, 1,
                              nulstr, &info->sky_path, NULL, &status)) {
            fits_report_error(stderr,status);
            return 0;
        }
        colnum=3;
        if (fits_read_col_str(fits, colnum, row, firstelem, 1,
                              nulstr, &info->seg_path, NULL, &status)) {
            fits_report_error(stderr,status);
            return 0;
        }


        info++;
    }

    return 1;
}

static long get_namelen_max(fitsfile *fits)
{
    long src_name_len = get_array_col_size(fits,"image_path");
    long sky_name_len = get_array_col_size(fits,"sky_path");
    long seg_name_len = get_array_col_size(fits,"seg_path");
    long namelen_max=src_name_len;

    if (sky_name_len > src_name_len) {
        namelen_max=sky_name_len;
    }
    if (seg_name_len > src_name_len) {
        namelen_max=seg_name_len;
    }

    return namelen_max;
}
static struct meds_info_cat *read_image_info(fitsfile *fits)
{
    int status=0;
    if (fits_movnam_hdu(fits, BINARY_TBL, "image_info", 0, &status)) {
        fits_report_error(stderr,status);
        return NULL;
    }

    long nrow=get_nrow(fits);
    if (!nrow) {
        return NULL;
    }

    long namelen_max=get_namelen_max(fits);
    struct meds_info_cat *self=meds_info_cat_new(nrow, namelen_max);


    if (!load_image_info(self, fits)) {
        self=meds_info_cat_free(self);
        return NULL;
    }
    return self;
}

//
// struct meds_meta
//

static char *fits_read_string_byname(fitsfile *fits, const char* colname,
                                     long row, int *status)
{

    int colnum=get_colnum(fits,colname);
    if (colnum==MEDS_DEFVAL) {
        (*status)=1;
        return NULL;
    }

    long len = get_array_col_size(fits,colname);
    char* nulstr=" ";
    LONGLONG firstelem=1;

    char *data=calloc(len, sizeof(char));
    if (fits_read_col_str(fits, colnum, row, firstelem, 1,
                nulstr, &data, NULL, status)) {
        fits_report_error(stderr,(*status));
        return 0;
    }
    return data;
}
static double fits_read_double_byname(fitsfile *fits, const char* colname, 
                                      long row, int *status)
{

    int colnum=get_colnum(fits,colname);
    if (colnum==MEDS_DEFVAL) {
        (*status)=1;
        return -9999;
    }

    double nulval=0;
    LONGLONG firstelem=1;

    double data=-9999;
    if (fits_read_col_dbl(fits, colnum, row, firstelem, 1,
                nulval, &data, NULL, status)) {
        fits_report_error(stderr,(*status));
        return -9999;
    }
    return data;
}


static struct meds_meta *read_meta(fitsfile *fits)
{
    int status=0;
    if (fits_movnam_hdu(fits, BINARY_TBL, "metadata", 0, &status)) {
        fits_report_error(stderr,status);
        return NULL;
    }

    struct meds_meta *self=calloc(1,sizeof(struct meds_meta));
    if (!self) {
        fprintf(stderr,"could not allocate struct meds_meta\n");
    }
    self->min_boxsize=MEDS_DEFVAL;
    self->max_boxsize=MEDS_DEFVAL;

    long row=1;
    // these will always exist

    self->cat_file = 
        fits_read_string_byname(fits, "cat_file", row,&status);
    self->coadd_file = 
        fits_read_string_byname(fits, "coadd_file", row,&status);
    self->coadd_srclist = 
        fits_read_string_byname(fits, "coadd_srclist", row, &status);
    self->cutout_file = 
        fits_read_string_byname(fits, "cutout_file", row, &status);

    // might not be present, if not will be null
    self->coaddcat_file = 
        fits_read_string_byname(fits, "coaddcat_file", row, &status);
    self->medsconf = 
        fits_read_string_byname(fits, "medsconf", row, &status);

    char *minbox_str = 
        fits_read_string_byname(fits, "min_boxsize", row, &status);
    if (status==0) {
        self->min_boxsize=atoi(minbox_str);
    }
    char *maxbox_str = 
        fits_read_string_byname(fits, "max_boxsize", row, &status);
    if (status==0) {
        self->max_boxsize=atoi(maxbox_str);
    }

    self->magzp_ref = 
        fits_read_double_byname(fits, "magzp_ref", row, &status);

    free(minbox_str);
    free(maxbox_str);

    return self;
}

void meds_meta_print(const struct meds_meta *self, FILE *stream)
{
    fprintf(stderr,"    cat_file:      %s\n", self->cat_file);
    fprintf(stderr,"    coadd_file:    %s\n", self->coadd_file);
    fprintf(stderr,"    coadd_srclist: %s\n", self->coadd_srclist);
    fprintf(stderr,"    cutout_file:   %s\n", self->cutout_file);

    if (self->coaddcat_file)
        fprintf(stderr,"    coaddcat_file: %s\n", self->coaddcat_file);
    if (self->medsconf)
        fprintf(stderr,"    medsconf:      %s\n", self->medsconf);
    if (self->min_boxsize > 0)
        fprintf(stderr,"    min_boxsize:   %d\n", self->min_boxsize);
    if (self->max_boxsize > 0)
        fprintf(stderr,"    max_boxsize:   %d\n", self->max_boxsize);
    if (self->magzp_ref > 0)
        fprintf(stderr,"    magzp_ref:     %lf\n", self->magzp_ref);
}

static struct meds_meta *meds_meta_free(struct meds_meta *self)
{
    if (self) {
        free(self->cat_file);
        free(self->coadd_file);
        free(self->coaddcat_file);
        free(self->coadd_srclist);
        free(self->cutout_file);
        free(self->medsconf);
        free(self);
        self=NULL;
    }
    return self;
}

//
// struct meds_cutout
//

struct meds_cutout *meds_cutout_free(struct meds_cutout *self)
{
    if (self) {
        if (self->rows) {
            if (self->rows[0]) {
                free(self->rows[0]);
            }
            self->rows[0]=NULL;

            free(self->rows);
            self->rows=NULL;
        }
        free(self);
        self=NULL;
    }
    return self;
}
struct meds_icutout *meds_icutout_free(struct meds_icutout *self)
{
    if (self) {
        if (self->rows) {
            if (self->rows[0]) {
                free(self->rows[0]);
            }
            self->rows[0]=NULL;

            free(self->rows);
            self->rows=NULL;
        }
        free(self);
        self=NULL;
    }
    return self;
}


static struct meds_cutout *cutout_from_ptr(double *ptr,
                                           long ncutout,
                                           long nrow,   // per cutout
                                           long ncol)   // per cutout
{
    struct meds_cutout *self=calloc(1, sizeof(struct meds_cutout));
    if (!self) {
        fprintf(stderr,"failed to allocate struct meds_cutout\n");
        exit(1);
    }

    self->ncutout=ncutout;

    self->mosaic_size = ncutout*nrow*ncol;
    self->mosaic_nrow = ncutout*nrow;
    self->mosaic_ncol = ncol;

    self->cutout_size = nrow*ncol;
    self->cutout_nrow=nrow;
    self->cutout_ncol=ncol;

    self->rows = calloc(self->mosaic_nrow,sizeof(double *));
    if (!self->rows) {
        fprintf(stderr,"could not allocate %ld image rows\n", self->mosaic_nrow);
        exit(1);
    }

    self->rows[0] = ptr;
    for(long i = 1; i < self->mosaic_nrow; i++) {
        self->rows[i] = self->rows[i-1] + self->cutout_ncol;
    }

    return self;
}
static struct meds_icutout *icutout_from_ptr(int *ptr,
                                             long ncutout,
                                             long nrow,   // per cutout
                                             long ncol)   // per cutout
{
    struct meds_icutout *self=calloc(1, sizeof(struct meds_icutout));
    if (!self) {
        fprintf(stderr,"failed to allocate struct meds_icutout\n");
        exit(1);
    }

    self->ncutout=ncutout;

    self->mosaic_size = ncutout*nrow*ncol;
    self->mosaic_nrow = ncutout*nrow;
    self->mosaic_ncol = ncol;

    self->cutout_size = nrow*ncol;
    self->cutout_nrow=nrow;
    self->cutout_ncol=ncol;

    self->rows = calloc(self->mosaic_nrow,sizeof(int *));
    if (!self->rows) {
        fprintf(stderr,"could not allocate %ld image rows\n", self->mosaic_nrow);
        exit(1);
    }

    self->rows[0] = ptr;
    for(long i = 1; i < self->mosaic_nrow; i++) {
        self->rows[i] = self->rows[i-1] + self->cutout_ncol;
    }

    return self;
}




//
// struct meds
//

struct meds *meds_open(const char *filename)
{
    int err=0;
    struct meds_cat *cat=NULL;
    struct meds_info_cat *image_info=NULL;
    struct meds_meta *meta=NULL;
    struct meds *self=NULL;

    fitsfile *fits=fits_open(filename);
    if (!fits) {
        return NULL;
    }

    cat=read_object_table(fits);
    if (!cat) {
        err=1;
        goto _meds_open_bail;
    }

    image_info=read_image_info(fits);
    if (!image_info) {
        err=1;
        goto _meds_open_bail;
    }
    meta=read_meta(fits);
    if (!meta) {
        err=1;
        goto _meds_open_bail;
    }

    self=calloc(1, sizeof(struct meds));
    if (!self) {
        fprintf(stderr,"could not allocate struct meds\n");
        exit(1);
    }

    self->meds_path=strdup(filename);
    self->fits=fits;
    self->cat=cat;
    self->image_info=image_info;
    self->meta=meta;

_meds_open_bail:
    if (err) {
        // if we get here, fits was opened but self was not allocated, no need
        // to free these functions are noop for NULL
        fits=fits_close(fits);
        cat=meds_cat_free(cat);
        image_info=meds_info_cat_free(image_info);
        meta=meds_meta_free(meta);
    }
    return self;
}

struct meds *meds_free(struct meds *self)
{
    if (self) {
        free(self->meds_path);
        self->cat=meds_cat_free(self->cat);
        self->image_info=meds_info_cat_free(self->image_info);
        self->meta=meds_meta_free(self->meta);
        self->fits=fits_close(self->fits);
        free(self);
    }
    return NULL;
}


void meds_print(const struct meds *self,FILE* stream)
{
    if (self) {
        fprintf(stream,"MEDS structure\n");
        fprintf(stream,"    meds path:     %s\n", self->meds_path);
        fprintf(stream,"    nobj:          %ld\n", self->cat->size);
        fprintf(stream,"    source images: %ld\n", self->image_info->size);
        fprintf(stream,"metadata\n");
        meds_meta_print(self->meta, stream);
    }
}

static int check_iobj(const struct meds_cat *self, long iobj)
{
    if (iobj < 0 || iobj >= self->size) {
        fprintf(stderr,"iobj %ld out of range [0,%ld)\n", 
                iobj, self->size);
        return 0;
    } else {
        return 1;
    }
}

long meds_get_size(const struct meds *self)
{
    return self->cat->size;
}
const char *meds_get_path(const struct meds *self)
{
    return self->meds_path;
}
const struct meds_meta *meds_get_meta(const struct meds *self)
{
    return self->meta;
}

const struct meds_obj *meds_get_obj(const struct meds *self, long iobj)
{
    const struct meds_cat *cat=self->cat;
    if (!check_iobj(cat, iobj)) {
        return NULL;
    }
    return &cat->data[iobj];
}
long meds_get_ncutout(const struct meds *self, long iobj)
{
    const struct meds_obj *obj=meds_get_obj(self, iobj);
    if (!obj) {
        return 0;
    }
    return obj->ncutout;
}

const struct meds_cat *meds_get_cat(const struct meds *self)
{
    return self->cat;
}

/*
static int check_icutout(const struct meds_obj *self, long icutout)
{
    long ncutout=self->ncutout;
    if (icutout < 0 || icutout >= ncutout) {
        fprintf(stderr,
                "icutout %ld out of range [0,%ld)\n",
                icutout, ncutout);
        return 0;
    } else {
        return 1;
    }
}
*/

static const struct meds_obj *check_iobj_icutout(const struct meds *self, 
                                                 long iobj,
                                                 long icutout)
{
    const struct meds_obj *obj=meds_get_obj(self, iobj);
    if (!obj) {
        return NULL;
    } else {
        long ncutout=obj->ncutout;
        if (icutout < 0 || icutout >= ncutout) {
            fprintf(stderr,
                    "icutout %ld out of range [0,%ld) for object %ld\n",
                    icutout, ncutout, iobj);
            return NULL;
        } else {
            return obj;
        }
    }
}

int meds_get_cutout_cen(const struct meds *self,
                        long iobj,
                        long icutout,
                        double *row,
                        double *col)
{
    const struct meds_obj *obj=check_iobj_icutout(self, iobj, icutout);
    if (!obj) {
        (*row)=MEDS_DEFVAL;
        (*col)=MEDS_DEFVAL;
        return 0;
    }
    (*row) = obj->cutout_row[icutout];
    (*col) = obj->cutout_col[icutout];
    return 1;
}

long meds_get_source_file_id(const struct meds *self,
                             long iobj,
                             long icutout)
{
    const struct meds_obj *obj=check_iobj_icutout(self, iobj, icutout);
    if (!obj) {
        return MEDS_DEFVAL;
    }
    return obj->file_id[icutout];
}

const struct meds_image_info *meds_get_source_info(const struct meds *self,
                                                   long iobj,
                                                   long icutout)
{
    long file_id = meds_get_source_file_id(self, iobj, icutout);
    if (file_id < 0) {
        return NULL;
    }
    return &self->image_info->data[file_id];
}

const char *meds_get_source_path(const struct meds *self,
                                 long iobj,
                                 long icutout)
{
    const struct meds_image_info *info=meds_get_source_info(self,iobj,icutout);
    if (!info) {
        return NULL;
    }
    return info->image_path;
}

const struct meds_distort 
*meds_get_distortion(const struct meds *self,
                     long iobj,
                     long icutout)
{
    const struct meds_obj *obj=check_iobj_icutout(self, iobj, icutout);
    if (!obj) {
        return NULL;
    }

    return obj->distortion;
}




// internal routine to read from the fits file
//
// input start_row is zero offset and gets converted to 1-offset further down

static double *get_cutout_data_dbl(const struct meds *self,
                                   const char *extname,
                                   long start_row, long npix)
{
    int status=0;
    double *pix=alloc_doubles(npix);

    if (fits_movnam_hdu(self->fits, IMAGE_HDU, (char*)extname, 0, &status)) {
        fits_report_error(stderr,status);
        return NULL;
    }

    if (!fits_load_sub_dbl(self->fits, start_row, npix, pix)) {
        free(pix);
        pix=NULL;
    }
    return pix;
}
static int *get_cutout_data_int(const struct meds *self,
                                   const char *extname,
                                   long start_row, long npix)
{
    int status=0;
    int *pix=alloc_ints(npix);

    if (fits_movnam_hdu(self->fits, IMAGE_HDU, (char*)extname, 0, &status)) {
        fits_report_error(stderr,status);
        return NULL;
    }

    if (!fits_load_sub_int(self->fits, start_row, npix, pix)) {
        free(pix);
        pix=NULL;
    }
    return pix;
}


// internal routine to get info for a cutout and read from
// the specified extension
static double *get_cutoutp_dbl(const struct meds *self,
                               const char *extname,
                               long iobj,
                               long icutout,
                               long *nrow,
                               long *ncol)
{
    double *pix=NULL;

    const struct meds_obj *obj=check_iobj_icutout(self, iobj, icutout);
    if (!obj) {
        goto _meds_get_cutoutp_dbl_bail;
    }

    long start_row=obj->start_row[icutout];

    *nrow=obj->box_size;
    *ncol=obj->box_size;
    long npix = (*nrow)*(*ncol);

    pix=get_cutout_data_dbl(self, extname, start_row, npix);

_meds_get_cutoutp_dbl_bail:
    if (!pix) {
        *nrow=0;
        *ncol=0;
    }
    return pix;
}
static int *get_cutoutp_int(const struct meds *self,
                            const char *extname,
                            long iobj,
                            long icutout,
                            long *nrow,
                            long *ncol)
{
    int *pix=NULL;

    const struct meds_obj *obj=check_iobj_icutout(self, iobj, icutout);
    if (!obj) {
        goto _meds_get_cutoutp_int_bail;
    }

    long start_row=obj->start_row[icutout];

    *nrow=obj->box_size;
    *ncol=obj->box_size;
    long npix = (*nrow)*(*ncol);

    pix=get_cutout_data_int(self, extname, start_row, npix);

_meds_get_cutoutp_int_bail:
    if (!pix) {
        *nrow=0;
        *ncol=0;
    }
    return pix;
}


// internal routine to get info for a mosaic and read from
// the specified extension
static double *get_mosaicp_dbl(const struct meds *self,
                               const char *extname,
                               long iobj,
                               long *ncutout,
                               long *nrow,
                               long *ncol)
{
    double *pix=NULL;

    const struct meds_obj *obj=check_iobj_icutout(self, iobj, 0);
    if (!obj) {
        goto _meds_get_mosaicp_dbl_bail;
    }

    long start_row=obj->start_row[0];

    *ncutout = obj->ncutout;
    *nrow    = obj->box_size;
    *ncol    = obj->box_size;

    long npix = (*nrow)*(*ncol)*(*ncutout);

    pix=get_cutout_data_dbl(self, extname, start_row, npix);

_meds_get_mosaicp_dbl_bail:
    if (!pix) {
        *ncutout=0;
        *nrow=0;
        *ncol=0;
    }
    return pix;
}
static int *get_mosaicp_int(const struct meds *self,
                            const char *extname,
                            long iobj,
                            long *ncutout,
                            long *nrow,
                            long *ncol)
{
    int *pix=NULL;

    const struct meds_obj *obj=check_iobj_icutout(self, iobj, 0);
    if (!obj) {
        goto _meds_get_mosaicp_int_bail;
    }

    long start_row=obj->start_row[0];

    *ncutout = obj->ncutout;
    *nrow    = obj->box_size;
    *ncol    = obj->box_size;

    long npix = (*nrow)*(*ncol)*(*ncutout);

    pix=get_cutout_data_int(self, extname, start_row, npix);

_meds_get_mosaicp_int_bail:
    if (!pix) {
        *ncutout=0;
        *nrow=0;
        *ncol=0;
    }
    return pix;
}




//
// public image reading routines
//

// cutouts as pointers
double *meds_get_cutoutp(const struct meds *self,
                         long iobj,
                         long icutout,
                         long *nrow,
                         long *ncol)
{
    return get_cutoutp_dbl(self,"image_cutouts",iobj,icutout,nrow,ncol);
}

double *meds_get_mosaicp(const struct meds *self,
                         long iobj,
                         long *ncutout,
                         long *nrow,
                         long *ncol)
{
    return get_mosaicp_dbl(self,"image_cutouts",iobj,ncutout,nrow,ncol);
}

double *meds_get_weight_cutoutp(const struct meds *self,
                                long iobj,
                                long icutout,
                                long *nrow,
                                long *ncol)
{
    return get_cutoutp_dbl(self,"weight_cutouts",iobj,icutout,nrow,ncol);
}

double *meds_get_weight_mosaicp(const struct meds *self,
                                long iobj,
                                long *ncutout,
                                long *nrow,
                                long *ncol)
{
    return get_mosaicp_dbl(self,"weight_cutouts",iobj,ncutout,nrow,ncol);
}

int *meds_get_seg_cutoutp(const struct meds *self,
                          long iobj,
                          long icutout,
                          long *nrow,
                          long *ncol)
{
    return get_cutoutp_int(self,"seg_cutouts",iobj,icutout,nrow,ncol);
}

int *meds_get_seg_mosaicp(const struct meds *self,
                          long iobj,
                          long *ncutout,
                          long *nrow,
                          long *ncol)
{
    return get_mosaicp_int(self,"seg_cutouts",iobj,ncutout,nrow,ncol);
}



// cutouts as meds_cutout structures
struct meds_cutout *meds_get_cutout(const struct meds *self,
                                    long iobj,
                                    long icutout)
{

    long nrow=0, ncol=0;
    double *pix=meds_get_cutoutp(self, iobj, icutout, &nrow, &ncol);
    if (!pix) {
        return NULL;
    }

    struct meds_cutout *cutout=cutout_from_ptr(pix, 1, nrow, ncol);
    return cutout;
}

struct meds_cutout *meds_get_mosaic(const struct meds *self, long iobj)
{

    long ncutout=0, nrow=0, ncol=0;
    double *pix=meds_get_mosaicp(self, iobj, &ncutout, &nrow, &ncol);
    if (!pix) {
        return NULL;
    }

    struct meds_cutout *cutout=cutout_from_ptr(pix, ncutout, nrow, ncol);
    return cutout;
}

struct meds_cutout *meds_get_weight_cutout(const struct meds *self,
                                           long iobj,
                                           long icutout)
{

    long nrow=0, ncol=0;
    double *pix=meds_get_weight_cutoutp(self, iobj, icutout, &nrow, &ncol);
    if (!pix) {
        return NULL;
    }

    struct meds_cutout *cutout=cutout_from_ptr(pix, 1, nrow, ncol);
    return cutout;
}

struct meds_cutout *meds_get_weight_mosaic(const struct meds *self, long iobj)
{

    long ncutout=0, nrow=0, ncol=0;
    double *pix=meds_get_weight_mosaicp(self, iobj, &ncutout, &nrow, &ncol);
    if (!pix) {
        return NULL;
    }

    struct meds_cutout *cutout=cutout_from_ptr(pix, ncutout, nrow, ncol);
    return cutout;
}

struct meds_icutout *meds_get_seg_cutout(const struct meds *self,
                                        long iobj,
                                        long icutout)
{

    long nrow=0, ncol=0;
    int *pix=meds_get_seg_cutoutp(self, iobj, icutout, &nrow, &ncol);
    if (!pix) {
        return NULL;
    }

    struct meds_icutout *cutout=icutout_from_ptr(pix, 1, nrow, ncol);
    return cutout;
}

struct meds_icutout *meds_get_seg_mosaic(const struct meds *self, long iobj)
{

    long ncutout=0, nrow=0, ncol=0;
    int *pix=meds_get_seg_mosaicp(self, iobj, &ncutout, &nrow, &ncol);
    if (!pix) {
        return NULL;
    }

    struct meds_icutout *cutout=icutout_from_ptr(pix, ncutout, nrow, ncol);
    return cutout;
}

