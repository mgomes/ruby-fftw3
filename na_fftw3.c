/*
  na_fftw3.c

  FFT using FFTW Ver.3 (www.fftw.org)

    (C) Takeshi Horinouchi
    NO WARRANTY.

*/

#include <ruby.h>
#include "narray.h"
#include <fftw3.h>

VALUE rb_mFFTW3;
VALUE mNumRu;

static VALUE
#ifdef FFTW3_HAS_SINGLE_SUPPORT
na_fftw3_double(int argc, VALUE *argv, VALUE self)
  /* to be called by na_fftw3 */
#else
na_fftw3(int argc, VALUE *argv, VALUE self)
  /* to be called directly */
#endif
{
  VALUE val, vdir;
  struct NARRAY *a1, *a2;
  int i, dir, *shape, *bucket;
  fftw_plan p;
  fftw_complex *in, *out;
  volatile VALUE v1, v2;

  if (argc<2){
    rb_raise(rb_eArgError, "Usage: fftw(narray, direction [,dim0,dim1,...])");
  }
  val = argv[0];
  vdir = argv[1];

  dir = NUM2INT(vdir);
  if ( dir != 1 && dir != -1 ){
    rb_raise(rb_eArgError, "direction should be 1 or -1");
  }
  v1 = na_cast_object(val, NA_DCOMPLEX);
  GetNArray(v1,a1);
  v2 = na_make_object( NA_DCOMPLEX, a1->rank, a1->shape, CLASS_OF(v1) );
  GetNArray(v2,a2);

  shape = ALLOCA_N(int, a2->rank);
  for (i=0; i<a2->rank; i++){
      shape[i] = a2->shape[a2->rank-1-i];
  }
  in = (fftw_complex*)a1->ptr;
  out = (fftw_complex*)a2->ptr;

  if (argc==2) {
      /* apply FFT to all dimensions */
      p = fftw_plan_dft( a2->rank, shape, 
			 in, out, dir, FFTW_ESTIMATE );
  } else {
      /* apply FFT to selected dimensions (by using the Guru interface) */
      { /* introduce a new scope for additonal local variables */
	  int fft_rank, howmany_rank, ib, j, jf, je, dim;
	  fftw_iodim *fft_dims, *howmany_dims;
	  int *dimids;
	  fft_rank = argc - 2;
	  fft_dims = ALLOCA_N(fftw_iodim, fft_rank);
	  dimids = ALLOCA_N(int, fft_rank);
	  howmany_rank = fft_rank + 1;
	  howmany_dims = ALLOCA_N(fftw_iodim, howmany_rank);
	  
	  for (i=2;i<argc;i++){
	      dim = NUM2INT(argv[i]);
	      if (dim<0) dim += a2->rank;  /* negative: count from the end */
	      if (dim<0 || dim>=a2->rank){
		  rb_raise(rb_eArgError, "dimension < 0 or >= rank");
	      }
	      dimids[i-2] = a2->rank - 1 - dim;
	      if ( i>2 && dimids[i-2] == dimids[i-3] ){
		  rb_raise(rb_eArgError, "redundant -- a same dimension is reppeated");
	      }
	  }
	  
	  /* bukcet sort in increasing order */
	  bucket = ALLOCA_N(int,a2->rank);
	  for(j=0; j<a2->rank; j++) bucket[j] = 0; /* initialize */
	  for(i=0; i<fft_rank; i++) bucket[ dimids[i] ] = 1;
	  for(j=0,i=0; j<a2->rank; j++) {
	      if (bucket[j]==1){
		  dimids[i] = j;
		  i++;
	      }
	  }

	  for(j=0; j<fft_rank; j++){
	      fft_dims[j].n = shape[ dimids[j] ];
	      fft_dims[j].is = 1;
	      for (i=dimids[j]+1 ; i<a2->rank ; i++){
		  fft_dims[j].is *= shape[i];
	      }
	      fft_dims[j].os = fft_dims[j].is;
	      /* printf("fft_ %d  n:%d  is:%d\n",j,
		       fft_dims[j].n,fft_dims[j].is);*/
	  }
	  for(j=0; j<=fft_rank; j++){
	      howmany_dims[j].n = 1;
	      jf = (j==0) ? 0 : (dimids[j-1]+1) ;
	      je = (j==fft_rank) ? a2->rank : (dimids[j]) ;
	      for (i=jf; i<je; i++){
		  howmany_dims[j].n *= shape[i];
	      }
	      howmany_dims[j].is = 1;
	      if (j<fft_rank){
		  for (i=dimids[j]; i<a2->rank; i++){
		      howmany_dims[j].is *= shape[i];
		  }
	      }
	      howmany_dims[j].os = howmany_dims[j].is;
	      /* printf("how_ %d  n:%d  is:%d\n",j,
		        howmany_dims[j].n,howmany_dims[j].is); */
	  }

	  p = fftw_plan_guru_dft( fft_rank, fft_dims, 
				  howmany_rank, howmany_dims,
				  in, out, dir, FFTW_ESTIMATE );

      }
  }

  fftw_execute(p);
  fftw_destroy_plan(p);

  return v2;
}

#ifdef FFTW3_HAS_SINGLE_SUPPORT

/* sourse code generation of na_fftw3_float:
  Copy na_fftw3_double, and replace 
     fftw --> fftwf
     DCOMPLEX --> SCOMPLEX
 */
static VALUE
na_fftw3_float(int argc, VALUE *argv, VALUE self)
{
  VALUE val, vdir;
  struct NARRAY *a1, *a2;
  int i, dir, *shape, *bucket;
  fftwf_plan p;
  fftwf_complex *in, *out;
  volatile VALUE v1, v2;

  if (argc<2){
    rb_raise(rb_eArgError, "Usage: fftw(narray, direction [,dim0,dim1,...])");
  }
  val = argv[0];
  vdir = argv[1];

  dir = NUM2INT(vdir);
  if ( dir != 1 && dir != -1 ){
    rb_raise(rb_eArgError, "direction should be 1 or -1");
  }
  v1 = na_cast_object(val, NA_SCOMPLEX);
  GetNArray(v1,a1);
  v2 = na_make_object( NA_SCOMPLEX, a1->rank, a1->shape, CLASS_OF(v1) );
  GetNArray(v2,a2);

  shape = ALLOCA_N(int, a2->rank);
  for (i=0; i<a2->rank; i++){
      shape[i] = a2->shape[a2->rank-1-i];
  }
  in = (fftwf_complex*)a1->ptr;
  out = (fftwf_complex*)a2->ptr;

  if (argc==2) {
      /* apply FFT to all dimensions */
      p = fftwf_plan_dft( a2->rank, shape, 
			 in, out, dir, FFTW_ESTIMATE );
  } else {
      /* apply FFT to selected dimensions (by using the Guru interface) */
      { /* introduce a new scope for additonal local variables */
	  int fft_rank, howmany_rank, ib, j, jf, je, dim;
	  fftw_iodim *fft_dims, *howmany_dims;
	  int *dimids;
	  fft_rank = argc - 2;
	  fft_dims = ALLOCA_N(fftw_iodim, fft_rank);
	  dimids = ALLOCA_N(int, fft_rank);
	  howmany_rank = fft_rank + 1;
	  howmany_dims = ALLOCA_N(fftw_iodim, howmany_rank);
	  
	  for (i=2;i<argc;i++){
	      dim = NUM2INT(argv[i]);
	      if (dim<0) dim += a2->rank;  /* negative: count from the end */
	      if (dim<0 || dim>=a2->rank){
		  rb_raise(rb_eArgError, "dimension < 0 or >= rank");
	      }
	      dimids[i-2] = a2->rank - 1 - dim;
	      if ( i>2 && dimids[i-2] == dimids[i-3] ){
		  rb_raise(rb_eArgError, "redundant -- a same dimension is reppeated");
	      }
	  }
	  
	  /* bukcet sort in increasing order */
	  bucket = ALLOCA_N(int,a2->rank);
	  for(j=0; j<a2->rank; j++) bucket[j] = 0; /* initialize */
	  for(i=0; i<fft_rank; i++) bucket[ dimids[i] ] = 1;
	  for(j=0,i=0; j<a2->rank; j++) {
	      if (bucket[j]==1){
		  dimids[i] = j;
		  i++;
	      }
	  }

	  for(j=0; j<fft_rank; j++){
	      fft_dims[j].n = shape[ dimids[j] ];
	      fft_dims[j].is = 1;
	      for (i=dimids[j]+1 ; i<a2->rank ; i++){
		  fft_dims[j].is *= shape[i];
	      }
	      fft_dims[j].os = fft_dims[j].is;
	      /* printf("fft_ %d  n:%d  is:%d\n",j,
		       fft_dims[j].n,fft_dims[j].is);*/
	  }
	  for(j=0; j<=fft_rank; j++){
	      howmany_dims[j].n = 1;
	      jf = (j==0) ? 0 : (dimids[j-1]+1) ;
	      je = (j==fft_rank) ? a2->rank : (dimids[j]) ;
	      for (i=jf; i<je; i++){
		  howmany_dims[j].n *= shape[i];
	      }
	      howmany_dims[j].is = 1;
	      if (j<fft_rank){
		  for (i=dimids[j]; i<a2->rank; i++){
		      howmany_dims[j].is *= shape[i];
		  }
	      }
	      howmany_dims[j].os = howmany_dims[j].is;
	      /* printf("how_ %d  n:%d  is:%d\n",j,
		        howmany_dims[j].n,howmany_dims[j].is); */
	  }

	  p = fftwf_plan_guru_dft( fft_rank, fft_dims, 
				  howmany_rank, howmany_dims,
				  in, out, dir, FFTW_ESTIMATE );

      }
  }

  fftwf_execute(p);
  fftwf_destroy_plan(p);

  return v2;
}

static VALUE
na_fftw3(int argc, VALUE *argv, VALUE self)
{
  VALUE val;
  volatile VALUE v1;
  struct NARRAY *a1;

  if (argc<2){
    rb_raise(rb_eArgError, "Usage: fftw(narray, direction [,dim0,dim1,...])");
  }
  val = argv[0];
  v1 = na_to_narray(val);
  GetNArray(v1,a1);
  if(a1->type <= NA_SFLOAT || a1->type == NA_SCOMPLEX ){
      return( na_fftw3_float(argc, argv, self) );
  } else {
      return( na_fftw3_double(argc, argv, self) );
  }

}

#endif

void
 Init_fftw3()
{
  mNumRu = rb_define_module("NumRu");
  rb_mFFTW3 = rb_define_module_under(mNumRu, "FFTW3");
  rb_define_module_function(rb_mFFTW3, "fft", na_fftw3, -1);
}
