=module NumRu::FFTW3

Fast Fourier Transforms by using ((<FFTW|URL:http://www.fftw.org>)) Ver.3.

Takeshi Horinouchi

(C) Takeshi Horinouchi / GFD Dennou Club,
2003

NO WARRANTY

==Features

* Uses ((<NArray|URL:http://www.ruby-lang.org/en/raa-list.rhtml?name=NArray>)).
* Multi-dimensional complex FFT. (Real data are coerced to complex).
* Supports both double and single float transforms.
* Not normalized as in FFTW

==Features yet to be introduced

* Sine / cosine transforms
* User choice of optimization levels (i.e., FFTW_MEASURE etc in
  addition to FFTW_ESTIMATE).
* Multi-threaded FFT3 support -- don't know whether it's really feasible.

==Installation

* Install ((<FFTW|URL:http://www.fftw.org>)) Ver.3.

  * NOTE: 
    To activate the single-float transform, you have to install FFTW3 with
    the single-float compilation, in addition to the default double-float
    version. This can be done by configuring FFTW3 with
    the --enable-float option, and install it again. The single-float
    version will coexist with the double-float version.
    If you do not install the single-float version, FFT is always done
    with the double precision, which is not bad if you are not time- and 
    memory-conscious.

* Install ((<NArray|URL:http://www.ruby-lang.org/en/raa-list.rhtml?name=NArray>)).

* Then, install this library as follows (replace "version" with
  the actual version number):

    % tar xvzf fftw3-version.tar.gz
    % cd fftw3-version
    % ruby extconf.rb
    % make
    % make site-install
  Or
    % make install
  (If you are using Ruby 1.8, make install is the same make site-install.)

==How to use

See the following peice of code. (Install this library and copy and
paste the following to the interactive shell irb).

  require "narray"
  require "numru/fftw3"
  include NumRu

  na = NArray.float(8,6)   # float -> will be corced to complex
  na[1,1]=1

  # <example 1>
  fc = FFTW3.fft(na, -1)/na.length  # forward 2D FFT and normalization
  nc = FFTW3.fft(fc, 1)       # backward 2D FFT (complex) --> 
  nb = nc.real                # should be equal to na except round errors  

  # <example 2>
  fc = FFTW3.fft(na, -1, 0) / na.shape[0]  # forward FFT with the first dim

  # <example 3>
  fc = FFTW3.fft(na, -1, 1) / na.shape[1]  # forward FFT with the second dim

==API Reference

===Module methods

---fft(narray, dir [,dim,dim,...])

    Complex FFT.

    The 3rd, 4th,... arguments are optional.

    ARGUMENTS
    * narray (NArray or NArray-compatible Array) : array to be
      transformed. If real, coerced to complex before transformation.
      If narray is single-precision and the single-precision
      version of FFTW3 is installed (before installing this module),
      this method does a single-precision transform. 
      Otherwise, a double-precision transform is used.
    * dir (-1 or 1) : forward transform if -1; backward transform if 1.
    * optional 3rd, 4th,... arguments (Integer) : Specifies dimensions 
      to apply FFT. For example, if 0, the first dimension is
      transformed (1D FFT); If -1, the last dimension is used (1D FFT);
      If 0,2,4, the first, third, and fifth dimensions
      are transformed (3D FFT); If entirely omitted, ALL DIMENSIONS
      ARE SUBJECT TO FFT, so 3D FFT is done with a 3D array.
 
    RETURN VALUE
    * a complex NArray

    NOTE
    * As in FFTW, return value is NOT normalized. Thus, a consecutive
      forward and backward transform would multiply the size of
      data used for transform. You can normalize, for example,
      the forward transform FFTW.fft(narray, -1, 0, 1)
      (FFT regarding the first (dim 0) & second (dim 1) dimensions) by
      dividing with (narray.shape[0]*narray.shape[1]). Likewise,
      the result of FFTW.fft(narray, -1) (FFT for all dimensions)
      can be normalized by narray.length.



