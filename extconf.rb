require "mkmf"

dir_config('narray',$sitearchdir,$sitearchdir)
dir_config('fftw3','/usr/local')

if ( ! ( have_header("narray.h") && have_header("narray_config.h") ) ) then
   print <<-EOS
   ** configure error **  
   Header narray.h or narray_config.h is not found. If you have these files in 
   /narraydir/include, try the following:

   % ruby extconf.rb --with-narray-include=/narraydir/include

EOS
   exit(-1)
end

if ( ! ( have_header("fftw3.h") && have_library("fftw3") ) ) then
   print <<EOS
   ** configure error **
   Header fftw3.h or the compiled fftw3 library is not found.
   If you have the library installed under /fftw3dir (that is, fftw3.h is
   in /fftw3dir/include and the library in /fftw3dir/lib/),
   try the following:

   % ruby extconf.rb --with-fftw3-dir=/fftw3dir

   Alternatively, you can specify the two directory separately
   with --with-fftw3-include and --with-fftw3-lib.
EOS
   exit(-1)
end

if have_library("fftw3f")
  $CFLAGS += ' -DFFTW3_HAS_SINGLE_SUPPORT'
end

if /cygwin|mingw/ =~ RUBY_PLATFORM
   have_library("narray") || raise("ERROR: narray library is not found")
end

create_makefile("numru/fftw3")
