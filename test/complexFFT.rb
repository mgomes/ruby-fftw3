require 'rubygems'
require 'narray'
require "numru/fftw3"
include NumRu

print "\n**TEST** all dimensions\n\n"

na = NArray.float(8,4).fill(1)  # will be corced to complex
na[1,1]=5
p na
fc = FFTW3.fft(na, -1)/na.length 
p fc
p fc.real

p FFTW3.fft(fc, 1).real

print "\n**TEST** single float (treated as single if lib fftw3f exits)\n"
print " --- see http://www.fftw.org/fftw3_doc/Precision.html for more info\n\n"
na = NArray.sfloat(8,4).indgen!
fc = FFTW3.fft(na, -1)/na.length 
p fc
p FFTW3.fft(fc, 1).real

print "\n**TEST** dimension selection\n\n"

fc = FFTW3.fft(na, -1, 0)/na.shape[0]
p fc
p FFTW3.fft(fc, 1, 0).real
fc = FFTW3.fft(na, -1, 1)/na.shape[1]
p fc
p FFTW3.fft(fc, 1, 1).real

na = NArray.float(4,3,8,3)
na[1,1,1,0]= 1
p( fc=FFTW3.fft(na, -1, 0,2) / (na.shape[0]*na.shape[2]) )
p( fc=FFTW3.fft(na, -1, 1) / na.shape[1] )
p( fc=FFTW3.fft(na, -1, 0,1,2)  / (na.shape[0]*na.shape[1]*na.shape[2]) )
p FFTW3.fft(fc, 1, 0,1,2).real


