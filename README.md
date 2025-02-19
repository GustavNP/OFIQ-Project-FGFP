# Open Source Face Image Quality (OFIQ) with Fine-Grained Face Parsing (FGFP)

The OFIQ project has been modified to compute Fine-Grained Face Parsing (FGFP) masks. It works by dividing the Face Skin region of the existing face parsing algorithm into more smaller sub-regions.

The Fine-Grained Face Parsing algorithm is implemented in the file [OFIQLib/src/OFIQImpl.cpp](OFIQLib/src/OFIQImpl.cpp).



README from original repository:

# Open Source Face Image Quality (OFIQ)

The __OFIQ__ (Open Source Face Image Quality) is a software library for computing quality 
aspects of a facial image. OFIQ is written in the C/C++ programming language.
OFIQ is the reference implementation for the ISO/IEC 29794-5 international
standard; see [https://bsi.bund.de/dok/OFIQ-e](https://bsi.bund.de/dok/OFIQ-e).

## License
Before using __OFIQ__ or distributing parts of __OFIQ__ one should have a look
on OFIQ's license and the license of its dependencies: [LICENSE.md](LICENSE.md)
  
## Getting started
For a tutorial on how to compile and operate OFIQ, see [here](BUILD.md).

## Reference manual
A full documentation of __OFIQ__ including compilation, configuration and a comprehensive doxygen documentation of 
the C/C++ API is contained in the reference manual:
see [doc/refman.pdf](doc/refman.pdf).

