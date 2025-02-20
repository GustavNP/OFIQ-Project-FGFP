# Open Source Face Image Quality (OFIQ) with Fine-Grained Face Parsing (FGFP)

The OFIQ project has been modified to compute Fine-Grained Face Parsing (FGFP) masks. It works by dividing the Face Skin region of the existing face parsing algorithm into more smaller sub-regions.

The Fine-Grained Face Parsing algorithm is implemented in the file [OFIQLib/src/OFIQImpl.cpp](OFIQLib/src/OFIQImpl.cpp).

The project has only been tested in Windows.

For the Fine-Grained Face Parsing algorithm to work, three folders have to be manually added to the installation folder "\OFIQ-Project-FGFP\install_x86_64\Release\bin\". The three folders that have to be added manually:
  - aligned_images
  - FGFP_images
  - score_files
    
In "aligned_images", cropped and aligned images are saved after being computed in OFIQ. In "FGFP_images", Fine-Grained Face Parsing masks are saved. In "score_files" UQS and CQM values are saved in csv files.

The project is installed in the same way as the original OFIQ project. See [BUILD.md](BUILD.md) for details and dependencies, otherwise look here below.

Build commands (for Windows):
  -  (if not already installed) pip install conan==2.0.17
  - cd \<path-to-OFIQ-Project-FGFP\>\\\\scripts
  - .\\build.cmd

Running OFIQ-FGFP:
  - cd \<path-to-OFIQ-Project-FGFP\>\\\\install_x86_64\\\\Release\\\\bin
  - Single image:
    - .\\\\OFIQSampleApp -c ..\\\\..\\\\..\\\\data\\\\ -i ..\\\\..\\\\..\\\\data\\\\tests\\\\images\\\\b-01-smile.png
  - All images in directory:
    - .\\\\OFIQSampleApp -c ..\\\\..\\\\..\\\\data\\\\ -i ..\\\\..\\\\..\\\\data\\\\tests\\\\images -o quality-scores-and-measures.csv

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

