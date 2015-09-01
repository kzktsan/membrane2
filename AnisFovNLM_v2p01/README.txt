-----------------------------------------------------------------------

  Image Filtering Software based on Nonlocal Foveated Self-Similarity
                Release ver. 2.01  (October 4, 2014)

-----------------------------------------------------------------------

Copyright (c) 2011-2014 Tampere University of Technology. 
All rights reserved.
This work should be used for nonprofit purposes only.

Authors:                    Alessandro Foi
                            Giacomo Boracchi

web page:                   http://www.cs.tut.fi/~foi/FoveatedNL/


-----------------------------------------------------------------------
 Contents
-----------------------------------------------------------------------

The package implements the method(s) published in [3] and contains
the following files:

*) demo_AnisFovNLM.m        : main demo script
*) AnisFovNLM.m             : Anisotropic Foveated Nonlocal Means filter
*) makeFovPSFs.m            : builds the constrained set of blurring PSFs that compose the foveation operator
*) makeAllPSFsTiled         : draws the mosaic of blurring PSFs that compose the fovation operator
*) Man_fragment256.tif      : small image for running the demo
*) README.txt               : (this file)
*) LEGAL_NOTICE.txt         : legal information related to this files
*) FovNLM_v1p12.zip         : older codes from [1], ver. 1.12 (isotropic only)


-----------------------------------------------------------------------
 Requirements
-----------------------------------------------------------------------

*) Matlab ver. 6.5 or later


-----------------------------------------------------------------------
 Instructions
-----------------------------------------------------------------------

Instructions, syntax, and usage examples can be found in the comment
headers of the Matlab files.


-----------------------------------------------------------------------
 Change log
-----------------------------------------------------------------------

v2.00  (9 October 2014)
 + Separate loop for windowing

v2.00  (4 October 2014)
 + Anisotropic codes released

v1.12  (16 February 2012)
 + Improved comments and notation

v1.11  (13 February 2012)
 + Improved comments

v1.1   (8 February 2012)
 + Code optimization and first public release
 
v1.0  (28 December 2011) 
 + Minor fixes (padding etc.)

v0.9  (30 July 2011)
 + Initial code

 
-----------------------------------------------------------------------
 References
-----------------------------------------------------------------------
 
[1] A. Foi and G. Boracchi, "Foveated self-similarity in nonlocal image filtering",
    Proc. IS&T/SPIE EI2012 - Human Vision and Electronic Imaging XVII, 8291-32,
    Burlingame (CA), USA, Jan. 2012.
[2] A. Foi and G. Boracchi, "Anisotropically Foveated Nonlocal Image Denoising",
    Proc. IEEE Int. Conf. Image Process. (ICIP 2013), pp. 464-468, Melbourne, Australia, Sep. 2013.
[3] A. Foi and G. Boracchi, "Foveated nonlocal self-similarity", preprint, Oct. 2014.

 All available online at http://www.cs.tut.fi/~foi/
 
 
-----------------------------------------------------------------------
 Disclaimer
-----------------------------------------------------------------------

Any unauthorized use of these routines for industrial or profit-
oriented activities is expressively prohibited. By downloading 
and/or using any of these files, you implicitly agree to all the 
terms of the TUT limited license, as specified in the document
LEGAL_NOTICE.txt (included in this package) and online at
http://www.cs.tut.fi/~foi/FoveatedNL/legal_notice.html


-----------------------------------------------------------------------
 Feedback
-----------------------------------------------------------------------

If you have any comment, suggestion, or question, please do
contact   Alessandro Foi  at  firstname.lastname@tut.fi


