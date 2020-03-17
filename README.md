# Andicam Optical & Near-Infrared Photometry Pipeline For The Point-like Sources

This respository contains a pipeline for point source photometry on CCD images collected with the ANDICAM instrument mounted on the CTIO 1.3m telescope ([DePoy et al. 2003](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/4841/1/A-Novel-Double-Imaging-Camera-ANDICAM/10.1117/12.459907.short)).

SMARTSOpt.py - The optical photometry module

SMARTSIR.py - The near-infrared photometry module

SOptfull.py - The wrapper for the optical photometry

SIRfull.py - The wrapper for the near-infrared photometry


### Prerequisites

Before using the pipeline, you need to install the following python packages: os, glob, string, linecache, subprocess, numpy, pandas, pyraf, astropy, alipy.
To install them with pip: 

```
pip install os, glob, string, linecache, subprocess, numpy, pandas, pyraf, astropy, alipy
```

You will also need the "daophot" stellar photometry package of Peter B. Stetson if you wanna do psf-photometry. This package is only available by its developer upon request after signing an agreement.


### Steps for the optical photometry

1. Gather the images into a directory.
2. Align the images with respect to a reference frame.
3. Do the photometry
4. Collect the photometry results

### Steps for the infrared photometry
1. Gather the images into a directory
2. Process the dithered images (make a sky noise image, denoise each of the dithered images, and then align and combine them to amplify the signal)
3. Align the combined images with respect to ta reference frame
4. Do the photometry
5. Collect the photometry results

The usage of the pipeline for optical and infrared photometry can be found in the examples folder.
