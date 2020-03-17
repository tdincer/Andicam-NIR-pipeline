# Andicam Optical & Near-Infrared photometry pipeline for point-like sources

This respository contains a pipeline for point source photometry on CCD images collected with the ANDICAM instrument onbard the CTIO 1.3m telescope.

SMARTSOpt.py - the main optical photometry module

SMARTSIR.py - the main near-infrared photometry module

SOptfull.py - the wrapper for the optical photometry

SIRfull.py - the wrapper for the near-infrared photometry

## Prerequisites

Before using the pipeline, you need to install the following python packages: os, glob, string, linecache, subprocess, numpy, pandas, pyraf, astropy, alipy.

```
pip install os, glob, string, linecache, subprocess, numpy, pandas, pyraf, astropy, alipy
```

You will also need the "daophot" stellar photometry package of Peter B. Stetson if you wanna do psf-photometry. This package is only available by its developer upon request after signing an agreement.


## Getting Started
#### Steps for the optical photometry

1. Gather the images into a directory.
2. Align the images with respect to a reference frame.
3. Do the photometry
4. Collect the photometry results

#### Steps for the infrared photometry
1. Gather the images into a directory
2. Process the dithered images (make a sky noise image, denoise each of the dithered images, and then align and combine them to amplify the signal)
3. Align the combined images with respect to ta reference frame
4. Do the photometry
5. Collect the photometry results
