# Andicam Optical & NIR photometry pipeline for point-like sources

This respository contains a pipeline for point source photometry on images collected with the ANDICAM instrument onbard the CTIO 1.3m telescope. Both aperture and psf techniques are possible. 

SMARTSOpt.py - the main optical photometry module

SMARTSIR.py - the main near-infrared photometry module

SOptfull.py - the wrapper for the optical photometry

SIRfull.py - the wrapper for the near-infrared photometry

## Getting Started

### Prerequisites

os, glob, alipy, string, linecache, subprocess, numpy, pandas, pyraf, astropy, and daophot (optional for psf-photometry).
