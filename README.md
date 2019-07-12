# Andicam Optical/NIR photometry pipeline for point-like sources

A repository to aid the point source photometry on images collected with the ANDICAM onboard the CTIO 1.3m telescope. Both aperture and psf techniques are possible. 

SMARTSOpt.py - the main optical photometry module

SMARTSIR.py - the main near-infrared photometry module

SOptfull.py - the wrapper for the optical photometry

SIRfull.py - the wrapper for the near-infrared photometry

These modules should not be used without a good knowledge of all steps in the code, especially the near-infrared one due to the complexity in the background construction.
