KEPLER TIME SERIES OF DETRENDED, NORMALIZED PHOTOMETRY

INTRODUCTION

    This Readme file describes the ASCII files containing time series of
detrended, normalized stellar fluxes of individual target stars as observed
by NASA's Kepler Mission.  These files have names ending in "dtr.txt".
These notes apply to version 1.0 HLSP data (identified as such in the 
filenames). These notes were compiled by Tim Brown, 18 Dec 2009.

    For a general description of the Kepler Mission and its aims, see
Borucki et al., Science (2010, in press) and references therein.

    Each dtr.txt file corresponds to a single target, which is identified
in the filename by its index number in the Kepler Input Catalog (KIC).
Each file begins with metadata and target-specific data encoded in FITS-like
keywords, as described below.  The latter part of the file is a 2-column
table containing the dates and corresponding detrended normalized flux
measurements for the target.

DATA PROCESSING

    Raw CCD images obtained by the spacecraft are extensively processed to
yield the tabulated detrended relative fluxes.  Time series data that are
archived elsewhere in the MAST will be reduced using slightly
different algorithms than for the early release described herein.  The
functional descriptions of the two processes are quite similar, however.
The processing chain is described in Caldwell et al., ApJ Lett. (2010, in
press), and in Jenkins et, al., Ap J. Lett (2010a,b in press).  The time
series described here are from the Kepler quarterly roll positions Q0 and Q1;
the governing release for these data is #2, and the corresponding release
notes on the MAST (KSCI-19042) largely apply.

     Processing of the "Long Cadence" data includes the following steps:
*  Data are averaged on-board over intervals of ~1766 s by co-adding single-
   integration CCD frames.
*  At the end of the ~30 minute long cadence period, pre-specified pixels for
   each target are selected from the coadd.
*  A mean bias level is subtracted.
*  The pixels are requantized so that the effective noise due to quantization
   is a constant percentage of the intrinsic noise for all signal levels.
*  The pixels are encoded using a lossless Huffman technique, and stored in 
   packets in preparation for downlink.

     After downlink:
*  We remove cosmic ray signals on a pixel-by-pixel basis.
*  Compute raw fluxes for each target as a weighted sum of the pixel 
   intensities for that target ("Simple Aperture Photometry").
*  Detrend and normalize these raw fluxes by forming the ratio with a
   1-day running median centered on each time sample.  For time series
   with identified transit signals, samples that fall within transits are
   masked out of the median average.

    The timestamps associated with each time sample are very precise, but
at present the correction to HJD accounts only for the boresight coordinates
of the spacecraft.  Because of Kepler's large field of view, the HJD for
a particular target may be different from this value by as much as about 30 s.
In a future version of the data analysis, we will correct this inaccuracy and
also report Barycentric Julian Date (BJD), rather than the current HJD.

    The uncertainties in the measured fluxes depend upon crowding of the
stellar images, and are difficult to model accurately.  The relevant parameters
are, however, largely independent of time.  The best estimate of photometric
uncertainty is thus derived by examination of the light-curve itself.  We
provide such an estimate in the metadata quantity FLUXRMS.  This is
the RMS of the detrended normalized time series, after removing
points that are discrepant by more than 5 sigma.  Here sigma is 
estimated as IQR/1.349, where IQR is the separation between the
1st and 3rd quartile points in a sorted list of normalized fluxes covering
the entire data set.  The constant 1.349 in this formula is the expected IQR
for a Gaussian distribution having unit variance.

DESCRIPTION OF TABLE METADATA

    FILENAME:  The name of the MAST data file containing the current table.
    TARGNAME:  The index number of the target in the Kepler Input Catalog (KIC).
        KIC data may be retrieved from the following url:
        http://archive.stsci.edu/kepler/kepler_fov/help/search_help.html
    TELESCOP:  The Kepler spacecraft contains one telescope: "Kepler".
    INSTRUME:  The Kepler spacecraft contains one instrument: "Photometer".
    RA_TARG:   [deg] The target's right ascension, epoch J2000, obtained
        from the KIC.
    DEC_TARG:  [deg] The target's declination, epoch J2000, obtained from
        the KIC.
    RMAGTARG:  The target's Sloan r-magnitude, obtained from the KIC.
    G-R_TARG:  The target's Sloan g-r color, obtained from the KIC.
    DATEMIN:   The HJD of the earliest data sample.
    DATEMAX:   The HJD of the latest data sample.
    CAD_TIME:  [s] The time cadence interval, estimated as the median time 
        between samples.
    TIME_FIL:  The fraction of time in the interval that is covered by 
        observations.

    FLUXRMS:   The RMS of the detrended normalized time series, described above.
    NPTS:      The number of time samples in the series.

DESCRIPTION OF TABLE COLUMNS

    HJD:  Heliocentric Julian Date of the center time of the integration
        for the corresponding time sample.  As noted above, in version 1
        this value for each target is subject to a nearly-constant offset
        of up to +/- 30 s.   
    Detr_Norm_Flux: Detrended noramlized flux from the target, averaged over
        the duration of the time sample.  Because of the normalization, this
        quantity must have a median value that is very close to unity.
