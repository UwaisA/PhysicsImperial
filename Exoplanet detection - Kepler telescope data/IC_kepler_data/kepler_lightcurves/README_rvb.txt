KEPLER-RELATED TIME SERIES OF RADIAL VELOCITIES AND BISECTOR SPANS

INTRODUCTION

    This Readme file describes the ASCII files containing time series of
radial velocity and spectral line bisector span for ground-based
spectroscopic observations of individual target stars that have been observed
by NASA's Kepler Mission.  These files have names ending in "rvb.txt".
These notes apply to version 1.0 HLSP data (identified as such in the
filenames). These notes were compiled by Tim Brown, 18 Dec 2009.

    For a general description of the Kepler Mission and its aims, see
Borucki et al., Science (2010, in press) and references therein.

    Each rvb.txt file corresponds to a single target, which is identified
in the filename by its index number in the Kepler Input Catalog (KIC).
Each file begins with metadata and target-specific data encoded in FITS-like
keywords, as described below.  The latter part of the file is a 6-column
table containing the dates and corresponding radial velocities, bisector
spans, uncertainties, and telescope/spectrograph combinations used for each 
observation of the target.

DATA PROCESSING

    The data processing to yield radial velocity and line bisector estimates
is complex, and varies to account for characteristics of the different
telescope/spectrograph combinations.  

     For the Keck/HIRES system most used in the current context, a description
of the reduction may be found in Johnson et al., PASP 121, 1104 (2009).
For this case, reduction of the CCD images from the spectrograph begins by
correcting the images for bias and fixed-pattern noise, removing radiation
events, correcting for flat-field effects, correcting for instrumental
stray light and background sky light, and extracting 1-dimensional spectra
for each echelle order.  Each order is then broken into small segments,
each of which contains both lines from the target star and calibration lines
from an iodine absorption cell or a Th-Ar emission lamp.  The segments are
anlyzed individually to yield estimates of the spectrograph point-spread
function, as well as Doppler shifts between the stellar and calibration lines.
Radial velocities are derived from differences between the stellar and
calibration Doppler shifts.  Uncertainties are estimated from scatter in
the values obtained from the many independent spectrum segments. 
    
    For NOT-FIES spectra, reduction of the CCD images begins by
correcting the images for bias, median combining three or more
consecutive exposures to effectively remove radiation events, correcting
for flat-field and pixel-to-pixel effects, correcting for instrumental
stray light, and optimally extracting one-dimensional spectra for each
echelle order. A ThAr frame is taken both before and after the stellar
exposure in order to create a separate wavelength solution for each
stellar spectrum from the co-added ThAr frames.
A multi-order cross correlation is performed to determine the radial
velocities, using most of the orders, leaving out orders with a poor
wavelength solution due to few ThAr lines and orders infected by
telluric absorption. The spectra are moved to a common log-lambda
wavelength scale by linear interpolation and all the spectra are
combined to a co-added template spectrum by shifting the spectra to a
common zero point and co-adding them. Each spectrum is then correlated
against this co-added template, and the resulting cross correlation
functions from the different orders are co-added, thus automatically
weighing the orders by their flux level. The center of the cross
correlation function is found by fitting the peak of the cross
correlation function. Uncertainties are estimated from the
order-to-order scatter of the center of the cross correlation function.

    The line bisector span (averaged over lines) is computed using the
techniques described by Queloz et al., A&A 379, 279 (2001) and by
Torres et al., ApJ 619, 558 (2005).  Bisector span uncertainties are formal,
and are derived from the fitting process that yields the bisector spans.

    The timestamps associated with each time sample are the Barycentric
Julian Date (BJD) of the flux-weighted centroid of the observation.

DESCRIPTION OF TABLE METADATA

    FILENAME:  The name of the MAST data file containing the current table.
    TARGNAME:  The index number of the target in the Kepler Input Catalog (KIC).
        KIC data may be retrieved from the following url:
        http://archive.stsci.edu/kepler/kepler_fov/help/search_help.html
    RA_TARG:   [deg] The target's right ascension, epoch J2000, obtained
        from the KIC.
    DEC_TARG:  [deg] The target's declination, epoch J2000, obtained from
        the KIC.
    RMAGTARG:  The target's Sloan r-magnitude, obtained from the KIC.
    G-R_TARG:  The target's Sloan g-r color, obtained from the KIC.
    DATEMIN:   The BJD of the earliest data sample.
    DATEMAX:   The BJD of the latest data sample.
    NPTS:      The number of time samples in the series.

DESCRIPTION OF TABLE COLUMNS

    BJD:  Barycentric Julian Date of the flux-weighted centroid of the 
        integration for the corresponding time sample.
    Radial_Velocity:  [m/s] Radial velocity of the target.  The zero point
        is arbitrary and chosen to yield a small time-averaged radial velocity.
    RV_Uncertainty:  [m/s] Formal uncertainty of the radial velocity.
    Bisector_Span:  [m/s] Difference between the line bisector measured at
        40% and 90% line depth, averaged over all spectrum lines. 
    BS_Uncertainty:  [m/s] Formal uncertainty in the bisector span.
    Telesco/Spectro:  Names of the telescope and spectrograph that were used
        to obtain the observation.
