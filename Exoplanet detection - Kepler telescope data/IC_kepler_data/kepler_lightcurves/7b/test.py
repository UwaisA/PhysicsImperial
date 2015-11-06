import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

from numpy import loadtxt
from scipy.optimize import curve_fit

lines = loadtxt("hlsp_exo_kepler_phot_KPLR5780885_kep_v1.0_dtr.txt", comments="#", unpack=False)

# plotting stuff
#rcParams['figure.figsize'] = 10, 6

plt.plot(lines[:,0],lines[:,1],'o-b')
plt.xlabel('Heliocentric Julian Date / Day', fontsize=16)
plt.ylabel('De-trended normalised flux', fontsize=16)
plt.title('Transit data for Kepler 4-b',fontsize=16)
plt.show()
