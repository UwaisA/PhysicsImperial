import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt
from scipy.optimize import curve_fit
import os

pi = np.pi
radJup = 6.991*10**7 #in metres, for conversions
radSun = 6.955*10**8 #in metres, for conversions
path = "IC_kepler_data/kepler_lightcurves/7b/hlsp_exo_kepler_phot_KPLR5780885_kep_v1.0_dtr.txt" #change this to the path of transit data
R=1.843*radSun #change this to radius of parent star (this value is Kepler 7)
edge1 = edge2 = 0 #x-coord of edges of a transit
rList = [] #radii predicted by each transit
vList = [] #velocity predicted by each transit

#retreiving data from file
mydir = os.path.dirname(os.path.realpath(__file__))
datafilepath = os.path.join(mydir, path)
lines = loadtxt(datafilepath, comments="#", unpack=False)
x = lines[:,0]
y = lines[:,1]

x = (x-x[0])*24 #converting x to hours

#draws a line graph using the parameters passed to it
def plot(x, y, colour = 'ro-', xlabel = "Phase(Hours)", ylabel = "Normalised Light Intensity"):
    plt.plot(x,y, colour)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title("Light Curve")

######################################################################
############################## ANALYSIS ##############################
######################################################################

def func(x, a, b, c):
    #a ~ transit depth, b ~ 2/transit time, c = centre of transit on graph
    return 1-a*np.exp(-(b*(x-c))**8)

def singleTransAnalysis(x, y, locOfTrans):
    popt, pcov = curve_fit(func, x, y, p0=[1-np.min(y), 1, locOfTrans])
    #justification for p0 initial values for curve fitting:
    # a - y coord of base of transit will be near lowest point on graph
    # b - no justification, 2 hours is a reasonable transit time 
    # c - x coord of centre of transit to be analysed
    
    global edge1, edge2
    #plotting
    plot(x, y, colour = 'ro-')
    smoothX = np.linspace(x[0], x[-1], num = len(x)*10)
    plt.plot(smoothX,func(smoothX,popt[0],popt[1],popt[2]), 'g-')
    plt.show()
    
    #x coords of edges of BASE of transit on the graph, where 90% or more of exo
    #intersects with star
    edge1 = popt[2] - ((-np.log(0.9))**(1./8))/popt[1]
    edge2 = popt[2] + ((-np.log(0.9))**(1./8))/popt[1]
    baseVals = np.arange(np.argwhere(x > edge1)[0], np.argwhere(x > edge2)[0], 1)
    drop = 1-np.mean(y[baseVals])
    
    #edges changed to where 50% of exo intersects with star as sharpness of
    #transit edge drops can be inaccurate, using 50% can eliminate most of the
    #error caused by this
    edge1 = popt[2] - ((-np.log(0.5))**(1./8))/popt[1]
    edge2 = popt[2] + ((-np.log(0.5))**(1./8))/popt[1]
    print popt[1]
    width = edge2-edge1 #transit width
    #adding current transit's predictions to the list
    rList.append(((drop*R**2)**0.5)/radJup)
    vList.append(2*R/(width*3600))

#averages the results of all transits
def multTransits():
    global edge1, edge2
    transLim = (np.min(y)+1)/2.
    lowPoint = np.min(y) #min value on entire graph
    while(lowPoint < transLim): #repeats iterative process until all transits have been analysed
        singleTransAnalysis(x, y, x[np.argmin(y)]) #analyses transit at 3rd parameter
        cen = (edge1 + edge2)/2.
        w = edge2-edge1
        #setting all of the transit's points to 1 so they are not reanalysed
        y[np.argwhere(x > cen-w)[0]:np.argwhere(x > cen+w)[0]] = 1.
        lowPoint = np.min(y)

multTransits()

print "radius of exoplanet:", np.mean(rList), "Jupiter radii"
print "volume of exoplanet:", 4/3*pi*(np.mean(rList)*11)**3, "Cubic Earth Radii"
print "velocity of exoplanet:", int(round(np.mean(vList))), "metres/second"

#desmos graph link - https://www.desmos.com/calculator/vgjby8ly1a