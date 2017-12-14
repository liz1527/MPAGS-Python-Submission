#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:04:52 2017

Code created for the MPAGS Coursework by Lizzie Elmer.

This code takes an imput of 3 fits catalogues from UDS observations. It will 
restrict the region this data comes from to a subset of the UDS field and 
extract the magnitude values for 7 year stacks of data. It then applies a crude
correction for variations in seeing across the years before calculating a 
measure of variability (median absolute deviation) for each object and plotting 
this against the average magnitude of the object.
The plot produced is also interactive - if you click on one of the blue crosses
it will bring up the light curve of that object.

I have attached pngs of the main plot and an example lightcurve.

@author: ppxee
"""

### Import required libraries ###
import matplotlib.pyplot as plt #for plotting
from astropy.io import fits #for handling fits
import numpy as np #for handling arrays
from astropy.stats import median_absolute_deviation
plt.close('all') #close any open plots

### Define user functions ###
def mag5_stacks(tbdata):
    ''' Function that takes a catalogue of magnitude data from the sextracor 
    output and makes a np array containing only the 3 arcsec aperture data for 
    each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        flux = an array with 8 columns containing flux values for each year '''
        
    flux = np.stack(([tbdata['MAG_APER_5_05B'],
                tbdata['MAG_APER_5_07B'], tbdata['MAG_APER_5_08B'],
                tbdata['MAG_APER_5_09B'], tbdata['MAG_APER_5_10B'], 
                tbdata['MAG_APER_5_11B'], tbdata['MAG_APER_5_12B']]), axis=1)
    return flux

def magerr5_stacks(tbdata):
    ''' Function that takes a catalogue of magnitude data from the sextracor 
    output and makes a np array containing only the 3 arcsec aperture error 
    data for each epoch
    Input:
        tbdata = the original combined catalogue of flux data 
    Output:
        fluxerr = an array with 8 columns containing flux error values for each
        year '''
        
    fluxerr = np.stack([tbdata['MAGERR_APER_5_05B'],
                tbdata['MAGERR_APER_5_07B'],tbdata['MAGERR_APER_5_08B'],
                tbdata['MAGERR_APER_5_09B'],tbdata['MAGERR_APER_5_10B'], 
                tbdata['MAGERR_APER_5_11B'],tbdata['MAGERR_APER_5_12B']], axis=1)
    return fluxerr

def psf_correct(baseflux, initflux, avgtype):
    ''' Function that applies a fractional PSF correction to initflux, based on
    the constant required to correct an epochs average flux to the overall 
    average flux in the baseflux array.
    i.e. if doing a correction based on all objects for objects in the chandra
    data then baseflux = flux and initflux = fluxchan.
    If doing a correction based on the stellar fluxes for all objects in the 
    field then baseflux = sflux and initflux = flux. 
    Basetype must be defined as either 'all' or 'star' so that the mean value
    can be used for all objects and the median value can be used for stellar 
    objects 
    Inputs:
        baseflux = flux array to base the corrections on
        initflux = flux array to be corrected
        basetype = either 'mean' or 'median' which dictate if the mean or 
                    median value is used for the correction 
    Output:
        Flux array with values crudely corrected for differences in seeing
        (average flux should now be the same for each epoch). '''
    if avgtype == 'mean':
        avgfluxperepoch = np.mean(baseflux, axis=0)#for UDS
        avgflux = np.mean(baseflux)
        const = avgflux/avgfluxperepoch
    elif avgtype == 'median':
        avgfluxperepoch = np.median(baseflux, axis=0)#for UDS
        avgflux = np.median(baseflux)
        const = avgflux/avgfluxperepoch
    else:
        print('Invalid basetype')
        return
    return initflux * const[None,:]

def err_correct(flux, fluxerr, fluxnew):
    ''' Function that applies a correction to the array of error values that 
    matches the correction applied to the corresponding array of fluxes.
    Inputs:
        flux = initial flux array before any corrections were applied
        fluxerr = initial flux err array
        fluxcorr = array of fluxes that have been corrected
    Output:
        Flux error array with values crudely corrected '''
        
    return fluxnew * (fluxerr/flux)


def lightcurve5(ob, fitsdata)  :
    ''' Function that plots the light curve of an object in terms of its flux 
    in an aperture 5 pixels across (i.e. 3 arcsec in diameter) 
    Inputs:
        ob = the ID of the object that you want the lightcurve from
        fitsdata = the original catalogue of data that the curve will be 
                    plotted from 
    Output:
        None '''
        
    #Create flux array for chandra area
    tbdata = chandra_only(fitsdata)
    fluxwhole = mag5_stacks(tbdata)
    
    #Get data for the object called
    mask = fitsdata['NUMBER_05B'] == ob
    obdata = fitsdata[mask]
    if not obdata: #Reject if no object number matches the input value
        print('error- invalid object number')
        return
    #Create arrays of flux values and error values
    flux = mag5_stacks(obdata)
    fluxerr = magerr5_stacks(obdata)

    # normalise and correct for seeing
    fluxcorr = psf_correct(fluxwhole, flux, 'mean')
    fluxerrcorr = err_correct(flux, fluxerr, fluxcorr)
    fluxcorr = np.squeeze(fluxcorr)
    fluxerrcorr = np.squeeze(fluxerrcorr)
    
    #set up time variable for plot
    t = np.array([1,3,4,5,6,7,8])
    years = ('05B', '06B', '07B', '08B', '09B', '10B', '11B', '12B')
    #Plot graph in new figure
    plt.figure()
    plt.errorbar(t, fluxcorr, yerr=fluxerrcorr, fmt = 'ro')
    plt.xticks(np.linspace(1,8,8), years)
    plt.xlabel('Semester')
    plt.ylabel('Magnitude of object')
    plt.title('Light curve for object number %i' % ob)
    return

def chandra_only(tbdata):
    ''' Function that restricts the objects included in analysis to only 
    those within the Chandra footprint 
    Input:
        tbdata = original catalogue of data 
    Output:
        newtbdata = new catalogue of data which only includes objects within 
                    the chandra footprint '''
    ### Restrict objects to those in the Chandra field ###
    mask1 = tbdata['DELTA_J2000_05B'] < -4.93 #max Dec
    mask2 = tbdata['DELTA_J2000_05B'] > -5.403 #min Dec
    mask3 = tbdata['ALPHA_J2000_05B'] < 34.72 #max RA
    mask4 = tbdata['ALPHA_J2000_05B'] > 34.07 #min RA
    mask = mask1 * mask2 * mask3 * mask4
    newtbdata = tbdata[mask]
    return(newtbdata)

def normalise(flux):
    ''' Normalise each objects flux to its average value
    Input:
        flux = array of object flux values 
    Output:
        array of object flux values normalised to the average flux of the 
        object '''
    avgflux = np.mean(flux, axis=1)
    return flux / avgflux[:,None]

def onpickchanonly(event):
    ''' Function that plots the lightcurve of an object when it is clicked on 
    the vairiability v flux plot '''
    
    combined = fits.open('mag_flux_table_best.fits')
    tbdata = combined[1].data
    tbdata = chandra_only(tbdata)
    
    ### remove values that are +/-99 ###
    fluxn = mag5_stacks(tbdata)
    fluxn[fluxn == 99] = np.nan
    mask = ~np.isnan(fluxn).any(axis=1)
    tbdata = tbdata[mask]
    
    ob = tbdata['NUMBER_05B'][event.ind] #Define the object number from the index of the selected point
    if len(ob) > 1: #Reject selection if more than one object has been selected
        print('Too many objects selected')
        return
    print('Object identified')
    lightcurve5(ob, tbdata) #Plot the lightcurve from the lightcurve function

def flux_variability_plot(flux, fluxchan, starflux=[],
                          normalised = False, stars=False):
    ''' Function to plot the variability vs mean flux plot using the MAD 
    statistic. Optional input of fluxes of stars or normalising the plot.
    Inputs:
        flux = array of flux values for UDS objects
        fluxchan = array of flux values for chandra objects
        starflux = optional array of fluxes for stars
        normalised = True or False (default) depending if the fluxes should be
                        normalised to the objects average flux 
        stars = True or False (default) depending on if stars should be added
                to the plot
    Output:
        fig = figure handle to allow clicking for light curves to be enabled if
                required '''
    
    fig = plt.figure()
    avgfluxperob = np.mean(flux, axis=1) #for UDS
    avgfluxchanperob = np.mean(fluxchan, axis=1) #for non-stellar chandra
    if stars==True:
        savgfluxperob = np.mean(starflux, axis=1) #for stars

    ### Check if normalisation is true and normalise if necessary ###
    if normalised == True:
        flux = normalise(flux)
        fluxchan = normalise(fluxchan) 
        if stars == True:
            starflux = normalise(starflux)
    ### Find out which plot type is specified and calculate appropriate statistic ###
    vary = median_absolute_deviation(flux, axis=1)
    varychan = median_absolute_deviation(fluxchan, axis=1)
          
    ### Plot the variability v mean as appropriate ###
    if stars==True:
        varystar = median_absolute_deviation(starflux, axis=1)
        plt.plot(savgfluxperob, varystar, 'm*', mfc = 'none', markersize = 10,
                 label='Secure Star') 
    line, = plt.plot(avgfluxperob, vary, 'b+', label='UDS Source', picker=2)
    plt.plot(avgfluxchanperob, varychan, 'ro', mfc = 'none', markersize = 10,
             label='Chandra Source') #no picker as will be selected in the UDS point

        
    ### Apply required plot charateristics ###
    plt.yscale('log')
    plt.xlabel('Mean Magnitude')
    plt.ylabel('MAD')
    plt.legend()
    
    return fig

def no99(flux):
    ''' Function to remove objects that have a magnitude value of 99 in any of
    the year stacks. This is because that means a negative flux was found at 
    that point by SExtractor 
    Input:
        flux = array of magnitudes with all objects
    Output:
        flux[mask] = array of magnitudes with only those objects that have 
                        actual values in all year stacks. '''
    flux[flux == 99] = np.nan
    mask = ~np.isnan(flux).any(axis=1)
    return flux[mask]


### Open the fits files and get data ###
tbdata = fits.open('mag_flux_table_best.fits')[1].data
chandata = fits.open('chandra_mag_flux_table_best.fits')[1].data
sdata = fits.open('stars_mag_flux_table.fits')[1].data

### Restrict objects to those in the Chandra field ###
# Field is restricted so can compare the number of variable sources that are
# and aren't visible in x-rays
tbdata = chandra_only(tbdata)
sdata = chandra_only(sdata)

## Create arrays of flux values from each year stack ###
fluxn = mag5_stacks(tbdata) # for all objects
fluxchann = mag5_stacks(chandata) # for chandra objects
sfluxn = mag5_stacks(sdata) # for secure stars

### Remove values that = 99 as these are points where the flux is negative ###
fluxn = no99(fluxn)
fluxchann = no99(fluxchann)
sfluxn = no99(sfluxn)

### Apply a crude correction for psf variations from year to year ###
fluxcorrn = psf_correct(fluxn, fluxn, 'mean') 
fluxchancorrn = psf_correct(fluxn, fluxchann, 'mean') 
sfluxcorrn = psf_correct(fluxn, sfluxn, 'mean') 

### Create plot of MAD vx mean magnitude ###
fig = flux_variability_plot(fluxcorrn, fluxchancorrn, starflux=sfluxcorrn,
                            normalised=True, stars=True)

### Allow for lightcurves to pop up when an object is selected ###
fig.canvas.mpl_connect('pick_event', onpickchanonly)
