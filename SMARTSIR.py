#!/usr/bin/env python

"""
This is a python module developed to reduce the SMARTS near-infrared data taken in the dithered mode.
"""

_author_ = 'Tolga Dincer'
_credits_ = ['Imran Hasan', 'Emily MacPherson']
_version_ = '1.4.2'
_laststableversion_ = '1.4.1'
_email_ = ['tolga.dincer@yale.edu', 'tolgadincer@gmail.com']
_maintainer_ = 'Tolga Dincer'
_status_ = 'Under Development'  # Prototype, Development, Production

import os
import sys
import glob
import alipy
import string
import linecache
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyraf import iraf
from astropy.time import Time
from astropy.io import fits, ascii


# load iraf packages
iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)
iraf.imred
iraf.noao.obsutil()
#iraf.nproto()

def filefindcp(inputdir = ['/net/xrb/ir'], source=['GX339-4'], band = ['J', 'H', 'K'],
               daterange = ['2017-01-01', '2017-02-28'], workingdir='./'):
    """
    INPUTS:
    -------
    inputdir: list
        The location of the images. Images can be located at several different directories.
    output: str
        Output file name. The list of the images that meet with the criteria will be written to this file.
    source: list
        The sources to be searched.
    band: list
        The filters to be searched.
    daterange: list
        The date range of the images to be searched. Date range must be given as the SMARTS date,
        NOT the observation date.
    workingdir: str
        The folder where the files will be copied into.
    """

    allflist = []
    for einputdir in inputdir:
        allflist = allflist+glob.glob(einputdir+'/binir*[0-9]*.*[0-9]*.fits.gz')

    temp_allflist = []
    for efile in allflist:
        dateefile = float(efile.split('/')[-1][5:11])
        if dateefile >= float(''.join(daterange[0][2:].split('-')))\
                and dateefile <= float(''.join(daterange[1][2:].split('-'))):
            temp_allflist.append(efile)

    allflist = temp_allflist

    temp_inflist = []
    for efile in allflist:
        # Read the header to get the following information: the object name, JD, filename, and filter.
        fitsfile = fits.open(efile)
        source_f = fitsfile[0].header['object']
        jd_f = fitsfile[0].header['jd']
        smartsdate_f = fitsfile[0].header['filename'][2:8]
        band_f = fitsfile[0].header['irfltid']
        location = workingdir+efile.split('/')[-1]

        inflist = [source_f, jd_f, band_f, location, smartsdate_f]
        print inflist
        # Select the fits images that match the selection criteria
        if source_f in source and band_f in band:
            print efile + ' matches with the selection criteria!'
            os.system('cp -v '+efile+' ' + workingdir)
            if efile[-3:] == '.gz':
                os.system('gunzip '+workingdir+'/'+efile.split('/')[-1])
            temp_inflist.append(inflist)

    headinf = pd.DataFrame(temp_inflist, columns=['Source', 'JulianDate', 'Filter', 'File', 'SmartsDate'])
    headinf.to_csv('obsdithlist.csv')


def readobsmenu(obsmenufile='obsmenu.list'):
    """
    Reads obsmenu.list file into pandas DataFrame and returns it. The file contains
    smartsdate, file, band, date, time and object information in tab separated format.
    If no input is given, it will try to read obsmenu.list file in the working directory.
    Input:
    ------
    obsmenufile: str
        This is the name of the file to be read. This can also be a full path to the file
        if run from a different location.
    Returns:
    --------
    obsmenu: pandas DataFrame.
    """
    if os.path.isfile(obsmenufile):  # if obsmenufile exists then read it, else print a message.
        obsmenu = pd.read_csv(obsmenufile, sep='\t', names=['smartsdate', 'file', 'band', 'date', 'time', 'object'])
        return obsmenu
    else:
        print obsmenufile+' does not exist!'


def readfilteredobsmenu(filteredobsmenufile):
    """
    Reads filteredobsmenu file into pandas DataFrame and returns the DataFrame. The file contains
    smartsdate, file, band. date, time and object information in csv separated format.
    Input:
    ------
    filteredobsmenu: str
        This is the name of the file to be read. This can also be a full path to the file.
    Returns:
    --------
    filteredobsmenu: pandas DataFrame.
    """
    if os.path.isfile(filteredobsmenufile):  # if filteredobsmenufile exists then read it, else print a message.
        filteredobsmenu = pd.read_csv(filteredobsmenufile, header=0)
        return filteredobsmenu
    else:
        print filteredobsmenufile+' does not exist!'


def listobs(expr):
    """
    Saves a list of all observations in the working directory into obsmenu.list file.
    INPUT:
    ------
    expr: str
        Expression to filter the data. e.g.: expr = 'OWNER=="YALE-03A-0001"' or expr = 'OWNER=="STANDARD"'
        The keyword OWNER is from the header of the fits file.
    """
    # Clean intermediate step products if there are any leftovers from previous run.
    tobedeleted = ['smartsdates', 'objectmenutemp.list', 'objectmenu.list']
    for i in tobedeleted:
        if i in glob.glob(i):
            os.system('rm '+i)

    # Find all the sources that agrees with the selection criteria and write them into obsmenutemp.list file
    iraf.hselect.unlearn()
    iraf.hselect(images='binir*', fields='$I,IRFLTID,DATE-OBS,TIME-OBS,OBJECT', expr=expr,
                 Stdout='obsmenutemp.list')
    # Here create the smarts dates from file names.
    os.system('awk "{print $1}" obsmenutemp.list | cut -c 6-11 > smartsdates')
    # Paste the smarts dates into the first column
    os.system('paste smartsdates obsmenutemp.list > obsmenu.list')
    # Remove intermedate step file
    os.system('rm smartsdates obsmenutemp.list')

    print "created obsmenu.list file."


def filterobslist(obslist, source, daterange):
    """
    Filters the obsmenu.list, filters it and saves the filtered list with the filteredobsmenu.csv file name
    INPUT:
    ------
    obslist: str
        This is obsmenu.list file
    source: list of str
        A list of targets. e.g. ['GX339-4', 'LMCX-3']
    datarange: list of str
        Date range. e.g. ['02-10-2015','03-10-2015']
    """
    if os.path.isfile('filteredobsmenu.csv'):  # if filteredobsmenu.csv exist, delete it.
        os.system('rm filteredobsmenu.csv')
        #print 'removed previous filteredobsmenu.csv'

    obsmenu = readobsmenu(obslist)  # read thhe obsmenu.list file

    # check the emptiness of the daterange. if it's not empty then apply the date query.
    # if it is empty then go on with the obsmenu.list information in hand.
    if daterange != []:
        dateexpr = "(" + string.replace(daterange[0][2:], '-', '') + \
                   ". <= smartsdate <= " + string.replace(daterange[1][2:], '-', '') + ".)"
        tempobsmenu = obsmenu.query(dateexpr)
    else:
        tempobsmenu = obsmenu

    # filter by target name
    objexpr = []
    if source != []:
        for i in source:
            objexpr.append('(object == "' + i + '")')

        totobjexpr = ' | '.join(objexpr)
        tempobsmenu = tempobsmenu.query(totobjexpr)

    sortedfilteredlist = tempobsmenu.sort_values(['object', 'band', 'date'])  # sort the table

    sortedfilteredlist.to_csv('filteredobsmenu.csv', index=False)  # write the table into file.
    print 'Filteredobsmenu.csv file created.'


def flatList_ir(inputdir, daterange):
    """
    This function finds the images of the infrared flats for a given date range. Then copies
    the images into the appropriate directories.

    The function creates the file 'irflats.csv'.
    irflats.csv would look like:
     ,Source,SMARTSpath,Date,JulianDate,Filter,File
    0,/net/xrb/ir/binir150803.0111.fits,2015-08-04,2457238.60834,K,./data_ir_temp/K/binir150803.0111.fits

    Keyword arguments:
    @param inputdir: List of directories to be searched. e.g. ['/net/xrb/ir','/net/xrb/ir/processed']
    @type inputdir: list of str

    @param daterange: The date range for the flats to be searched. The dates are the SMARTS dates
        meaning that the dates that are in the image names. e.g. ['2015-08-01','2015-09-01']
    @type daterange: list of string
    """
    allflist = []
    for einputdir in inputdir:
        allflist = allflist + glob.glob(einputdir + '/ir*[0-9]*.flat*.fits*')

    temp_allflist = []
    for efile in allflist:
        dateefile = float(efile.split('/')[-1][2:8])
        if ((dateefile >= float(''.join(daterange[0][2:].split('-')))) and
                (dateefile <= float(''.join(daterange[1][2:].split('-'))))):
            temp_allflist.append(efile)

    allflist = temp_allflist

    # create directory structure
    if not os.path.isdir('./flat_ir'):
        os.system('mkdir flat_ir')

    # convert daterange into the Julian date
    t = Time(daterange).jd

    temp_inflist = []

    for efile in allflist:
        # read the header of the flat and get the necessary information.
        date_strip = efile.split('/')[-1][2:8]
        date_f = float(efile.split('/')[-1][2:8])
        jd_f = Time('20' + date_strip[0:2] + '-' + date_strip[2:4] + '-' + date_strip[4:6]).jd
        band_f = efile.split('/')[-1][13].upper()
        location = './flat_ir/' + 'n_' + efile.split('/')[-1][:-3]

        inflist = [efile, date_f, jd_f, band_f, location]

        # apply the time criteria and eliminate the flats outside the given date range
        if jd_f >= t[0] and jd_f <= t[1]:
            # copy the flat image to the appropriate direcotory e.g. "./flat_ir/V"
            os.system('cp -v ' + efile + ' ./flat_ir/')
            temp_inflist.append(inflist)

    # write the flat information into the output file
    headinf = pd.DataFrame(temp_inflist, columns=['SMARTSpath', 'Date', 'JulianDate', 'Filter', 'File'])
    headinf.to_csv('irflats.csv')

    os.chdir('./flat_ir')
    for zfits in glob.glob('*fits.gz'):
        print zfits
        expr = 'gunzip %s' % (zfits)
        p = subprocess.Popen(expr, shell=True)
        p.communicate()

        flatfile = zfits[:-3]
        newflatfilename = "n_" + flatfile
        iraf.imstat.unlearn()
        flatinfo = iraf.imstat(flatfile, Stdout=1)
        meanflat = float(flatinfo[1].split()[2])

        iraf.imarith.unlearn()
        iraf.imarith(flatfile, '/', meanflat, newflatfilename, divzero=1.)

    os.chdir('../')

    return


def dithsummary():
    filelisting = pd.read_csv('obsdithlist.csv')
    grdate = filelisting.groupby(['SmartsDate', 'Filter', 'Source'])

    allinfo = []
    for name, group in grdate:
        date = name[0]
        color = name[1]
        obj = name[2]
        noffile = len(group.File)
        infoforone = [obj, date, color, noffile]

        allinfo.append(infoforone)

    headinf = pd.DataFrame(allinfo, columns=['Source', 'SmartsDate', 'Filter', 'NFile'])
    headinf.to_csv('logofdith.csv')


def nearestflat(jdlist, color):
    """
    This function calculates an average time for a set of dithering images,
    and using this time finds the flat image that was observed closest in time.
    Returns a string that gives the path to this flat.
    """
    # find the average observation time
    meantime = np.nanmean(jdlist)
    # load up the data frame of flats
    flatdf = pd.read_csv('irflats.csv')
    # grab the flats with the matching filters
    sflatdf = flatdf[flatdf.Filter == color]
    ### print 'Possible flats: ', sflatdf
    # subtract the average time form all the flat times, take absolute value, find index of min,
    # use this index to grab the file name with this observation date
    sflatdf = sflatdf.reset_index(drop=True)
    luckyflat = sflatdf.loc[np.argmin(np.abs(sflatdf.JulianDate.values.astype(float) - meantime))].File
    return luckyflat


def imconstrcut(ref_ir, obsstrategy, target, masking):

    listofrawimages = pd.read_csv('obsdithlist.csv')
    listofrawimages.File = [x[2:-3] for x in listofrawimages.File]
    grpdrawimages = listofrawimages.groupby(['SmartsDate', 'Filter'])

    # Loop over dates
    for name, group in grpdrawimages:
        bandpass = name[1]  # find the filter
        filelist = group.File  # make a list of all raw (dithered) images

        if len(filelist)/obsstrategy[bandpass] < 1:  # if there is not enough file in the group for the corresponding strategy, then skip the group
            print 'Not enough images for the corresponding strategy!'
            print 'So skipping the group with following files:'
            print filelist
        elif len(filelist)/obsstrategy[bandpass] >= 1:
            splittedfilelist = np.array_split(filelist, len(filelist)/obsstrategy[bandpass])
            # print 'length of splittedlist: ', len(splittedfilelist)

            ditgroupcounter = 1
            for ditgroup in splittedfilelist:
                binirimages = list(ditgroup)
                print binirimages
                scaleimages(binirimages)
                writecountingstatstoafile(binirimages, ditgroupcounter, bandpass, target[0])

                skyim = makesky(binirimages, 0)  # 0 means do not use object masking to make the sky image.
                sbinirimages = subtractskyrescaling(binirimages,
                                                     skyim)  # Returns the names of the skysubtracted images.
                flatim = nearestflat(group.JulianDate, bandpass)  # This will find the flat nearest in time.
                fsbinirimages = flatfield(sbinirimages, flatim, null=1)

                if masking == 1:
                    mfsbinirimages = objmask(fsbinirimages)  # This will create the mfsbinirimages.
                    cbinirimages = cpimages(binirimages)  # This will make a copy of the binirimages with names "cbinir...".
                    addobjmasktoheader(cbinirimages, mfsbinirimages)
                    skyim2 = makesky(cbinirimages, 1)  # Create a new sky from object-masked images.
                    scbinirimages = subtractskyrescaling(cbinirimages, skyim2)
                    fscbinirimages = flatfield(scbinirimages, flatim, null=1)
                elif masking == 0:
                    fscbinirimages = cpimages_nomasking(fsbinirimages) # This will change the filenames sfbinirimages to scfbinirimages.


                align(fscbinirimages, bandpass, ref_ir)
                print 'Aligning finished!'

                alignedfiles = []
                inputfiles = map(lambda x: './data_ir_temp/' + bandpass + '/fsc' + x[0:16] + '_gregister.fits',
                                 filelist)
                print bandpass
                print inputfiles

                for aligned in inputfiles:
                    # print "Aligned files: ", aligned
                    if os.path.isfile(aligned):
                        alignedfiles.append(aligned)

                if len(alignedfiles) >= 3:
                    print "combining the images"
                    combimage(alignedfiles, bandpass, ditgroupcounter, target[0])
                else:
                    print "not enough files were aligned"

                ditgroupcounter += 1
                os.system('rm cbin* fsbin* sky* mfsbin* fscbin* sbin* scbin*')
                rmstr = 'rm ./data_ir_temp/%s/*_gregister.fits' % bandpass
                os.system(rmstr)


def makesky(images, withmask):
    imagesinirafformat = str(images).replace("', '", ",").replace("['", "").replace("']", "")
    if withmask == 0:
        skyname = 'sky0.' + images[0][6:12] + '.fits'  # so the format will be sky0.yymmdd.fits
        iraf.imcombine.unlearn()
        iraf.imcombine(imagesinirafformat, skyname, combine='median', reject='minmax', nlow=0, nhigh=1)#, scale=scales)
    elif withmask == 1:
        skyname = 'sky1.' + images[0][7:13] + '.fits'  # so the format will be sky1.yymmdd.fits
        iraf.imcombine.unlearn()
        iraf.imcombine(imagesinirafformat, skyname, combine='median', masktype='!OBJMASK')#, reject='minmax', nlow=0, nhigh=1)

    return skyname


def subtractsky(images, sky):
    for im in images:
        iraf.imarith.unlearn()
        iraf.imarith(im, '-', sky, 's' + im)


def subtractskyrescaling(images, sky):
    for im in images:
        iraf.imstat.unlearn()
        imstat = iraf.imstat(im, field='image, mean, midpt, mode, stddev', Stdout=1)
        #meanim = np.mean([float(x.split()[1]) for x in imstat[1:]])
        iraf.imstat.unlearn()
        skystat = iraf.imstat(sky, field='image, mean, midpt, mode, stddev', Stdout=1)
        #meansky = np.mean([float(x.split()[1]) for x in skystat[1:]])
        #addresfact = meanim - meansky
        #addresfact=0.
        iraf.imarith.unlearn()
        #iraf.imarith(sky, '-' , addresfact, 'res'+sky)
        #iraf.imarith(im, '-', 'res'+sky, 's'+im)
        iraf.imarith(im, '-', sky, 's'+im)
        #removetempsky = 'rm %s' % ('res'+sky)
        #os.system(removetempsky)

    skysubtim = ['s'+x for x in images]

    return skysubtim


def objmask(images):
    # This function creates object masks for the input images and writes OBJMASK keyword to the header of the given
    # images. Masking the objects is important when constructing sky/background images.
    for im in images:
        iraf.nproto.objmasks.unlearn()
        iraf.nproto.objmasks(images=im, objmasks='m'+im, hsigma=1.5)

    mimages = ['m'+im for im in images]
    return mimages


def writecountingstatstoafile(images, ditgroupcounter, bandpass, target):
    imagesinirafformat = str(images).replace("', '", ",").replace("['", "").replace("']", "")
    date = images[0][5:11]
    outputfilename = '%s.%s.%s.%s.txt' % (target, date, bandpass, ditgroupcounter)
    iraf.imstat.unlearn()
    iraf.imstat(imagesinirafformat, field='image, mean, midpt, mode, stddev, npix', Stdout=outputfilename)


def imagevar(target, bandpass):
    alldates = np.array([])
    allmeans = np.array([])
    allstds = np.array([])

    filestruct = '%s*%s*txt' % (target, bandpass)
    files = glob.glob(filestruct)
    for ff in files:
        data = ascii.read(ff)
        mean = np.mean(data['MEAN'])
        std = np.std(data['MEAN'])
        date = float(ff.split('.')[1])
        alldates = np.append(alldates, date)
        allmeans = np.append(allmeans, mean)
        allstds = np.append(allstds, std)

    percentage = allstds/allmeans
    plt.scatter(alldates, percentage)
    print alldates
    print allmeans
    print allstds
    print percentage*100.
    figname = '%s.ratio.eps' % (target)
    plt.savefig(figname)


def cpimages_nomasking(images):
    for im in images:
        cpcommand = 'cp %s %s' % (im, 'sc'+im[1:])
        os.system(cpcommand)
    cpim = ['sc'+im[1:] for im in images]
    return cpim


def cpimages(images):
    for im in images:
        cpcommand = 'cp %s %s' % (im, 'c'+im)
        os.system(cpcommand)
    cpim = ['c'+im for im in images]
    return cpim


def combimage(images,bandpass, ditgroupcounter, target):
    imagesinirafformat = str(images).replace("', '", ",").replace("['", "").replace("']", "")
    finalimage = images[0][0:17]+target+'.'+bandpass+'.'+images[0].split('/')[-1][8:14]+'.'+str(ditgroupcounter)+'.fits'
    iraf.imcombine.unlearn()
    iraf.imcombine(imagesinirafformat, finalimage, combine='median', mclip='yes', reject='none', lsigma=3, hsigma=3, nlow=1, nhigh=1, nkeep=1)


def addobjmasktoheader(images, masks):
    for i in np.arange(len(images)):
        mm = masks[i]+'[pl]'
        iraf.hedit(images[i], 'OBJMASK', mm, add='yes', verify='no')


def flatfield(images, flatfile, null=0):
    # This function divides the flat-fields any given images, and returns the name assigned to the flat-fielded image(s).
    # Images have to be in python list format.

    if type(flatfile) == type([]):
        flatfile = flatfile[0]
        
    iraf.imarith.unlearn()
    if null==0:
        for im in images:
            iraf.imarith(im, '/', flatfile, 'f'+im)
    else:
    	for im in images:
    		iraf.imarith(im, '/', 1, 'f'+im)
            
    ffimages = ['f' + x for x in images]

    return ffimages


def scaleimages(images):

    ims = str(images).replace("', '", ", ").replace("']", "").replace("['", "")
    iraf.imstat.unlearn()
    imstat = iraf.imstat(ims, field='image, mean, stddev', Stdout=1)
    means = [float(x.split()[1]) for x in imstat[1:]]
    stds = [float(x.split()[2]) for x in imstat[1:]]
    meanofmeans = np.mean(means)
    upperlimits = np.array(meanofmeans)+2.*np.array(stds)
    lowerlimits = np.array(meanofmeans)-2.*np.array(stds)

    means2 = np.array([])
    stds2 = np.array([])
    for index, im in enumerate(images):
        iraf.imstat.unlearn()
        lower = lowerlimits[index]
        upper = upperlimits[index]
        imstat = iraf.imstat(im, field='image, mean, stddev', Stdout=1, lower=lower, upper=upper)
        means2 = np.append(means2, float(imstat[1].split()[1]))
        stds2 = np.append(stds2, float(imstat[1].split()[2]))
    meanofmeans = np.mean(means2)

    for index, im in enumerate(images):
        scale = meanofmeans/means[index]
        scaledimage = 'sc' + im
        iraf.imarith.unlearn()
        iraf.imarith(im,'*', scale, result=scaledimage)
        oim = 'o' + im
        expr = 'cp %s %s' % (im, oim)
        os.system(expr)
        expr = 'mv %s %s' % (scaledimage, im)
        os.system(expr)

def align(filelist, bandpass, ref_ir):
    """
    This function aligns the dithering images of the same day using the first image as a reference.
    The output folder is defined as outdrc = './data_ir_temp/'+color

    Keyword Arguments:
    @param filelist: A list of the name of dithering images.
    @param color: A string that shows the bandpass used to obtain the images.
    This parameter has no effect on the aligning, it is only used to create the path of the aligned images.
    """
    ref_image = ref_ir[bandpass]
    #ref_image = filelist[0]
    # get a listing of the path for each image for a given filter
    # #everything below here i copied from the alipy demo http://obswww.unige.ch/~tewes/alipy/tutorial.html
    identifications = alipy.ident.run(ref_image, filelist, visu=False)
    # That's it !
    # Put visu = True to get visualizations in form of png files (nice but much slower)
    # On multi-extension data, you will want to specify the hdu (see API doc).

    # The output is a list of Identification objects, which contain the transforms :
    for id in identifications:  # list of the same length as images_to_align.
        if id.ok == True:  # i.e., if it worked
            print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
            # id.trans is a alipy.star.SimpleTransform object. Instead of printing it out as a string,
            # you can directly access its parameters :
            # print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
            # print id.trans.matrixform()
            # print id.trans.inverse() # this returns a new SimpleTransform object
        else:
            print "%20s : no transformation found !" % (id.ukn.name)

    # Minimal example of how to align images :
    outputshape = alipy.align.shape(ref_image)
    # This is simply a tuple (width, height)... you could specify any other shape.
    outdrc = './data_ir_temp/' + bandpass
    for id in identifications:
        if id.ok:
            # Variant 2, using geomap/gregister, correcting also for distortions :
            alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, verbose=False,
                                  shape=outputshape, makepng=False, outdir=outdrc)
    return
