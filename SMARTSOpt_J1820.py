#!/usr/bin/env python

"""
This is the standard SMARTS module for the analysis of the optical and infrared photometric images.
"""

import os
import glob
import alipy
import string
import linecache
import subprocess
import numpy as np
import pandas as pd
from pyraf import iraf
from astropy.time import Time
from astropy.io import fits, ascii


_author_ = 'Tolga Dincer'
_credits_ = ['Emily MacPherson', 'Imran Hasan']
_version_ = '1.4.2'
_laststableversion_ = '1.4.1'
_email_ = ['tolgadincer@gmail.com']
_maintainer_ = 'Tolga Dincer'
_status_ = 'Production'  # Prototype, Development, Production

# load iraf packages
iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)
iraf.imred(_doprint=0)
iraf.noao.obsutil(_doprint=0)
iraf.prcacheOff()

def synchfolders(dir1, dir2):
    # Find files
    # dir1
    optlist1 = glob.glob(dir1+'/*opt.list')[0]
    obsmenulist1 = glob.glob(dir1+'/obsmenu.list')[0]
    obslist1 = glob.glob(dir1+'/obslist.csv')[0]
    photdoneobsmenu1 = glob.glob(dir1+'/photdoneobsmenu.csv')[0]
    flagalignedfilteredobsmenu1 = glob.glob(dir1+'/flagalignedfilteredobsmenu.csv')[0]
    filteredobsmenu1 = glob.glob(dir1+'/filteredobsmenu.csv')[0]
    aligned1 = glob.glob(dir1+'/aligned')[0]
    raws1 = glob.glob(dir1+'/*raw.csv')

    #dir2
    optlist2 = glob.glob(dir2+'/*opt.list')[0]
    obsmenulist2 = glob.glob(dir2+'/obsmenu.list')[0]
    obslist2 = glob.glob(dir2+'/obslist.csv')[0]
    photdoneobsmenu2 = glob.glob(dir2+'/photdoneobsmenu.csv')[0]
    flagalignedfilteredobsmenu2 = glob.glob(dir2+'/flagalignedfilteredobsmenu.csv')[0]
    filteredobsmenu2 = glob.glob(dir2+'/filteredobsmenu.csv')[0]
    aligned2 = glob.glob(dir2+'/aligned')[0]
    raws2 = glob.glob(dir2+'/*raw.csv')

    # Synchronize them
    #optlist1 = os.system('cp -v '+efile+' '+ workingdir)
    #awk '{print $0}' J1709.Rraw.csv && awk 'NR!=1{print $0}' J1709.Iraw.csv
    




def filefindcp(inputdir = ['/net/xrb/ccd', '/net/xrb/ccd/todo', '/net/xrb/ccd/reduced'], output = 'obsfile.list',
               source=['GX339-4'], band = ['I', 'V', 'B', 'R'], daterange = ['2015-08-16', '2015-08-18'], workingdir='./'):
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
        The date range of the images to be searched. Date range must be given as the SMARTS date, NOT the observation date.
    workingdir: str
        The folder where the files will be copied into.
    """
    allflist = []
    for einputdir in inputdir:
        allflist = allflist+glob.glob(einputdir+'/rccd[0-9]*.fits.gz')
        allflist = allflist+glob.glob(einputdir+'/rccd[0-9]*.fits')

    #allflist = [line.strip() for line in open("allf2", 'r')]

    print daterange


    temp_allflist = []
    for efile in allflist:
        dateefile = float(efile.split('/')[-1][4:10])
        if dateefile >= float(''.join(daterange[0][2:].split('-'))) and dateefile <= float(''.join(daterange[1][2:].split('-'))):
            temp_allflist.append(efile)

    allflist = temp_allflist

    temp_inflist = []
    f = open(output, 'w')
    for efile in allflist:
        # read the header information
        fitsfile = fits.open(efile)
        if efile[-3:] == '.gz':
            efilewogz = efile[:-3]
        else:
            efilewogz = efile
        source_f = fitsfile[0].header['object']
        date_f = fitsfile[0].header['date-obs']
        jd_f = fitsfile[0].header['jd']
        band_f = fitsfile[0].header['ccdfltid']
        location = workingdir+efile.split('/')[-1]

        inflist = [source_f, efile, date_f, jd_f, band_f, location]
        print inflist
    
        # decide on the selection criteria
        if source_f in source and band_f in band:
            print efile
            os.system('cp -v '+efile+' '+ workingdir)
            if efile[-3:] == '.gz':
                os.system('gunzip '+workingdir+'/'+efile.split('/')[-1])  # 'gunzip ./data_temp/'+band_f+'/'+efile.split('/')[-1])
            f.write(efilewogz+'\n') # write the image name into the output.
            temp_inflist.append(inflist)

    f.close()

    headinf = pd.DataFrame(temp_inflist,columns=['Source', 'FName', 'Date', 'JDate', 'Filter', 'Location'])
    headinf.to_csv('obslist.csv')


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
    iraf.hselect(images='rccd*', fields='$I,CCDFLTID,DATE-OBS,TIME-OBS,OBJECT', expr=expr,
                 Stdout='obsmenutemp.list')
    # Here create the smarts dates from file names.
    os.system('awk "{print $1}" obsmenutemp.list | cut -c 5-10 > smartsdates')
    # Paste the smarts dates into the first column
    os.system('paste smartsdates obsmenutemp.list > obsmenu.list')
    # Remove intermedate step file
    os.system('rm smartsdates obsmenutemp.list')

    print "created obsmenu.list file."


def statsobs():
    """
    Note to myself: This part needs to be developed.
    """
    isfile = os.path.isfile('filtered.csv')
    
    if isfile:
        print 'filtered.csv verified.'
    else:
        print 'filtered.csv does not exist!'
        return



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
        print 'removed previous filteredobsmenu.csv'

    obsmenu = readobsmenu(obslist)  # read thhe obsmenu.list file

    # check the emptiness of the daterange. if it's not empty then apply the date query.
    # if it is empty then go on with the obsmenu.list information in hand.
    if daterange != []:
        dateexpr = "(" + string.replace(daterange[0][2:], '-', '') + \
                   ". <= smartsdate <= " + string.replace(daterange[1][2:], '-', '') + ".)"
        print dateexpr
        tempobsmenu = obsmenu.query(dateexpr)
    else:
        tempobsmenu = obsmenu

    # filter by target name
    objexpr = []
    if source != []:
        for i in source:
            print i
            objexpr.append('(object == "'+i+'")')
            print objexpr
        
        totobjexpr = ' | '.join(objexpr)
        print totobjexpr
        tempobsmenu = tempobsmenu.query(totobjexpr)    

    sortedfilteredlist = tempobsmenu.sort_values(['object', 'band', 'date'])  # sort the table

    sortedfilteredlist.to_csv('filteredobsmenu.csv', index=False)  # write the table into file.
    print 'created filteredobsmenu.csv file'


def getrefimage(target, band, pathtorefimages):  # new version
    """
    This function finds the reference image for a given target and bandpass.
    INPUTS:
    -------
    target: str
        Target name. It has to be the same as in the header. e.g. 'GX 339-4'
    band: str
        Name of the bandpass. e.g: 'I'
    pathtorefimages: str
        Path to the reference images. e.g. '/net/xrb/tables'
    Returns:
    --------
        the path to the reference image for a given target and bandpass.
    """
    refimage = glob.glob(pathtorefimages+'/ref_images/'+target+'/'+band+'/*.fits')[0]

    return refimage


def flagaligned(filteredobsmenu, filenames):
    filteredobsmenu = readfilteredobsmenu(filteredobsmenu)
    
    if len(filenames) != 0:
        for i in filteredobsmenu.file:
            if i in filenames:
                filteredobsmenu.loc[filteredobsmenu.file == i, 'aligned'] = int(1)
            else:
                filteredobsmenu.loc[filteredobsmenu.file == i, 'aligned'] = int(0)
    else:
        filteredobsmenu.aligned = int(0)
    flagalignedfilteredobsmenu = filteredobsmenu       
 
    return flagalignedfilteredobsmenu


def flagphotdone(inputforphot, outputofphot):
    """
    INPUTS:
    -------
    inputforphot: pandas data frame 
        the list that is given for the photometry. it is in the format of obsmenu.list
    outputofphot: list
        list of filenames to be compared.
    RETURNS: pandas data frame
    --------
    """
    filenames = outputofphot
    print 'filenames: ', filenames
    if len(filenames) != 0:  # if number of elements in filenames list is not empty
        if outputofphot[0][-3:] == 'als':
            filenames = [x[:-4] for x in outputofphot]
        elif outputofphot[0][-6:] == 'mag.MC':
            filenames =  [x[:-7] for x in outputofphot]
        for i in inputforphot.file:            
            print 'inputforphot: ', i
            if i in filenames:
                print 'Yes, the photometry has been performed on ' + i 
                inputforphot.loc[inputforphot.file == i, 'photdone'] = 1
            else:
                print 'Nope, the photometry has not been performed on ' + i 
                inputforphot.loc[inputforphot.file == i, 'photdone'] = 0
    else:
        print 'There is no photometry done! Whats happening dude?'
        inputforphot.photdone = int(0)
    photdoneobsmenu = inputforphot       
 
    return photdoneobsmenu


def makeonedarr(filelist):
    """
    Transforms multidimensional list to a single list.
    INPUT:
    ------
    filelist: list in list. e.g.: [[1,2,3],[a,b,c,d],[4,5]] ---> [1,2,3,a,b,c,d,4,5]
    Returns:
    --------
    1-d array.
    """
    onedarr = []
    for i in filelist:
        for j in i:
            onedarr.append(j)
    return onedarr


def alignprocessedwrapper(inputobsmenu, pathtorefimages):
    """
    INPUT:
    inputobsmenu: It is a csv file that contains observation info. 
                  smartsdate, file, band, date, time, object  
    """
    # Read inputobsmenu file. It is filtered file and in csv format.
    filteredobsmenu = readfilteredobsmenu(inputobsmenu)

    # Group the observations by source name
    objectfiltergrouped = filteredobsmenu.groupby(['object', 'band'])

    # create aligned directory if it doesn't exist
    if os.path.isdir('./aligned') == False:
        os.system('mkdir aligned')

    tempcounter = 0
    tempfilenames = []
    for name, group in objectfiltergrouped:
        target = name[0]
        band = name[1]
        imagelist = group.file.values
        refimage = getrefimage(target, band, pathtorefimages)
        print refimage
        outdir = './aligned'
        
        counter, filenames = alignprocessed(imagelist, refimage, outdir)
        tempcounter = tempcounter + counter
        tempfilenames.append(filenames)

    filenames = makeonedarr(tempfilenames)
    if len(filenames) == 0:
        os.sys.exit('found 0 aligned files! exiting...')
    else:
        pass

    flagalignedfilteredobsmenu = flagaligned('filteredobsmenu.csv', filenames)

    print "Number of images aligned: {0}".format(tempcounter)
    print 'Images aligned: ', filenames

    if os.path.isfile('flagalignedfilteredobsmenu.csv'):
        os.system('rm flagalignedfilteredobsmenu.csv')
    flagalignedfilteredobsmenu.to_csv('flagalignedfilteredobsmenu.csv', index=False)

    gregtoa()

    return tempcounter, filenames


def findphottype(obj):
    database = {'Swift_J1753.5-0127': {'phot': 'psf'},
    			'J1820': {'phot': 'psf'},
                'J1357': {'phot': 'psf'},
                'J1535-571': {'phot': 'psf'},
                'XTEJ1118': {'phot': 'aper'},
				'AqlX-1': {'phot': 'aper'},       # Psf for optical, aperture for infrared
                'GX339-4': {'phot': 'psf'},
                '4U1543-47': {'phot': 'psf'},
                'GROJ1655-40': {'phot': 'psf'},  # This source is named as something else in earlier data and it causes problem.
                'GS1354-645': {'phot': 'aper'},
                'A0620_YALO' : {'phot': 'aper'},
                'A0620-00' : {'phot': 'aper'},
                'MAXI_J1305-704' : {'phot': 'psf'},
                'CenX-4': {'phot': 'aper'},
                'LMCX-3': {'phot': 'aper'},
                'J1709': {'phot': 'aper'},
                'CXOJ0924': {'phot': 'psf'},
                'LANDOLT-TPheD': {'phot': 'aper'},
                'LANDOLT-RU149': {'phot': 'aper'},
                'LANDOLT-PG1047': {'phot': 'aper'},
				'1124-186': {'phot': 'aper'}}
    photto = database[obj]['phot']

    return photto


def gregtoa():
    files = glob.glob('./aligned/*_gregister.fits')
    if len(files) > 0:
        for i in files:
            newfilenamei = i.replace('_gregister.fits', '.fits')
            newfilename = newfilenamei.replace('rccd', 'arccd')
            expr = 'mv ' + i + ' ' + newfilename
            os.system(expr)


def performphot(flagalignedfilteredobsmenu, pathtorefimages):
    flagalignedfilteredobsmenu = readfilteredobsmenu(flagalignedfilteredobsmenu)
    flagalignedfilteredobsmenualigned = flagalignedfilteredobsmenu.query('aligned==1')  # Choose only the aligned
    # images to work on them.
    listgrouped = flagalignedfilteredobsmenualigned.groupby(['object', 'band'])

    for group, name in listgrouped:
        obj = name.object.iloc[0]
        images = np.array(name.file)
        band = name.band.iloc[0]
        phottype = findphottype(obj)

        print obj
        print band
        print images
        
        if phottype == 'aper':
            print 'performing aperture photometry!'
            performaper(obj, images, band, pathtorefimages)
        elif phottype == 'psf':
            os.chdir('aligned')
            performpsf(obj, images, band, pathtorefimages)
            print 'performing psf photometry!'
            os.chdir('../')

    inputforphot = flagalignedfilteredobsmenu
    outputofphotpath1 = glob.glob('./aligned/*.fits.mag.MC')
    outputofphotpath2 = glob.glob('./aligned/*.fits.als')
    outputofphotpath = outputofphotpath1 + outputofphotpath2
    outputofphot = [x.split('/')[-1].replace('arccd', 'rccd') for x in outputofphotpath] 
    photdoneobsmenu = flagphotdone(inputforphot, outputofphot)

    if os.path.isfile('photdoneobsmenu.csv'):
        os.system('rm photdoneobsmenu.csv')
    photdoneobsmenu.to_csv('photdoneobsmenu.csv', index=False)


# Here I am.
def makerawtables(pathtorefimages):
    photdonemenufile = 'photdoneobsmenu.csv'
    photdonemenu = readfilteredobsmenu(photdonemenufile)
    photdonelist = photdonemenu.query('photdone==1')
    photdonelistgrouped = photdonelist.groupby(['object', 'band'])

    for group, name in photdonelistgrouped:
        obj = name.object.iloc[0]
        band = name.band.iloc[0]
        images = np.array(name.file)
        magfiles = ['./aligned/'+x.replace('rccd', 'arccd').replace('.fits', '.fits.mag.MC') for x in images]
        #alsfiles = ['./aligned/'+x.replace('rccd', 'arccd').replace('.fits', '.fits.als') for x in images]
        photcoordsfile = pathtorefimages+'/ref_images/'+obj+'/'+band+'/photcoords.lis'
        photcoords = ascii.read(photcoordsfile)
        nsources = len(photcoords)  # Number sources
        print 'Obj:', obj
        print 'Band: ', band
        print 'Mag Files: ', magfiles
        magtocsv(obj, band, magfiles, nsources)  # This will read the magnitudes from mag.MC files and
        # put them in a [source].[band]raw.csv file. e.g.: GX 339-4.Iraw.csv


# Here I am.
def magtocsv(obj, band, magfiles, nsources):
    """
    INPUTS:
    -------
    magfiles: str
        a list of *mag.MC files with full path. This list should contain a list for only one source and one filter.
    the keys for this dictionary will be the columns for the data frame.
    each key starts off with an empy python list as its value.
    we will read each output file one by one and append its data to the appropriate lists
    so that each element in the list for a key will represents a row in that column
    """
    # Create the empty dictionary with standards keys and appropriate number of mag and merr
    columns = ['object', 'fname', 'date', 'juliandate', 'airmass']
    rowdict = {'object': [], 'fname': [], 'date': [], 'juliandate': [], 'airmass': []}
    for i in np.arange(nsources):
        magno = 'mag' + str(i+1)
        merrno = 'merr' + str(i+1)
        rowdict[magno] = []
        rowdict[merrno] = []
        columns.append(magno)
        columns.append(merrno)

    # Grab the numbers from magfiles and put them in the rowdict       
    for i in magfiles:
        # use the handy astropy ascii object to read photometry files
        photdata = ascii.read(i)
        print photdata
        # the astropy.ascii automatically masks any bad data values. I'd prefer to keep the bad values as np.nans
        if np.ma.is_masked(photdata['MAG']):
            photdata['MAG'][photdata['MAG'].mask] = np.nan
        if np.ma.is_masked(photdata['MERR']):
            photdata['MERR'][photdata['MERR'].mask] = np.nan
        # for append the data from this photometry file to the dictionary
        rowdict['object'].append(obj)
        rowdict['fname'].append(photdata['IMAGE'][0])
        rowdict['date'].append(float(photdata['IMAGE'][0][5:11]))
        rowdict['juliandate'].append(photdata['OTIME'][0])
        rowdict['airmass'].append(round(photdata['XAIRMASS'][0], 3))

        for j in np.arange(len(photdata['MAG'])):
            magno = 'mag'+str(j+1)
            merrno = 'merr'+str(j+1)
            rowdict[magno].append(photdata['MAG'][j])
            rowdict[merrno].append(photdata['MERR'][j])
    
    # create a pandas dataframe out of the dictionary
    df = pd.DataFrame(rowdict)
    # change the data type of the date column to int. cant translate from string to int above for some reason
    df['date'] = df['date'].astype(int)

    csvname = obj+'.'+band+'raw.csv'
    df.to_csv(csvname, columns=columns, index=False, float_format='%.3f')
    

def alignprocessed(images_to_align, ref_image, outdir):
    """
    align the processed ccd images using alipy, a python module I found that
    calls source extractor and iraf imalign tasks
    INPUT:
    target: a string indicating the target you want to reduced
    fltid: the andicam filter for which these data were taken in.
    ref_image: a string indicating the path to the refence image used to align the other images
    outdir: a string indicating the path to save the new aligned fits images to
    OUTPUT:
    new fits images that are aligned to the reference images. 'geregister' is depended
    to the file name of each input image
    """
    # everything below here i copied from the alipy demo http://obswww.unige.ch/~tewes/alipy/tutorial.html
    # this line uses source extractor to identify stars in the images
    identifications = alipy.ident.run(ref_image, images_to_align, visu=False)
    # That's it
    # Put visu=True to get visualizations in form of png files (nice but much slower)
    # On multi-extension data, you will want to specify the hdu (see API doc).

    # The output is a list of Identification objects, which contain the transforms :
    # These tell you how rotate and translate the images so they are aligned to the referene image
    for id in identifications:  # list of the same length as images_to_align.
        if id.ok:  # i.e., if it worked
            print "%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio)
            # id.trans is a alipy.star.SimpleTransform object. Instead of printing it out as a strin
            # you can directly access its parameters :
            # print id.trans.v # the raw data, [r*cos(theta)  r*sin(theta)  r*shift_x  r*shift_y]
            # print id.trans.matrixform()
            # print id.trans.inverse() # this returns a new SimpleTransform object
        else:
            print "%20s : no transformation found !" % (id.ukn.name)
            # return 0

    # Minimal example of how to align images :
    outputshape = alipy.align.shape(ref_image)
    # This is simply a tuple (width, height)... you could specify any other shape.

    # finally, for each image where a transform was found, create a new image where the data are transformed
    # to be alinged with the reference image, and print it to 'outdir', the last argument of this function

    counter = 0.
    filenames = [] 
    for id in identifications:
        if id.ok == True:
            # Variant 2, using geomap/gregister, correcting also for distortions :
            alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, verbose=False,
                                  shape=outputshape, makepng=False, outdir=outdir)
    
            counter = counter + 1.
            filenames.append(id.ukn.filepath) # Take the name of the images aligned.
        else:
            counter = counter

    return counter, filenames


def performaper(obj, alignedimages, band, refimages):
    """
    We now have many aligned images and a coordinate list, we are ready to do aperture photometry. First we need to set
    the task parameters for all the IRAF tasks that work with the photometry task for reference check out what suzanne
    did here http://www.ctio.noao.edu/noao/content/13-m-smarts-photometric-calibrations-bvri. Setting pyraf parameters
    is somewhat confusing, for more details head here :
    http://www.stsci.edu/institute/software_hardware/pyraf/pyraf_faq#4.6
    """
    # centerpars
    centerpars = iraf.centerpars.getParDict()
    centerpars['calgorithm'].set('centroid')
    centerpars['cbox'].set('10')
    centerpars['cthreshold'].set('0')
    centerpars['minsnratio'].set('1')
    centerpars['cmaxiter'].set('10')
    centerpars['maxshift'].set('5')
    centerpars['clean'].set('no')
    centerpars['rclean'].set('1')
    centerpars['rclip'].set('2')
    centerpars['kclean'].set('3')
    centerpars['mkcenter'].set('yes')

    # datapars
    datapars = iraf.datapars.getParDict()
    datapars['scale'].set('1')
    datapars['fwhmpsf'].set('5')
    datapars['emission'].set('yes')
    datapars['sigma'].set('INDEF')
    datapars['datamin'].set('-500')
    datapars['datamax'].set('30000')
    datapars['noise'].set('poisson')
    datapars['readnoise'].set('6.5')
    datapars['epadu'].set('2.3')
    datapars['exposure'].set('EXPTIME')
    datapars['airmass'].set('SECZ')
    datapars['filter'].set('CCDFLTID')
    datapars['obstime'].set('JD')

    # findpars
    findpars = iraf.findpars.getParDict()
    findpars['threshold'].set('4.0')
    findpars['nsigma'].set('1.5')
    findpars['ratio'].set('1')
    findpars['theta'].set('0')
    findpars['sharplo'].set('0.2')
    findpars['sharphi'].set('1.')
    findpars['roundlo'].set('-1.')
    findpars['roundhi'].set('1.')
    findpars['mkdetections'].set('no')

    # skypars
    skypars = iraf.fitskypars.getParDict()
    skypars['salgorithm'].set('mode')
    skypars['annulus'].set('25.')
    skypars['dannulus'].set('7.')
    skypars['skyvalue'].set('0.')
    skypars['smaxiter'].set('10')
    skypars['sloclip'].set('0.')
    skypars['shiclip'].set('0.')
    skypars['snreject'].set('50')
    skypars['sloreject'].set('3.')
    skypars['shireject'].set('3.')
    skypars['khist'].set('3.')
    skypars['binsize'].set('0.1')
    skypars['smooth'].set('no')
    skypars['rgrow'].set('0.')
    skypars['mksky'].set('no')

    # photpars
    photpars = iraf.photpars.getParDict()
    photpars['weighting'].set('constant')
    photpars['apertures'].set('19')
    photpars['zmag'].set('25.')
    photpars['mkapert'].set('no')
 
    pathalignedimages = ["./aligned/"+x.replace('rccd', 'arccd') for x in alignedimages]    

    for im in pathalignedimages:
        coords = refimages+"/ref_images/"+obj+"/"+band+"/photcoords.lis"
        iraf.phot(im, fwhmpsf=5, datamin=-500., datamax=40000., readnoise=6.5, epadu=2.3, exposur="EXPTIME",
                  airmass="SECZ", obstime="JD", filter="CCDFLTID", annulus=25, dannulu=7, apertur="19", verify='no',
                  output=im+".mag.MC", coords=coords, verbose='no')


def performpsf(obj, alignedimages, band, refimages):
    alignedimages = ['a'+x.split('/')[-1] for x in alignedimages]
    makePSFPhotDF(obj, alignedimages, band, refimages, fwhmthresh=14.0, pickle=True, csv=True)




def avgFWHM(image, coords):
    """
    use the iraf task psfmeasure to compute the average fwhm of a few stars
    INPUT: 
    image: the fits image you want to measure the fwhm for
    coords: a text file listing the x,y pixel coordinates of stars in image to use to compute the fwhm.
    I'd recommend choosing ~10, bright stars, evenly sampled across the frame.
    OUTPUT: a string of the average fwhm also your graphics terminal will open up with a gross iraf graphics terminal
    If you are interested you can look at it to see info like ellipticity and average fwhm of the image.
    If you run this in batch mode and consistently are seeing fwhm values >~ 5, it might mean you need to choose
    different stars in your coords file, or check that the coordinates are accurate.
    """

    # we need to create a text file that has the letter 'q' in it
    # this gets fed to psfmeasure as the graphcur input. it forces iraf to quit out of the terminal
    # this way you can compute the fwhm with out the user having to click on anything
    with open('graphcur', 'w') as f:
        f.write('q')

    # use the iraf task psfmeasure to get the fwhm
    # the output is not machine readable friendly, but basically we are after the last thing it prints out
    # this is the avg fwhm measured
    stdOutput = iraf.psfmeasure(image, imagecur=coords, display='no', graphcur='graphcur', Stdout=1, logfile='')
    #FWHM = stdOutput[-1].split()[-1]

    psfstarlines = stdOutput[4:-2]
    linessplitted = [x.split() for x in psfstarlines]
    FWHM = np.median(np.array([x[3] for x in linessplitted]).astype(float))

    # clean up the graphcur textfile we made
    iraf.delete('graphcur')

    return FWHM


def prepDAOPHOT(fwhm):
    """
    create the input tables daophot needs to do its job
    you need a priori info about the fwhm and detector characteristics to use daophot well
    this program figures that out and prepares the tables for you to do psf photometry on
    INPUT:
    fwhm: the full width half max of an image
    OUTPUT:
    daophot.opt: a text file that lists the parameters for daophot
    photo.opt: a text file that lits the parameters daophot uses to do aperture photometry
    allstar.opt: a text file that allstar uses to do psf photometry
    """
    with open('daophot.opt', 'w') as dao:
        # use the fwhm to determine the fwhm, fit, and psf
        dao.write('FWHM=' + str(fwhm) + '\n')
        dao.write('FIT=' + str(fwhm) + '\n')
        dao.write('PSF=' + str(float(fwhm) * 3.0) + '\n')
        # put in andicam characteristics
        dao.write('READ=6.5\n')
        dao.write('GAIN=2.3\n')
        # set the threshold to 3.5 sigma
        dao.write('TH=3.0\n')
        # use gausian analytic model
        dao.write('AN=1\n')
        # set the low signal as 10 sigma so we dont read into the noise
        dao.write('LOWBAD=10\n')
        # andicam is linear up to this many adu
        dao.write('HIBAD=45000\n')
        # don't ask for user input
        dao.write('WATCH=-2.0\n')
        # the psf should remain the same across the detector
        dao.write('VAR=0\n')

    with open('photo.opt', 'w') as photo:
        # set the photometry radius to the fwhm
        photo.write('A1=' + str(fwhm) + '\n')
        # inner sky radius to 3*fwhm
        photo.write('IS=' + str(float(fwhm) * 3.0) + '\n')
        # outer sky radius to 4*fwhm
        photo.write('OS=' + str(float(fwhm) * 4.0) + '\n')

    with open('allstar.opt', 'w') as allstar:
        # use the fwhm as the psf fitting radius
        allstar.write('fit=' + str(fwhm) + '\n')
        # inner skyr adius to 3*fwhm
        allstar.write('isky=' + str(float(fwhm) * 3.0) + '\n')
        # outer sky radius to 4+fwhm
        allstar.write('osky=' + str(float(fwhm) * 4.0) + '\n')
        # dont ask for user input
        allstar.write('watch=0\n')
        # allow all star to redetermine centroids to give a better quality of fit
        allstar.write('redet=1\n')

    return


def daophotwrap(filename):
    """
    creates a shell script that will call daophot to run find, pick, and psf
    on the image "filename"
    INPUT:
    filename: a string of the filename you want to run daophot on
    OUTPUT:
    daophotGo.sh: a shell script that wraps around daophot
    """
    with open('daophotGo.sh', 'w') as f:
        f.write('#!/bin/tcsh\n')
        f.write('daophot <<__DAOPHOT-END__\n')
        f.write('attatch ' + filename + '\n')
        f.write('find\n')
        f.write('1,1\n')
        f.write(filename + '.coo\n')
        f.write('y\n')
        f.write('phot\n')
        f.write('photo.opt\n')
        f.write('\n')
        f.write(filename + '.coo\n')
        f.write(filename + '.ap\n')
        f.write('pick\n')
        f.write(filename + '.ap\n')
        # pick 20 stars, magnitude limit 20
        f.write('20,20\n')
        f.write(filename + '.lst\n')
        f.write('psf\n')
        f.write(filename + '.ap\n')
        f.write(filename + '.lst\n')
        f.write(filename + '.psf\n')
        f.write('exit\n')
        f.write('__DAOPHOT-END__')
        f.write('\n')

    # give executable permission to the shell script
    os.chmod('daophotGo.sh', 0755)
    return


def allstarWrap(filename):
    """
    creates a shells cript that will call allstar to do psf photometry on the image "filename"
    INPUT:
    filename: a string of the filename you want to do psfphotometry on
    OUTPUT:
    allstarGo.sh: a shell script that wraps around allstar
    """
    with open('allstarGo.sh', 'w') as f:
        f.write('#!/bin/tcsh\n')
        f.write('allstar <<__ALLSTAR-END__\n')
        f.write('\n')
        f.write(filename + '\n')
        f.write(filename + '.psf\n')
        f.write(filename + '.ap\n')
        f.write(filename + '.als\n')
        f.write(filename + '.rej\n')
        f.write(filename + '.sub\n')
        f.write('\n')
        f.write('__ALLSTAR-END__\n')
        f.write('\n')
        
    # give executable permissions to shell script
    os.chmod('allstarGo.sh', 0777)
    pass


def findSources(alsfile, coords, tol):
    """
    look in the output of allstars, alsfile, and find the entries of particular stars,
    whose x,y pixel coordinate positions are given in the textfile, coords. 
    Uses a vectorized nearest neighbor search
    this function was made to be called by others like makePSFPhotDF below, to help return a photometry table for a
    target by itself, this function only returns photometry information for one target
    INPUT:
    ------
    alsfile: the file name of the allstar output. typically this has a .als multi-extension
    coords: the file name of a a text file that lists the coordinates of the stars you are interested in
    IMPORTANT: The program is currently hard coded to expect a third colum with junk values. This is a legacy aspect.
    When using imexamine to print out coordinates, a third column with count number was also printed. this program is
    hard coded to read in the third column, and subsequently ignore it. It is vital that photcoods contains a third
    column with junk values to work, oherwise the program will crash. If you follow the steps in the wiki and use
    imexamine and press the 'x' key to to generate this list, that will already be the case, and you will be fine.
    Consult the wiki for further details on how to use this.
    tol: tolerance distance in pixels. If we cannot find any stars with distances below the tol
    we don't remember count its daophot data and instead return nans
    OUTPUT:
    -------
    photBios: a python dictionary with photometry info for these sources 
    """
    # first read in the psf photometry output file to a pandas dataframe
    alsPhot = pd.read_table(alsfile, skiprows=[0, 1, 2],
                            names=['ID', 'x', 'y', 'mag', 'err', 'sky', 'nit', 'chi', 'sharp'], sep='\s+')

    # now read the coords file into a dataframe
    coords = pd.read_table(coords, names=['x', 'y', 'junk'], sep='\s+')

    # count up how many sources are in the coords list
    numSources = len(coords)

    # prepare a numpy array to hold the dataframe indicies for the stars in the coords list
    starIlocs = np.zeros(numSources)

    # loop over the stars in the coords dataframe. for every entry, find its nearest neighbor
    # in the alsPhot dataframe
    for i in xrange(numSources):
        coordsView = coords.iloc[i]
        x = coordsView.x
        y = coordsView.y
        # calculate the distances between this star in coordsView and all the stars in the als dataframe
        distances = np.sqrt(np.square(alsPhot.x - x) + np.square(alsPhot.y - y))
        # if the smallest distance is within the tolerance range, keep this star's index number in
        # the als data frame in the starIloc array
        if distances.min() < tol:
            print long(np.sqrt(np.square(alsPhot.x - x) + np.square(alsPhot.y - y)).argmin())
            starIlocs[i] = np.sqrt(np.square(alsPhot.x - x) + np.square(alsPhot.y - y)).argmin()
        # otherwise just put a nan there
        else:
            starIlocs[i] = np.nan

    # prepare a python dict to keep the data and meta data from the als file on the stars we found
    photBios = {}

    # loop over the stars we found, look them up in the alsPhot dataframe, and save their info
    # in the photBios dict.
    for i in xrange(numSources):
        # if we were able to find a star with distance within tol save its alsPhot data in photBios
        if np.isfinite(starIlocs[i]):
            alsPhotView = alsPhot.iloc[long(starIlocs[i])]
            for j in alsPhotView.index:
                photBios['s' + str(i) + '_' + j] = alsPhotView[j]
        # otherwise just fill everything with nans
        else:
            for j in alsPhot.columns:
                photBios['s' + str(i) + '_' + j] = np.nan

    return photBios


def makePSFPhotDF(obj, alignedimages, band, refimages, fwhmthresh=8.0, pickle=True, csv=True):
    """
    pathalignedimages = ["./aligned/"+x.replace('rccd', 'arccd') for x in alignedimages]

    for im in pathalignedimages:
        coords = refimages+"/ref_images/"+obj+"/"+band+"/photcoords.lis"
        iraf.phot(im, fwhmpsf=5, datamin=-500., datamax=30000., readnoise=6.5, epadu=2.3, exposur="EXPTIME",
                  airmass="SECZ", obstime="JD", filter="CCDFLTID", annulus=25, dannulu=7, apertur="9", verify='no',
                  output=im+".mag.MC", coords=coords, verbose='no')
    """
    row_list = []
    for im in alignedimages:
        photcoords = refimages+"/ref_images/"+obj+"/"+band+"/photcoords.lis"
        psfcoords = refimages+"/ref_images/"+obj+"/"+band+"/psfcoords.lis"
        print 'makePSFPhotDF is running for ' + im + ' in ' + band + '-band'
        print(os.getcwd())

        if len(im) < 31:
            hdulist = fits.open(im)
            header = hdulist[0].header
            hdulist.close()
            # get some useful info about the image from its header
            # occasionally the header info is corrupted, so we add error handling
            # just put in np.nan if there is a key error
            try:
                JD = header['JD']
            except KeyError:
                JD = np.nan
            try:
                airmass = header['SECZ']
            except KeyError:
                airmass = np.nan
            try:
                filename = header['FILENAME']
            except KeyError:
                filename = np.nan
            try:
                obsdate = header['DATE-OBS']
            except KeyError:
                obsdate = np.nan
            try:
                exptime = header['EXPTIME']
            except KeyError:
                exptime = np.nan

            # use the psfcoords.lis file to measure the fwhm of im
            fwhm = avgFWHM(im, psfcoords)
            print 'FWHM:', fwhm
            # if the fwhm is lower than the fwhmthresh, proceed with the reduction
            if float(fwhm) < fwhmthresh and float(fwhm) > 1.0:
                # create the .opt files daophot and allstars use
                prepDAOPHOT(fwhm)
                # run daophot on this image
                daophotwrap(im)
                os.system('./daophotGo.sh')

                # sometimes daophot cannot create a psf.
                # in this case, it creates an empty .psf file
                # we must check the length of the .psf file here, and proceed if it succeeded
                psf_line_num = int(subprocess.check_output(['wc', '-l', im + '.psf']).split()[0])

                if psf_line_num > 0:
                    # run allstar on this image
                    allstarWrap(im)
                    os.system('./allstarGo.sh')

                    # all star might have encoutered a problem. make sure the .als file exists
                    # before we try to read into it
                    if len(glob.glob(im + '.als')) == 1:
                        # get the psfphotoemtry data of the stars in the photcoords.lis file
                        photBios = findSources(im + '.als', photcoords, float(fwhm))
                        # add the header info we found earlier to the photBios dict
                        # put in error handling, sometimes the fits headers will have bad values that cant be converted to floats
                        try:
                            photBios['JD'] = float(JD)
                        except ValueError:
                            photBios['JD'] = np.nan
                        try:
                            photBios['airmass'] = float(airmass)
                        except ValueError:
                            photBios['airmass'] = np.nan
                        try:
                            photBios['fwhm'] = float(fwhm)
                        except ValueError:
                            photBios['fwhm'] = np.nan
                        try:
                            photBios['exptime'] = float(exptime)
                        except ValueError:
                            photBios['exptime'] = np.nan
                        photBios['filename'] = filename
                        photBios['obsdate'] = obsdate
                        photBios['align_filename'] = im
                        # add this dict to the running list
                        row_list.append(photBios)
                    else:
                        print "could not find an .als file. " + im + " wont be added to photometry log"

                # if the psf file is empty print a message to the stdout saying we're skippnig psfphotometry
                else:
                    print "there is no psf for " + str(im)
                    print "skipping psf photometry for " + str(im)

            # if the fwhm was too big, skip the photometry steps and print out an error message
            else:
                print "the fwhm for " + str(im) + " is " + str(fwhm) + " which is bigger than the threshold " + str(
                    fwhmthresh)
                print "skipping photometry for " + str(im)

        else:
            print "the file " + im + " has " + str(len(im)) + " characters"
            print "daophot will only accept files with 30 characters or less"
            print im + " will not be reduced"

        # clean up the intermediate files before going on to the next file
        # flag = os.system('rm *.opt *.sh *.ap *.coo *.lst *.nei *s.fits')
        flag = os.system('rm *.sh *.ap *.coo *.lst *.nei')
    # shove everything ito a pandas data frame and return it
    photdf = pd.DataFrame(row_list)
    csvfilename = obj+band+'photdf.csv'
    pklfilename = obj+band+'photdf.pkl'
    if csv:
        photdf.to_csv(csvfilename)
    if pickle:
        photdf.to_pickle(pklfilename)

    return photdf


def mineALS(filelist, fltid, photcoords='photcoords.lis', pickle=True, csv=True):
    '''
    this function loops over all the .als files that are produced from psf photometry
    and extracts data for particular stars which are specified by their pixel coordinates
    in the input file photcoords

    INPUT:
    filelist: a python list containing the file names of the .als files you want to work on
    fltid: a string specifying the filter that this data were taken in. 
	all the data need to be from the same filter. example values are 'B','V','I'.
	this is important because it tells the program the path to the fits image that made
	these .als files. 
    photcoords: a text file that has the x y coordinates of sources you are intereted in looking up
	in the .als files. it is set to 'photcoords.lis' by default.
	if you want to create different dataframes with which use different stars, you can specify 
	different photcoords files

    OUTPUT:
    photdf: a pandas dataframe. the columns are from the fits header and the column headers from the .als file
	each row is from a different .als file.
    photdf_temp.pkl: if the keyword pickle is set to True, a pickled version of the dataframe is saved to disk
	with the file name 'photdf.pkl'
    photdf_temp.csv: if the keyword csv is set to True, the dataframe is saved to disk as a comma delimited file.
    '''

    row_list = []

    # loop over the .als files. recall that these are the output files from psf photometry
    for f in filelist:
        # we have to read the fwhm from the head of the als file.
        # it is the last item in the second column. we can use the python module linecache
        # read the second line, strip the leading and ending white space, split the strings, and grab the last one
        fwhm = linecache.getline(f, 2).strip().split()[-1]

        # get the psfphotoemtry data of the stars in the photcoords file
        # recall that findSources will look at the photcoords file for pixel coordinates
        # and then look up these coordinates in the .als file. if it finds a close match
        # it will pull out all the psf photometry data at that coordinate position from the .als file
        photBios = findSources(f, photcoords, float(fwhm))

        # now that we are storing important observation data in the python dictionary photBios
        # lets also add in the fwhm and name of the als file we are working on
        try:
            photBios['fwhm'] = float(fwhm)
        except ValueError:
            photBios['fwhm'] = np.nan

        photBios['align_filename'] = f[:-4]
        # add this dict to the running list
        row_list.append(photBios)

        # now we have to open up the header of the fits image that produced this als file
        # so we can read important data like observation time and airmass
        hdulist = fits.open('../' + fltid.upper() + '_align/' + f[:-4])
        header = hdulist[0].header
        hdulist.close()

        # now that he have the header stored in memory, we can lookup important keyword value pairs
        # we will also do some error handling, incase there is bad or corrupted header data
        # in the case of bad values, we will put an np.nan

        # grab julian date
        try:
            JD = header['JD']
        except KeyError:
            JD = np.nan
        try:
            photBios['JD'] = float(JD)
        except ValueError:
            photBios['JD'] = np.nan

        # grab airmass
        try:
            airmass = header['SECZ']
        except KeyError:
            airmass = np.nan
        try:
            photBios['airmass'] = float(airmass)
        except ValueError:
            photBios['airmass'] = np.nan

        # grab exposure time
        try:
            exptime = header['EXPTIME']
        except KeyError:
            exptime = np.nan
        try:
            photBios['exptime'] = float(exptime)
        except ValueError:
            photBios['exptime'] = np.nan

        # grab filename. this is the name as it was written by the andicam computer
        # so it will be of the form ccdYYMMDD.xxxx
        try:
            photBios['filename'] = header['FILENAME']
        except KeyError:
            photBios['filename'] = np.nan
        # grab observation date
        try:
            photBios['obsdate'] = header['DATE-OBS']
        except KeyError:
            photBios['obsdate'] = np.nan

    # now that we are done looping over the als files,
    # lets stick our list into a pandas data frame
    photdf = pd.DataFrame(row_list)

    # if the keywords were set, save the data frame
    # by either pickling it or saving it as a csv
    if pickle:
        photdf.to_pickle('mineALS.pkl')
    if csv:
        photdf.to_csv('mineALS.csv')

    return photdf
