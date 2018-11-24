import SMARTSOpt_J1820 as Opt
Opt.pd.set_option('display.width', 1000)


# Input parameters
daterange = ['2018-10-23', '2018-10-31'] # ['2018-03-16', '2018-07-01']
source = ['J1820']
expr = 'OWNER=="YALE-03A-0001"'
pathtorefimages = '/Volumes/HighEnergy1/dailylc'
inputdirs = ['/Volumes/HighEnergy1/archive/SMARTS13m/CCD/J1820']
output_opt = 'J1535-571_opt.list'



###
### Optical Pipeline
###

Opt.filefindcp(inputdir=inputdirs, output=output_opt, source=source, band=['B','R','I','V'], daterange=daterange, workingdir='./')

#Opt.listobs(expr)

#Opt.filterobslist(obslist='obsmenu.list', source=source, daterange=daterange)

#naligned, filenames = Opt.alignprocessedwrapper('filteredobsmenu.csv', pathtorefimages)

#Opt.performphot('flagalignedfilteredobsmenu.csv', pathtorefimages=pathtorefimages)

#Opt.makerawtables(pathtorefimages)
