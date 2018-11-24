import SMARTSIR as IR
IR.pd.set_option('display.width', 1000)


# Input parameters
daterange = ['2018-07-12', '2018-10-24']
source = ['J1820'] # 'J1820', 
expr = 'OWNER=="YALE-03A-0001"'
pathtorefimages = '/Volumes/HighEnergy1/dailylc'
inputdirs = ['/Volumes/HighEnergy1/archive/SMARTS13m/IR/J1820']
output_ir = 'MAXIJ1820_ir.list'
inputdirs_flat_ir = ['/Volumes/HighEnergy1/archive/SMARTS13m/IRFLATS']
daterange_flat_ir = ['2018-07-12', '2018-10-24']
ref_ir     = {'J': './reference/J/clean1.fits', 
              'H': './reference/H/clean1.fits',
              'K': './reference/K/final.fits'}
obsstrategy = {'J': 7, 'H': 7, 'K':7}



###
### Optical Pipeline
###

#IR.filefindcp(daterange=daterange, inputdir=inputdirs, source=source, band=['J', 'H'])

#IR.listobs(expr)

#IR.filterobslist(obslist='obsmenu.list', source=source, daterange=daterange)

#IR.flatList_ir(inputdirs_flat_ir, daterange_flat_ir)

#IR.dithsummary()

#IR.redwrap(ref_ir, obsstrategy, source)
#IR.imconstrcut(ref_ir, obsstrategy, source, 1)
#IR.imagevar(source[0], 'K')

IR.performphot(source)

#IR.makerawtables(pathtorefimages)


