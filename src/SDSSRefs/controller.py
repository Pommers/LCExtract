"""controller.py

Module
------
controller.py

Description
-----------
Module to control the selection of 'pseudo-random' stars from SDSS database for use in definition of transient facility
intrinsic scatter / error for different magnitudes.

"""
from numpy import arange

from LCExtract import config as LCconfig
from LCExtract.dataretrieve import AstroObjectClass, AODataClass
from SDSSRefs import config
from SDSSRefs.SDSS import SDSSdata, NoLCList, Sky
from SDSSRefs.outputStats import Stats
from SDSSRefs.resample import resample


def controller():
    statsOfStats = Stats()
    # set config for samples

    LCconfig.filterSelection = 'gr'
    archive = LCconfig.archives['z']
    fileInitialise = True
    if config.SDSSDataRefresh:
        for mag in arange(config.minMag, config.maxMag, config.stepMag):
            if config.verbose:
                print(f'Sampling magnitude {mag} from SDSS.')
            sdssDataHolder = SDSSdata(mag)
            status, objectsList = sdssDataHolder.getSamples()
            if status:
                if config.verbose:
                    print(f'Sample completed successfully.')
                sdssDataHolder.addColumns()
                sdssDataHolder.createFile(overwrite=True)
                noLClist = NoLCList()
                for i in objectsList:
                    if config.verbose == 'full':
                        print(f'Querying data for object {i.index + 1} ({i["Name"]}). ', end='')
                    objectHolder = AstroObjectClass(i['Name'], i['RA'], i['DEC'], i['Description'])
                    archiveDataHolder = AODataClass(objectHolder, archive)
                    # if no file exists or refresh set true get archive data for object
                    if not archiveDataHolder.fileExists(mag) or config.LCDataRefresh:
                        if not noLClist.inList(i["Name"]):
                            if archiveDataHolder.getData():
                                if config.verbose == 'full':
                                    print(f'Object {i.index + 1} lightcurve data found, incl. '
                                          f'{archiveDataHolder.filtersReturned} filter(s)')
                                if config.LCDataSave:
                                    # save data if requested
                                    archiveDataHolder.saveTable(mag)
                                archiveDataHolder.objectStatSave(i)
                            else:
                                noLClist.add(i["Name"], mag)
                                if config.verbose == 'full':
                                    print(f'No lightcurve data exists in archive for object {i.index + 1} - '
                                          f'Added {i["Name"]} to list.')
                        else:
                            if config.verbose == 'full':
                                print(f'\r{" ":66}\r', end='')
                                print(f'Object {i.index + 1} found in no lightcurve list.')
                    else:
                        # retrieve light curve data from file
                        if config.verbose == 'full':
                            print(f'\r{" ":66}\r', end='')
                            print(f'Object {i.index + 1} lightcurve data loaded from previous file.')
                        archiveDataHolder.loadTable()
                        archiveDataHolder.objectStatSave(i)

                fileInitialise = sdssDataHolder.saveSamples(fileInitialise)

                for f in 'gr':
                    statsOfStats.add_row(sdssDataHolder.getStats(f))
            print()

    resample()
