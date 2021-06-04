"""parse.py

Module
------
parse.py - a module containing general utilities used by the SDSSRefs class and app.

"""
import os
import pathlib

import pandas as pd

from SDSSRefs import config
from LCExtract import config as LCconfig


class ParseDir:
    table = pd.DataFrame()

    @classmethod
    def getSubDir(cls, dirPath: str) -> list:
        if not os.path.isdir(dirPath):
            print(f'Subdirectory \"{dirPath}\" does not exist.')
            exit(1)
        print(f'Found Subdirectory \"{dirPath}\".')
        return [name for name in os.listdir(dirPath) if os.path.isdir(os.path.join(dirPath, name))]

    @classmethod
    def getFiles(cls, subDir: str) -> list:
        return [name for name in os.listdir(subDir) if os.path.isfile(os.path.join(subDir, name))]

    @classmethod
    def getTable(cls, fileName: str):
        cls.table = pd.read_csv(fileName, header=0)
        return cls.table

    @classmethod
    def pruneMultiple(cls) -> bool:  # DO NOT USE
        oids = cls.table['oid'].unique()
        if len(oids) > 1:
            print('Multiple objects found. First object used.')
            cls.table = cls.table[cls.table['oid'] == oids[0]]
        return len(oids) > 1

    @classmethod
    def saveTable(cls, fileName: str):
        cls.table.to_csv(fileName, index=False, header=True)


def parseLC(fullLCTable: pd.DataFrame):
    # fullLCTable = pd.DataFrame()
    basePath = 'data/LC'
    fullFN = f'{basePath}/fullLCdata.csv'

    # if full LC data file already exists
    if os.path.isfile(fullFN):
        # then load into dataframe
        fullLCTable = pd.read_csv(fullFN, header=0)
    else:
        # otherwise, read in and consolidate from all individual LC files in mag subdirectories
        LCconfig.filterSelection = 'gr'
        archive = LCconfig.archives['z']
        # dirInLC = ParseDir.getSubDir('data/LC')
        subDir_Count = 0
        file_Count = 0
        LCdirs = ParseDir.getSubDir(basePath)
        dirTotal = len(LCdirs)
        for subDir in LCdirs:
            subDir_Count += 1
            if config.verbose:
                print(f'\nParsing subdirectory {subDir}, {subDir_Count} of {dirTotal}:')
            filePath = basePath + os.sep + subDir
            LCfiles = ParseDir.getFiles(filePath)
            for LCFile in LCfiles:

                if 'gr' in LCFile:
                    filePathName = filePath + os.sep + LCFile
                    if config.verbose == 'full':
                        print(f'Opening {filePathName}')
                    file_Count += 1
                    if len(fullLCTable):
                        fullLCTable = fullLCTable.append(ParseDir.getTable(filePathName))
                    else:
                        fullLCTable = ParseDir.getTable(filePathName)

        print(f'\nSubdirectories count: {subDir_Count}')
        print(f'\nTotal files count: {file_Count}')

        # save file for next time
        fullLCTable.to_csv(f'{basePath}/fullLCdata.csv', index=False, header=True)

    return fullLCTable


def parseSDSS(fullSDSSTable: pd.DataFrame):
    # fullLCTable = pd.DataFrame()
    basePath = 'data/SDSS'
    fullFN = f'{basePath}/fullSDSSdata.csv'

    # if full SDSS data file already exists
    if os.path.isfile(fullFN):
        # then load into dataframe
        fullSDSSTable = pd.read_csv(fullFN, header=0, sep=' ')
    else:
        # otherwise, read in and consolidate from all individual LC files in mag subdirectories
        LCconfig.filterSelection = 'gr'
        archive = LCconfig.archives['z']
        # dirInLC = ParseDir.getSubDir('data/LC')
        file_Count = 0

        filePath = basePath
        SDSSfiles = ParseDir.getFiles(filePath)
        for SDSSFile in SDSSfiles:

            filePathName = filePath + os.sep + SDSSFile
            if config.verbose == 'minimal':
                print(f'Opening {filePathName}')
            file_Count += 1
            if len(fullSDSSTable):
                fullSDSSTable = fullSDSSTable.append(ParseDir.getTable(filePathName))
            else:
                fullSDSSTable = ParseDir.getTable(filePathName)

        print(f'\nTotal files count: {file_Count}')

        # save file for next time
        fullSDSSTable.to_csv(fullFN, index=False, header=True)

    return fullSDSSTable