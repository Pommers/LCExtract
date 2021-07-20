"""
Module
------
images.py: Wrapper for image retrieval module

Summary
-------
    Allows selection of file or manual input.

    Allows filter selection.

    Iterates through required entries

    retrieving available data

    summarising and plotting output
"""

from datetime import datetime
from matplotlib import pyplot as plt, gridspec
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

# Self-authored package modules for inclusion
from LCExtract import config
from LCExtract.SDSS import SDSSdata
from LCExtract.config import archives
from LCExtract.dataretrieve import AstroObjectClass, AODataClass, addColumns
from LCExtract.entry import getObjects, setFilterUsage, setArchiveUsage
from LCExtract.ztf import refZTFobj, getZTFOidOffsetStats


def startup():
    print()
    print('Image retrieval module\n'
          '----------------------')
    print()


def showBaseImages(imageDataList, title, subtitle):
    rows = len(imageDataList)
    height = 3.5 * rows + 1

    fig = plt.figure(num=2, figsize=(11, height), dpi=300)
    fig.suptitle(title, fontsize=16)
    spec = gridspec.GridSpec(ncols=3, nrows=rows, figure=fig)
    # baseAxes = (1 * 100) + (len(hdu) * 10) + 1

    for r, hduR in enumerate(imageDataList):
        for c, hd in enumerate(hduR[0]):
            size = hduR[1][c]
            imageThumb(fig, r, c, hd, spec, f'{subtitle[r]} ({size} arcmin)')

    fig.tight_layout(h_pad=4, w_pad=6, rect=(0.07, 0.04, 0.98, 0.99))
    fig.savefig(f'data/plots/{title}_ref_img.png')
    if config.plotToScreen:
        fig.show()
    plt.close(fig=2)


def showOutlierImages(hduThumbs, title, subtitle):
    rows = int((len(hduThumbs)-1)/3)+1
    height = 3.5 * rows + 1

    fig = plt.figure(num=3, figsize=(11, height), dpi=300)
    fig.suptitle(title, fontsize=16)
    spec = gridspec.GridSpec(ncols=3, nrows=rows, figure=fig)
    # baseAxes = (1 * 100) + (len(hdu) * 10) + 1

    for i, ((hd, size), stday) in enumerate(hduThumbs):
        r = int(i/3)
        c = int(i % 3)
        imageThumb(fig, r, c, hd[0], spec, f'{subtitle}{stday} ({size[0]} arcmin)')

    fig.tight_layout(h_pad=4, w_pad=6, rect=(0.07, 0.04, 0.98, 0.99))
    fig.savefig(f'data/plots/{title}_outlier_img.png')
    if config.plotToScreen:
        fig.show()
    plt.close(fig=3)


def imageThumb(fig, r, c, hd, spec, subtitle):
    wcs = WCS(hd[0].header)
    raPix = hd[0].header['NAXIS1']
    raMin = wcs.wcs_pix2world(((-raPix, 0),), 0)[0][0]
    raMax = wcs.wcs_pix2world(((raPix, 0),), 0)[0][0]
    decPix = hd[0].header['NAXIS1']
    decMin = wcs.wcs_pix2world(((0, -decPix),), 0)[0][1]
    decMax = wcs.wcs_pix2world(((0, decPix),), 0)[0][1]
    # axisNum = baseAxes + c
    ax = fig.add_subplot(spec[r, c], projection=wcs)
    ax.set_title(subtitle, fontsize=12)
    ax.imshow(hd[0].data, origin='lower', cmap=plt.cm.gist_yarg)
    ax.grid(color='gray', ls='dotted')
    ra = ax.coords['ra']
    if raMin < raMax:
        ax.invert_xaxis()
    ra.set_axislabel('RA (deg)', minpad=0.9)
    ra.set_major_formatter('dd:mm:ss')
    ra.set_ticks(number=4)
    ra.set_ticklabel(exclude_overlapping=True, fontsize=8)
    # ra.set_separator(('d', "'", '"'))
    dec = ax.coords['dec']
    if decMin > decMax:
        ax.invert_yaxis()
    # if not c:
    dec.set_axislabel('Dec (deg)', minpad=-0.8)
    # else:
    #     dec.set_auto_axislabel(False)
    dec.set_ticks(number=4)
    dec.set_ticklabel(exclude_overlapping=True, fontsize=8, ha='center')


def image():
    config.filterSelection = 'gr'
    archiveList = 'z'

    archiveLCData = {}
    outliersFound = {}
    cointegration = {}
    refAOList = {}
    meanMagOffset = {}
    startup()
    objectsList = getObjects()
    # archiveList = setArchiveUsage()
    # setFilterUsage()
    addColumns(objectsList, config.filterSelection)

    for i in objectsList:
        AO = AstroObjectClass(i['Name'], i['RA'], i['DEC'], i['Description'])
        print(f"\n-------------------------------------------------------------\n"
              f"Object name: {AO.objectName} - summary statistics")
        for a in archiveList:
            archiveLCData[a] = AODataClass(AO, archives[a])
            if archiveLCData[a].getData():
                AO.incrPlotRows()
                AO.append_filtersToPlot(archiveLCData[a].filtersReturned)
                archiveLCData[a].objectOutput()
                if archives[a].name == 'ZTF':
                    # Save data
                    # archiveLCData[a].objectStatSave(i)
                    # check if archive data has data outside 3 sigma
                    outliersFound[a] = archiveLCData[a].outliersExist(archives[a], config.threshold)
                    # Get first row images - change to ref image?
                    for f in archiveLCData[a].filtersReturned:
                        refImageData = archiveLCData[a].getRefImage(f)
                        medianRow = archiveLCData[a].getMedianRow(f)
                        medianImageData = archiveLCData[a].getImages(medianRow, config.imgSize)

                        titles = [f'Reference image', f'Median image']
                        showBaseImages([refImageData, medianImageData], f'{AO.objectName} ({f}-band)', titles)

                        if outliersFound[a][f]:
                            objectsList[i.index]['ZTFOutliers'] = True
                            outliers = archiveLCData[a].table[
                                (archiveLCData[a].table['outlier'] != 'inside') &
                                (archiveLCData[a].table['filterID'] == f)
                            ]
                            if config.checkOutliers:
                                outlierImageData = []
                                # for each of the outlier datapoints in the LC sample
                                for i1, s in outliers.iterrows():
                                    outlierImageData.append([archiveLCData[a].getImages(i1, [config.imgSize[-1]]), f"{s['filefracday']}"])

                                    if False:  # Don't do this at the moment
                                        refs = refZTFobj(AO.pos,
                                                         s[archives[a].timeField],
                                                         s['filefracday'],
                                                         s[archives[a].magField])
                                        ZTFoffset, ZTFsd = getZTFOidOffsetStats(refs)
                                        print(f'{AO.objectName} - Outlier MJD{s[archives[a].timeField]}: '
                                              f'Offset = {ZTFoffset:-8.5f} '
                                              f'SD = {ZTFsd:-8.5f}')
                                showOutlierImages(outlierImageData, f'{AO.objectName} Outliers ({f}-band)', '')
                                sdssRefAOs = SDSSdata(archiveLCData[a].median, AO.pos)
                                if sdssRefAOs.getSamples():
                                    pass
                            # return a df of image data frames (idf) to check
                            # for idf in archiveLCData[a].outlierIDFList:
                            # refAOList[idf] = archiveLCData[a].getRefObjects(idf)
                            # get a df of ref objects in idf (with SDSS mag)
                            # meanMagOffset[idf] = archiveLCData[a].getIDFmagOffset(refAOList[idf])
                            # get mean mag offset for ref objects in IDF
            else:
                print(f'No data available or retrieved from {archives[a].name}')
                print()
        filtersToPlot = AO.get_filtersInAO()
        if AO.plotRowsTrue():
            fig, ax = AO.preparePlot(filtersToPlot)
            for a in archiveList:
                archiveLCData[a].plot(fig, ax, archives[a], filtersToPlot)
            AO.finalisePlot(fig, ax, filtersToPlot)

    # datetime object containing current date and time
    now = datetime.now()

    dataFileN = 'LCStats' + now.strftime('%Y%m%d_%H%M%S') + '.csv'

    with open(f'data/{dataFileN}', 'w') as f:
        objectsList.write(f, format='csv')
