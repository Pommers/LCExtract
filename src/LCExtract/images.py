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
import copy

import numpy as np
from astropy.visualization import ImageNormalize, SqrtStretch
from matplotlib import pyplot as plt, gridspec
from astropy.wcs import WCS
from astropy.stats import SigmaClip, sigma_clipped_stats
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from photutils import Background2D, MedianBackground

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
            imageThumb(fig, r, c, hd, spec, f"{subtitle[r]} ({size}')")
            if c == 1:  # if this is central image...
                estimateBackground(hd, subtitle[r], out=True)

    fig.tight_layout(h_pad=4, w_pad=6, rect=(0.07, 0.04, 0.98, 0.99))
    fig.savefig(f'data/plots/{title}_ref_img.png')
    if config.plotToScreen:
        fig.show()
    plt.close(fig=2)


def showOutlierImages(hduThumbs, title, subtitle):
    rows = int((len(hduThumbs) - 1) / 3) + 1
    height = 3.5 * rows + 1

    fig = plt.figure(num=3, figsize=(11, height), dpi=300)
    fig.suptitle(title, fontsize=16)
    spec = gridspec.GridSpec(ncols=3, nrows=rows, figure=fig)
    # baseAxes = (1 * 100) + (len(hdu) * 10) + 1

    for i, ((hd, size), stday) in enumerate(hduThumbs):
        r = int(i / 3)
        c = int(i % 3)
        imageThumb(fig, r, c, hd[0], spec, f"{subtitle}{stday} ({size[0]}')")

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
    # extract (or calc) sd of image thumbnail for display
    # sd hd[0].header['GSTDDEV']
    sd = np.nanstd(np.hstack(hd[0].data))
    try:
        see = hd[0].header['SEEING']
        see_text = f', FWHM={see:<.2f} pix'
    except KeyError:
        see_text = ''

    ax = fig.add_subplot(spec[r, c], projection=wcs)
    ax.set_title(f"{subtitle}\n\u03C3={sd:<.1f}{see_text}", fontsize=12)
    norm = ImageNormalize(stretch=SqrtStretch())
    ax.imshow(hd[0].data, origin='lower', cmap=plt.cm.plasma, norm=norm)  # was gist_yarg, inverted grayscale
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


def saturatedPixelFlag(imageData):
    headerData = imageData[0].header
    maxPixValue = max(np.hstack(imageData[0].data))
    print(f"Infobits:{headerData['INFOBITS']:>016b}, Max: {maxPixValue}")
    return False


def estimateBackground(imageData, id, out=False):
    sigma_clip = SigmaClip(sigma=3., maxiters=10)
    bkg_estimator = MedianBackground()
    bkg = Background2D(imageData[0].data, box_size=(12, 12), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

    if out:
        mean, median, std = sigma_clipped_stats(imageData[0].data, sigma=3.0)
        meanBkg, medianBkg, stdBkg = sigma_clipped_stats(bkg.background[0].data, sigma=3.0)

        print(f'{id}, Bkg Med     : {bkg.background_median}')
        print(f'{id}, Bkg RMS Med : {bkg.background_rms_median}')
        print(f'{id}, sig clip mean - img: {mean:>7.2f}, bkg: {meanBkg:>7.2f}, '
              f'ratio: {mean / meanBkg:>7.4f}')
        print(f'{id}, sig clip med. - img: {median:>7.2f}, bkg: {medianBkg:>7.2f}, '
              f'ratio: {median / medianBkg:>7.4f}')
        print(f'{id}, sig clip std. - img: {std:>7.2f}, bkg: {stdBkg:>7.2f}, '
              f'ratio: {std / stdBkg:>7.4f}')
        print()

    return bkg


def showBkgImages(hdu, bkg, dif, title):
    rows = 1
    height = 3.5 * rows + 1

    fig = plt.figure(num=4, figsize=(11, height), dpi=300)
    fig.suptitle(title, fontsize=16)
    spec = gridspec.GridSpec(ncols=3, nrows=rows, figure=fig)
    # baseAxes = (1 * 100) + (len(hdu) * 10) + 1

    imageThumb(fig, 0, 0, hdu, spec, f"Image")
    imageThumb(fig, 0, 1, bkg, spec, f"Background")
    imageThumb(fig, 0, 2, dif, spec, f"Subtraction")

    fig.tight_layout(h_pad=4, w_pad=6, rect=(0.07, 0.04, 0.98, 0.99))
    fig.savefig(f'data/plots/{title}_background_img.png')
    if config.plotToScreen:
        fig.show()
    plt.close(fig=4)


def image():
    config.filterSelection = 'gr'
    archiveList = 'z'

    archiveLCData = {}
    outliersFound = {}
    consecutive = {}
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
            # get lightcurve data from archive for object
            if archiveLCData[a].getData():
                AO.incrPlotRows()  # adjust number of rows to include in plot
                AO.append_filtersToPlot(archiveLCData[a].filtersReturned)  # adjust filters returned as appropriate
                archiveLCData[a].objectOutput()  # output summary stats for object to console
                if archives[a].name == 'ZTF':  # the rest is only for ZTF archive
                    # Save data
                    # archiveLCData[a].objectStatSave(i)
                    # check if archive data has data outside 3 sigma
                    outliersFound[a] = archiveLCData[a].outliersExist(archives[a], config.threshold)
                    # Get first row images - change to ref image?
                    for f in archiveLCData[a].filtersReturned:
                        refImageData = archiveLCData[a].getRefImage(f)
                        imageExcessSaturation = False
                        while True:
                            medianRow = archiveLCData[a].getMedianRow(f)
                            medianImageData = archiveLCData[a].getImages(medianRow, config.imgSize)
                            for c, medianImageHDU in enumerate(medianImageData[0]):
                                print(f'{archiveLCData[a].AO.objectName}... Filter:{f}, Median image:{c}, ', end='')
                                imageExcessSaturation |= saturatedPixelFlag(medianImageHDU)
                            if not imageExcessSaturation:
                                break

                        titles = [f'Reference image', f'Median image']
                        showBaseImages([refImageData, medianImageData], f'{AO.objectName} ({f}-band)', titles)

                        if outliersFound[a][f]:
                            objectsList[i.index]['ZTFOutliers'] = True
                            outliers = archiveLCData[a].table[
                                (archiveLCData[a].table['outlier'] != 'inside') &
                                (archiveLCData[a].table['filterID'] == f)
                                ]
                            inliers = archiveLCData[a].table[
                                (archiveLCData[a].table['outlier'] == 'inside') &
                                (archiveLCData[a].table['filterID'] == f)
                                ].nsmallest(9, 'magmedoffset')
                            consecutive[f] = archiveLCData[a].areOutliersConsecutive(outliers)
                            c = consecutive[f]
                            if c['response']:
                                print(f'Consecutive detections in {f}-band: {c["maxCount"]} max.  '
                                      f'Total {c["runs"]} run(s).')
                            else:
                                print(f'No {f}-band consecutive detections')
                            if config.checkOutliers:
                                inlierImageData = []
                                # for 9 of the inlier datapoints in the LC sample
                                for i1, s in inliers.iterrows():
                                    inlierImageData.append(
                                        [archiveLCData[a].getImages(i1, [config.imgSize[-1]]), f"{s['filefracday']}"])

                                showOutlierImages(inlierImageData,
                                                  f'{AO.objectName} 9 Inliers '
                                                  f'({f}-band)', '')

                                outlierImageData = []
                                # for each of the outlier datapoints in the LC sample
                                for i1, s in outliers.iterrows():
                                    outlierImageData.append(
                                        [archiveLCData[a].getImages(i1, [config.imgSize[-1]]), f"{s['filefracday']}"])

                                showOutlierImages(outlierImageData,
                                                  f'{AO.objectName} Outliers '
                                                  f'({f}-band)', '')

                                # Now estimate background levels on outliers (compare with ref and median?)
                                # Maybe use 1' outlier images in order to estimate object positions

                                imgSizeIndx = 1  # should be last image ~ 1'

                                inlierImageData2 = []
                                # for 9 of the inlier datapoints in the LC sample
                                for i2, s in inliers.iterrows():
                                    imageData = archiveLCData[a].getImages(i2, [config.imgSize[imgSizeIndx]])
                                    inlierImageData2.append([imageData, f"{s['filefracday']}"])

                                    HDU = imageData[0][0]
                                    # data = HDU[0].data

                                    title = f'{AO.objectName} Inliers, ' \
                                            f'{config.imgSize[imgSizeIndx]} arcmin ' \
                                            f'({f}-band)'

                                    print(f'{title}, limiting mag / zero point difference: '
                                          f'{archiveLCData[a].table.loc[i2, "limitzpdiff"]} mag')
                                    bkg = estimateBackground(HDU, title, out=True)
                                    HDUbkg = copy.deepcopy(HDU)
                                    HDUbkg[0].data = bkg.background
                                    HDUdif = copy.deepcopy(HDU)
                                    HDUdif[0].data = HDU[0].data - bkg.background
                                    showBkgImages(HDU, HDUbkg, HDUdif, f'{AO.objectName} Inlier {s["filefracday"]}, '
                                                                       f'{config.imgSize[imgSizeIndx]}arcmin '
                                                                       f'({f}-band)')

                                outlierImageData2 = []
                                # for each of the outlier datapoints in the LC sample
                                for i2, s in outliers.iterrows():
                                    imageData = archiveLCData[a].getImages(i2, [config.imgSize[imgSizeIndx]])
                                    outlierImageData2.append([imageData, f"{s['filefracday']}"])

                                    HDU = imageData[0][0]
                                    # data = HDU[0].data

                                    title = f'{AO.objectName} Outliers, ' \
                                            f'{config.imgSize[imgSizeIndx]} arcmin ' \
                                            f'({f}-band)'

                                    print(f'{title}, limiting mag / zero point difference: '
                                          f'{archiveLCData[a].table.loc[i2, "limitzpdiff"]} mag')
                                    bkg = estimateBackground(HDU, title, out=True)
                                    HDUbkg = copy.deepcopy(HDU)
                                    HDUbkg[0].data = bkg.background
                                    HDUdif = copy.deepcopy(HDU)
                                    HDUdif[0].data = HDU[0].data - bkg.background
                                    showBkgImages(HDU, HDUbkg, HDUdif, f'{AO.objectName} Outlier {s["filefracday"]}, '
                                                                       f'{config.imgSize[imgSizeIndx]}arcmin '
                                                                       f'({f}-band)')

                                outlierImageData3 = []
                                # for each of the outlier datapoints in the LC sample
                                for i3, s in outliers.iterrows():
                                    imageData = archiveLCData[a].getImages(i3,
                                                                           [config.imgSize[imgSizeIndx]],
                                                                           data=True)
                                    outlierImageData3.append([imageData, f"{s['filefracday']}"])

                                    HDU = imageData[0][0]
                                    # data = HDU[0].data

            else:
                print(f'No data available or retrieved from {archives[a].name}')
                print()
        filtersToPlot = AO.get_filtersInAO()
        if AO.plotRowsTrue():
            fig, ax = AO.preparePlot(filtersToPlot)
            for a in archiveList:
                archiveLCData[a].plot(fig, ax, archives[a], filtersToPlot)
            AO.finalisePlot(fig, ax, filtersToPlot, archiveList)

    # datetime object containing current date and time
    # now = datetime.now()
    #
    # dataFileN = 'LCStats' + now.strftime('%Y%m%d_%H%M%S') + '.csv'
    #
    # with open(f'data/{dataFileN}', 'w') as f:
    #     objectsList.write(f, format='csv')
