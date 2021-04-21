## LCExtract
### Quickstart
To run application from python console 
>`from LCExtract.LCExtract import LCExtract`
> 
>`LCExtract()`
### Description
Currently a standalone application written in **Python** to take astronomical object positions and search for photometry 
data over time which produces a _lightcurve_ for the object.

Queries data from [Panoramic Survey Telescope And Rapid Response System](https://panstarrs.ifa.hawaii.edu/pswww/) 
(Pan-STARRS) DR2 as hosted at the [Mikulski Archive for Space Telescopes](https://archive.stsci.edu) (MAST),
[Palomar Transient Factory](https://www.ptf.caltech.edu) (PTF) and the 
[Zwicky Transient Facility](https://www.ztf.caltech.edu) (ZTF), the latter of which are hosted at
[NASA/IPAC Infrared Science Archive](https://irsa.ipac.caltech.edu/frontpage/) (IRSA). Other time domain based facilities will be added in 
to the query over time to provide a more comprehensive history of the object's luminosity variation.

Options selectable by the user are
* File or Manual entry
  * File input requires `name`, `RA`, `DEC`, `description` CSV separation - see example file [test_objects.csv](https://github.com/Pommers/LCExtract/blob/master/data/test_objects.csv)
  * Manual entry allows named or positional coordinate entry
* Archive / Catalog 
  * Pan-STARRS
  * PTF 
  * ZTF
  * All - consolidate data from facilities
* Photometric filter selection
  * Filters `g`, `r`, `i`, `z`, `y`, `R` are default filters, but note
    * a subset is user-selectable
    * not all filters available within archive(s), e.g., ZTF uses only `g`, `r`, `i` filters
    * not all filters available for an object

Application will calculate a number of summary statistics from the data returned for an object and provide a plot of any archive data returned.
### Data files
Data files are expected to be located in the working directory, preferably under a subdirectory `data/`