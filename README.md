[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.14635.svg)](http://dx.doi.org/10.5281/zenodo.14635)

A regional perspective on Holocene fire-climate-human interactions in the Pacific Northwest of North America
====================

January 27, 2014
*Coded by*: Simon J. Goring
*Report bugs to*: Simon J. Goring - goring@wisc.edu

This code-base accompanies the paper:

Walsh MK, Marlon JR, Goring S, Brown KJ, Gavin D. A regional perspective on Holocene fire-climate-human interactions in the Pacific Northwest of North America. *Annals of the Association of American Geographers.* Accepted: AN-2014-0079

The code can be used to generate and replicate the fPCA analysis as performed in the paper.  This code-base includes raw data for the projects.  All files in the `data` folder are predicated with their Global Charcoal Database (GCD) code.  Some datasets are not yet accessioned to the GCD, but will (or should be) accessioned by the next major GCD release.  Metadata for sites is included in the file `PNW master sitelist_Jan 30.csv`.

In addition, several datasets, particularly shapefiles, are required to undertake this analysis, in particular plotting.  These files are found at:

* [http://wwwn.cdc.gov/epiinfo/script/shapefiles.aspx]()
* [http://www.geobase.ca/geobase/en/search.do?produit=cgb&language=en]()
* [http://www.cec.org/atlas/files/Terrestrial_Ecoregions_L3/TerrestrialEcoregions_L3_Shapefile.zip]()

This analysis is carried out in the file `charcoal.to.fpca_supplement.R`.  Helper functions are found in `get.fpca.R`.  To run the analysis the following R libraries are required:

* [fda](http://CRAN.R-project.org/package=fda)
* [rgdal](http://CRAN.R-project.org/package=rgdal)
* [ggplot2](http://had.co.nz/ggplot2/book)
* [reshape2](http://www.jstatsoft.org/v21/i12/)
* [gridExtra](http://CRAN.R-project.org/package=gridExtra)
* [mgcv](http://cran.r-project.org/web/packages/mgcv/index.html)

References:
--------------
Auguie B. 2012. gridExtra: functions in Grid graphics. R package version 0.9.1. [http://CRAN.R-project.org/package=gridExtra]()
  
Bivand R, Keitt T & Rowlingson B (2014). rgdal: Bindings for the Geospatial Data Abstraction Library. R package version 0.9-1. [http://CRAN.R-project.org/package=rgdal]()

Ramsay JO, Wickham H, Graves S & Hooker G. 2014. fda: Functional Data Analysis. R package version 2.4.4. [http://CRAN.R-project.org/package=fda]()

Wickham H. 2007. Reshaping Data with the reshape Package. *Journal of Statistical Software*, 21: 1-20. URL http://www.jstatsoft.org/v21/i12/.
  
Wickham H. 2009. *ggplot2: elegant graphics for data analysis*. Springer New York.

Wood, S.N. 2011. Fast stable restricted maximum likelihood and marginal likelihood estimation of semiparametric generalized linear models. *Journal of the Royal Statistical Society (B)* 73:3-36.

Wood, S.N. 2004. Stable and efficient multiple smoothing parameter estimation for
generalized additive models. *Journal of the American Statistical Association.* 99:673-686.

Wood, S.N. 2006. *Generalized Additive Models: An Introduction with R.* Chapman and
Hall/CRC.

Wood, S.N. 2003. Thin-plate regression splines. *Journal of the Royal Statistical Society (B)* 65:95-114.

Wood, S.N. 2000. Modelling and smoothing parameter estimation with multiple quadratic
penalties. *Journal of the Royal Statistical Society (B)* 62:413-428.
