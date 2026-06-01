# 3D-event-calculation-code
This file includes functions for 3D spatiotemporal identification and tracking of compound drought–heat events (CDHEs), as well as the associated vegetation loss assessment.
Note:
System requirements: The code requires R packages including ncdf4, dbscan, ggplot2, factoextra, zoo, raster, sf, and abind. All packages can be downloaded directly from CRAN, which is an open and free source. The code has been tested on Windows and Linux platforms.

Installation guide: The code can be directly run in the R language across Windows and Linux platforms. The R version is 4.3.0 or higher, and the associated packages can be installed in this R version. There is no need to spend time installing this code, just configure the required packages as listed above.

Usage: To run this code, you need gridded climate and vegetation data, as well as global basin shapefiles data. IDI integrates meteorological, hydrological, and agricultural drought information (precipitation, runoff, and soil moisture data). STI is calculated from near-surface temperature data. CDHI couples IDI with STI to identify CDHEs. SVI is derived from NDVI for vegetation loss assessment.
