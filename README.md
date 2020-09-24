# SustainabilityHotspots
Code repo for the manuscript: "Sustainability hotspots of changing global freshwater availability"


**In prep folder:** <br>
*(these scripts should be run first)* <br>
`Dim_data_res_harmonize.r` harmonizes all data input to 0.5 degrees <br> 
`GPCC_annualsums.r` calculates annual precipitation over 1985-2014 <br>
`GRACE_csv_convert.r` takes Rodell et al. (2018) source data and converts to raster in WGS84 <br>
`GRACE_precip_scaler.r` divides GRACE TWS trends by annual precipitation and scales to range of [-1, +1] for use as indicator <br>
`GriddedWithdrawals.r` de-compresses and sums water withdrawals for the year 2010 from Huang et al. (2018) <br>
`Mask_derive.r` creates mask that excludes Antarctica, oceans, and regions of earthquake interference <br>
`WGS84gridarea.r` calculates the area per 0.5 degree WGS grid cell based on the reference ellipsoid <br>

**Core scripts:** <br>
`gen_funs.r` creates custom geospatial functions to be used in `Hotspots.r` and `Dimensions.r` <br>
`Dimensions.r` derives the population, agricultural, economic, and environmental presence indicators, and multiplies each by the trend severity indicator to derive each dimension's hotspot indicator. <br>
`Hotspots.r` derives the inidivudal and overall sustainability hotspots based on analytical hierarchy process-derived weightings.
