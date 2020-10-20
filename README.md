Code repo for the manuscript: "Sustainability hotspots of changing global freshwater availability" (Huggins et al., in review).

**In Preprocessing folder:** <br>
*(Required to be executed before the core scripts)* <br>
`Dim_data_res_harmonize.r` harmonizes all data input to 0.5 degrees. <br> 
`GPCC_annualsums.r` calculates annual precipitation over 1985-2014 using GPCC data. <br>
`GRACE_csv_convert.r` takes Rodell et al. (2018) source data and converts to raster in WGS84. <br>
`GRACE_precip_scaler.r` divides GRACE TWS trends by annual precipitation and scales to range of [-1, +1] for use as indicator. <br>
`GriddedWithdrawals_EconomicSectors_2010.r` sums water withdrawals for electricity generation, manufaturing, livestock, and mining for the year 2010 from Huang et al. (2018). <br>
`Mask_derive.r` creates mask that excludes Antarctica, oceans, and regions of earthquake interference. <br>
<br>
**In main folder (core scripts):** <br>
`gen_funs.r` creates custom geospatial functions to be used in `Hotspots.r` and `Dimensions.r` <br>
These include: <br>
    `RasterAreaPercentiles` converts a distributed grid into area-weighted percentiles, with masking optional. <br>
    `RasterGridBuffer` buffers land mask by 1 grid cell to ensure no land areas are inadvertently masked. <br>
    `tmap_clipproj` prepares global rasters for plotting in tmap. <br>
    `hotspot_id_smoother` post-processes sustainability hotspots using two filters. <br>
    `WGS84_areaRaster` returns a raster with grid cell values representing WGS84 grid cell areas (up to 2.5 minute resolution). <br>
<br>
`Dimensions.r` derives the population, agricultural, economic, and environmental presence indicators, and multiplies each by the trend severity indicator to derive each dimension's hotspot indicator. <br>
`DimHistograms.r` plots histograms of population count, crop calories, GDP at PPP, and water sensitive environment surface area agaist the relative water availability treds. <br>
`Hotspots.r` derives the inidivudal and overall sustainability hotspots based on analytical hierarchy process-derived weightings. <br>
<br>
**In Shiny folder:** <br>
`app.r` base script used to control Shiny Web App currently hosted on https://waterhotspots.weebly.com/. <br>
