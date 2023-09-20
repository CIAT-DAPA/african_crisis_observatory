// Climate Security Observatory: Obtain climate derived indicators
// Alliance Bioversity-CIAT
// H. Achicanoy, 2021
// https://code.earthengine.google.com/1e7c8e56741f3fa1774fbae535ce6a62

// Load datasets

// Load datasets
var countries = ee.FeatureCollection('USDOS/LSIB/2017'); // Countries
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY'); // Precipitation
var wcover    = ee.ImageCollection("ESA/WorldCover/v100").first().select('Map'); // ESA WorldCover 10m v100
var aux       = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density").first().select('population_density');

var cwtdef = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filter(ee.Filter.date('1981-01-01','2020-12-01')).select('def'); // Climate water deficit
var tmax   = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filter(ee.Filter.date('1981-01-01','2020-12-01')).select('tmmx'); // Maximum temperature
var aet    = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filter(ee.Filter.date('1981-01-01','2020-12-01')).select('aet'); // Actual evapotranspiration
var srtm = ee.Image('CGIAR/SRTM90_V4').select('elevation'); // Elevation model 90 m
var landcover = ee.ImageCollection('MODIS/006/MCD12Q1').filterDate('2019-01-01').select('LC_Type1'); // Land cover 500 m // .last()
var ppdensity = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density").select('population_density'); // Population density
var hlth_accs = ee.Image('Oxford/MAP/accessibility_to_healthcare_2019'); // Accessibility to healthcare services 2019
var frct_srfc = ee.Image('Oxford/MAP/friction_surface_2019'); // Friction surface 2019
var lght_nght = ee.ImageCollection('NOAA/DMSP-OLS/NIGHTTIME_LIGHTS'); // Lightnigths
var npp = ee.ImageCollection("MODIS/006/MOD17A3HGF").select('Npp'); // Net Primary Production 500 m
var carbon_ct = ee.Image("OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02").select('b60'); // Soil organic carbon content 250 m
var clay_ct   = ee.Image("OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02").select('b60'); // Clay content 250 m
var ph_ct     = ee.Image("OpenLandMap/SOL/SOL_PH-H2O_USDA-4C1A2A_M/v02").select('b60'); // Soil pH in H2O 250 m 
var deforest  = ee.Image("UMD/hansen/global_forest_change_2020_v1_8").select('lossyear'); // Year of gross forest cover loss event. 2000-2020

// Countries list
var cList = ['Kenya','Senegal','Uganda','Nigeria','Sudan','Mali','Zimbabwe', 'Philippines']; // Countries list
// One specific country
var country = countries.filter(ee.Filter.inList('COUNTRY_NA',['Philippines']));
// Map.addLayer(country, {}, 'Country shp');

Map.centerObject(country, 5)
//////////////////////////////////////////////
// Median
//////////////////////////////////////////////

/////////////////////////
// Precipitation
/////////////////////////

// Function to calculate rainfall for 1 year
var yearlyRainfall = function(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');
  var filtered  = chirps.filter(ee.Filter.date(startDate, endDate));
  var total     = filtered.reduce(ee.Reducer.sum()).set({year: year, 'system:time_start':startDate});
  return total;
};

var years = ee.List.sequence(1981, 2020);
var yRain = ee.ImageCollection(years.map(yearlyRainfall));
var mRain = yRain.median();
var mRainClipped = mRain.clip(country);
var p90Rain = yRain.reduce(ee.Reducer.percentile([90]));
var p90RainClipped = p90Rain.clip(country);

// var palette = ['red', 'white', 'blue'];
// Map.addLayer(mRainClipped, {palette: palette, min: 0, max: 1000, bands: ['precipitation_sum']}, 'Median');

/////////////////////////
// Climate water deficit
/////////////////////////

// Function to calculate climate water deficit for 1 year
var yearlyCWDf = function(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');
  var filtered  = cwtdef.filter(ee.Filter.date(startDate, endDate));
  var median    = filtered.reduce(ee.Reducer.median()).set({year: year, 'system:time_start':startDate});
  return median;
};

var years = ee.List.sequence(1981, 2020);
var yCWDf = ee.ImageCollection(years.map(yearlyCWDf));
var mCWDf = yCWDf.median();
var mCWDfClipped = mCWDf.clip(country);
var avgCWDf = yCWDf.mean();
var avgCWDfClipped = avgCWDf.clip(country);

// var palette = ['red', 'white', 'blue'];
// Map.addLayer(mCWDfClipped, {palette: palette, min: 0, max: 2000, bands: ['def_median']}, 'Def median');

/////////////////////////
// Maximum temperature
/////////////////////////

// Function to calculate maximum temperature for 1 year
var yearlyTmax = function(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');
  var filtered  = tmax.filter(ee.Filter.date(startDate, endDate));
  var median    = filtered.reduce(ee.Reducer.median()).set({year: year, 'system:time_start':startDate});
  return median;
};

var years = ee.List.sequence(1981, 2020);
var yTmax = ee.ImageCollection(years.map(yearlyTmax));
var mTmax = yTmax.median();
var mTmaxClipped = mTmax.clip(country);
var p90Tmax = yTmax.reduce(ee.Reducer.percentile([90]));
var p90TmaxClipped = p90Tmax.clip(country);


// var palette = ['red', 'white', 'blue'];
// Map.addLayer(mTmaxClipped, {palette: palette, min: 200, max: 400, bands: ['tmmx_median']}, 'Tmax median');

/////////////////////////
// Actual evapotranspiration
/////////////////////////

// Function to calculate actual evapotranspiration for 1 year
var yearlyAET = function(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');
  var filtered  = aet.filter(ee.Filter.date(startDate, endDate));
  var median    = filtered.reduce(ee.Reducer.median()).set({year: year, 'system:time_start':startDate});
  return median;
};

var years = ee.List.sequence(1981, 2020);
var yAET = ee.ImageCollection(years.map(yearlyAET));
var mAET = yAET.median();
var mAETClipped = mAET.clip(country);
var avgAET = yAET.mean();
var avgAETClipped = avgAET.clip(country);

// var palette = ['red', 'white', 'blue'];
// Map.addLayer(mCWDfClipped, {palette: palette, min: 0, max: 2000, bands: ['def_median']}, 'Def median');

//////////////////////////////////////////////
// Coefficient of variation
//////////////////////////////////////////////

/////////////////////////
// Precipitation
/////////////////////////

var aRain = yRain.mean();
var sRain = yRain.reduce(ee.Reducer.stdDev());
var vRain = sRain.divide(aRain);
var vRainClipped = vRain.clip(country);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(vRainClipped, {palette: palette, min: -1, max: 1, bands: ['precipitation_sum_stdDev']}, 'CV');

/////////////////////////
// Climate water deficit
/////////////////////////

var aCWDf = yCWDf.mean();
var sCWDf = yCWDf.reduce(ee.Reducer.stdDev());
var vCWDf = sCWDf.divide(aCWDf);
var vCWDfClipped = vCWDf.clip(country);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(vCWDfClipped, {palette: palette, min: -1, max: 1, bands: ['def_median_stdDev']}, 'CV');

/////////////////////////
// Maximum temperature
/////////////////////////

var aTmax = yTmax.mean();
var sTmax = yTmax.reduce(ee.Reducer.stdDev());
var vTmax = sTmax.divide(aTmax);
var vTmaxClipped = vTmax.clip(country);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(vTmaxClipped, {palette: palette, min: -1, max: 0.5, bands: ['tmmx_median_stdDev']}, 'CV');

/////////////////////////
// Actual evapotranspiration
/////////////////////////

var aAET = yAET.mean();
var sAET = yAET.reduce(ee.Reducer.stdDev());
var vAET = sAET.divide(aAET);
var vAETClipped = vAET.clip(country);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(vCWDfClipped, {palette: palette, min: -1, max: 1, bands: ['def_median_stdDev']}, 'CV');

//////////////////////////////////////////////
// Trend
//////////////////////////////////////////////

// Adds a band containing image date as years since 1981.
function createTimeBand(img) {
  var year = ee.Date(img.get('system:time_start')).get('year').subtract(1981);
  return ee.Image(year).byte().addBands(img);
}

/////////////////////////
// Precipitation
/////////////////////////

var yRainTS = ee.ImageCollection(years.map(yearlyRainfall)).map(createTimeBand);
var tRain = yRainTS.reduce(ee.Reducer.sensSlope());
var tRainClipped = tRain.clip(country).select('slope');
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(tRainClipped, {palette: palette, min: -30, max: 30, bands: ['slope']}, 'Slope');

/////////////////////////
// Climate water deficit
/////////////////////////

var yCWDfTS = ee.ImageCollection(years.map(yearlyCWDf)).map(createTimeBand);
var tCWDf = yCWDfTS.reduce(ee.Reducer.sensSlope());
var tCWDfClipped = tCWDf.clip(country).select('slope');
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(tCWDfClipped, {palette: palette, min: -30, max: 30, bands: ['slope']}, 'Slope');

/////////////////////////
// Maximum temperature
/////////////////////////

var yTmaxTS = ee.ImageCollection(years.map(yearlyTmax)).map(createTimeBand);
var tTmax = yTmaxTS.reduce(ee.Reducer.sensSlope());
var tTmaxClipped = tTmax.clip(country).select('slope');
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(tTmaxClipped, {palette: palette, min: -2, max: 2, bands: ['slope']}, 'Slope');

/////////////////////////
// Actual evapotranspiration
/////////////////////////

var yAETTS = ee.ImageCollection(years.map(yearlyAET)).map(createTimeBand);
var tAET = yAETTS.reduce(ee.Reducer.sensSlope());
var tAETClipped = tAET.clip(country).select('slope');
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(tCWDfClipped, {palette: palette, min: -30, max: 30, bands: ['slope']}, 'Slope');

//////////////////////////////
// WOLRDCOVER ///////////////
////////////////////////////

var prj = aux.projection();

//var wcov2 =wcover.eq(20).and(wcover.eq())
// var wcov = wcover.eq(20).add(wcover.eq(30).add(wcover.eq(40).add(wcover.eq(50))))
var wcov = wcover.eq(30).add(wcover.eq(40).add(wcover.eq(50)));
var wcoverClipped = wcov.clip(country);
var wcoverMode = wcov
    // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.max(),
      maxPixels: 10000
    })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: prj
    });
var wcoverModeClipped = wcoverMode.clip(country);



///////////////////////////////////////////
/////// SOCIECONOMIC VARIABLE ////////////
/////////////////////////////////////////

// Accessibility to healthcare 2019
var access = hlth_accs.select('accessibility');
var accessClipped = access.clip(country);

// Friction surface 2019
var frct = frct_srfc.select('friction');
var frctClipped = frct.clip(country);

// Upper limit NPP
var nppUL = npp.max();
var nppULClipped = nppUL.clip(country);

// Soil organic carbon content
var carbonClipped = carbon_ct.clip(country);

// Clay content
var clayClipped = clay_ct.clip(country);

// Soil pH in H2O
var pHClipped = ph_ct.clip(country);

// Forest cover loss
var deforestClipped = deforest.clip(country);

/////////////////////////
// Land cover
/////////////////////////

var mLCvr = landcover.median();
var mLCvrClipped = mLCvr.clip(country);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(mLCvrClipped, {palette: palette, min: 1, max: 17, bands: ['LC_Type1']}, 'Land Cover');

/////////////////////////
// Population density
/////////////////////////

var mPopD = ppdensity.median();
var mPopDClipped = mPopD.clip(country);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(mPopDClipped, {palette: palette, min: 0, max: 1000, bands: ['population_density']}, 'Population density');

/////////////////////////
// Lightnights
/////////////////////////

var mlght = lght_nght.select('stable_lights').median();
var mlghtClipped = mlght.clip(country);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(mlghtClipped, {palette: palette, min: 0, max: 63, bands: ['stable_lights']}, 'Nightlights');

//////////////////////////////////////////////
// Coefficient of variation
//////////////////////////////////////////////

/////////////////////////
// Land cover
/////////////////////////

var aLCvr = landcover.mean();
var sLCvr = landcover.reduce(ee.Reducer.stdDev());
var vLCvr = sLCvr.divide(aLCvr);
var vLCvrClipped = vLCvr.clip(country);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(vLCvrClipped, {palette: palette, min: 0, max: 1, bands: ['LC_Type1_stdDev']}, 'Land cover CV');

/////////////////////////
// Population density
/////////////////////////

var aPopD = ppdensity.mean();
var sPopD = ppdensity.reduce(ee.Reducer.stdDev());
var vPopD = sPopD.divide(aPopD);
var vPopDClipped = vPopD.clip(country);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(vPopDClipped, {palette: palette, min: 0, max: 1, bands: ['population_density_stdDev']}, 'Pop. Density CV');

/////////////////////////
// Lightnights
/////////////////////////

var alght = lght_nght.select('stable_lights').mean();
var slght = lght_nght.select('stable_lights').reduce(ee.Reducer.stdDev());
var vlght = slght.divide(alght);
var vlghtClipped = vlght.clip(country);

/////////////////////////
// Net primary product
/////////////////////////

var anpp = npp.mean();
var snpp = npp.reduce(ee.Reducer.stdDev());
var vnpp = snpp.divide(anpp);
var vnppClipped = vnpp.clip(country);

//////////////////////////////////////////////
// Trend
//////////////////////////////////////////////

/////////////////////////
// Population density
/////////////////////////

var yppdensityTS = ppdensity.map(createTimeBand);
var tPopD = yppdensityTS.reduce(ee.Reducer.sensSlope());
var tPopDClipped = tPopD.clip(country).select('slope');
// print(tPopDClipped);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(tRainClipped, {palette: palette, min: -30, max: 30, bands: ['slope']}, 'Slope');

/////////////////////////
// Lightnights
/////////////////////////

var lght_nghtTS = lght_nght.select('stable_lights').map(createTimeBand);
var tlght = lght_nghtTS.reduce(ee.Reducer.sensSlope());
var tlghtClipped = tlght.clip(country).select('slope');
// print(tlghtClipped);
// var palette = ['red', 'white', 'blue'];
// Map.addLayer(tlghtClipped, {palette: palette, min: -5, max: 5, bands: ['slope']}, 'Slope');

/////////////////////////
// Net primary product
/////////////////////////

var nppTS = npp.map(createTimeBand);
var tnpp = nppTS.reduce(ee.Reducer.sensSlope());
var tnppClipped = tnpp.clip(country).select('slope');




//////////////////////////////////////////////
// Download processed files
//////////////////////////////////////////////

// Median: mRainClipped, mCWDfClipped, mTmaxClipped, mAETClipped
// Image to GDrive
Export.image.toDrive({
  image: mRainClipped.float(),
  description: 'Median_rainfall',
  folder: 'GEE data',
  fileNamePrefix: 'medn_prec',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: mCWDfClipped.float(),
  description: 'Median_CWDf',
  folder: 'GEE data',
  fileNamePrefix: 'medn_cwdf',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: mTmaxClipped.float(),
  description: 'Median_Tmax',
  folder: 'GEE data',
  fileNamePrefix: 'medn_tmax',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: mAETClipped.float(),
  description: 'Median_AET',
  folder: 'GEE data',
  fileNamePrefix: 'medn_aet',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// CV: vRainClipped, vCWDfClipped, vTmaxClipped, vAETClipped
// Image to GDrive
Export.image.toDrive({
  image: vRainClipped.float(),
  description: 'CV_rainfall',
  folder: 'GEE data',
  fileNamePrefix: 'cvar_prec',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: vCWDfClipped.float(),
  description: 'CV_CWDf',
  folder: 'GEE data',
  fileNamePrefix: 'cvar_cwdf',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: vTmaxClipped.float(),
  description: 'CV_Tmax',
  folder: 'GEE data',
  fileNamePrefix: 'cvar_tmax',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: vAETClipped.float(),
  description: 'CV_AET',
  folder: 'GEE data',
  fileNamePrefix: 'cvar_aet',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// Trend: tRainClipped, tCWDfClipped, tTmaxClipped, tAETClipped
// Image to GDrive
Export.image.toDrive({
  image: tRainClipped.float(),
  description: 'Trend_rainfall',
  folder: 'GEE data',
  fileNamePrefix: 'trnd_prec',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: tCWDfClipped.float(),
  description: 'Trend_CWDf',
  folder: 'GEE data',
  fileNamePrefix: 'trnd_cwdf',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: tTmaxClipped.float(),
  description: 'Trend_Tmax',
  folder: 'GEE data',
  fileNamePrefix: 'trnd_tmax',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: tAETClipped.float(),
  description: 'Trend_AET',
  folder: 'GEE data',
  fileNamePrefix: 'trnd_aet',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// NEW VARS 

Export.image.toDrive({
  image: avgPopDClipped.float(),
  description: 'average_Pop_density',
  folder: 'GEE data',
  fileNamePrefix: 'avg_popd',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: p90RainClipped.float(),
  description: 'p90_precipiation',
  folder: 'GEE data',
  fileNamePrefix: 'p90_prec',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: p90TmaxClipped.float(),
  description: 'p90_tempmax',
  folder: 'GEE data',
  fileNamePrefix: 'p90_tmax',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: avgAETClipped.float(),
  description: 'avg_AETClipped',
  folder: 'GEE data',
  fileNamePrefix: 'avg_aet',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: avgCWDfClipped.float(),
  description: 'avg_CWDfClipped',
  folder: 'GEE data',
  fileNamePrefix: 'avg_cwdf',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});


/// LAND cover

Export.image.toDrive({
  image: wcoverClipped.float(),
  description: 'WorldCover',
  folder: 'GEE data',
  fileNamePrefix: 'wcover',
  region: country,//.geometry().bounds().getInfo(),
  scale: 10,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

/// SOCIECONOMIC- VARIABLES


// Access variables
// Image to GDrive: accessibility
Export.image.toDrive({
  image: accessClipped.float(),
  description: 'Accessibility2healthcare',
  folder: 'GEE data',
  fileNamePrefix: 'acess',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
// Image to GDrive: accessibility
Export.image.toDrive({
  image: frctClipped.float(),
  description: 'FrictionSurface',
  folder: 'GEE data',
  fileNamePrefix: 'friction',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
// Image to GDrive: Upper bound NPP
Export.image.toDrive({
  image: nppULClipped.float(),
  description: 'UpperBoundNPP',
  folder: 'GEE data',
  fileNamePrefix: 'ub_npp',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
// Image to GDrive: Soil organic carbon content
Export.image.toDrive({
  image: carbonClipped.float(),
  description: 'CarbonContent',
  folder: 'GEE data',
  fileNamePrefix: 'soil_carbon',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
// Image to GDrive: Soil clay content
Export.image.toDrive({
  image: clayClipped.float(),
  description: 'ClayContent',
  folder: 'GEE data',
  fileNamePrefix: 'soil_clay',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
// Image to GDrive: Soil clay content
Export.image.toDrive({
  image: pHClipped.float(),
  description: 'SoilpH',
  folder: 'GEE data',
  fileNamePrefix: 'soil_pH',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: deforestClipped.float(),
  description: 'Deforestation',
  folder: 'GEE data',
  fileNamePrefix: 'deforest',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: wcoverModeClipped.float(),
  description: 'WorldCover',
  folder: 'GEE data',
  fileNamePrefix: 'wcover',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
})

// Median: mLCvrClipped, mPopDClipped, mlghtClipped
// Image to GDrive
Export.image.toDrive({
  image: mLCvrClipped.float(),
  description: 'Median_Land_cover',
  folder: 'GEE data',
  fileNamePrefix: 'medn_lcvr',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: mPopDClipped.float(),
  description: 'Median_Pop_density',
  folder: 'GEE data',
  fileNamePrefix: 'medn_popd',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: mlghtClipped.float(),
  description: 'Median_Nightlights',
  folder: 'GEE data',
  fileNamePrefix: 'medn_lght',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// CV: vLCvrClipped, vPopDClipped, vlghtClipped, vnppClipped
Export.image.toDrive({
  image: vLCvrClipped.float(),
  description: 'CV_Land_cover',
  folder: 'GEE data',
  fileNamePrefix: 'cvar_lcvr',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: vPopDClipped.float(),
  description: 'CV_Pop_density',
  folder: 'GEE data',
  fileNamePrefix: 'cvar_popd',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: vlghtClipped.float(),
  description: 'CV_Nightlights',
  folder: 'GEE data',
  fileNamePrefix: 'cvar_lght',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: vnppClipped.float(),
  description: 'CV_NPP',
  folder: 'GEE data',
  fileNamePrefix: 'cvar_npp',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// Trend: tPopDClipped, tlghtClipped, tnppClipped
Export.image.toDrive({
  image: tPopDClipped.float(),
  description: 'Trend_Pop_density',
  folder: 'GEE data',
  fileNamePrefix: 'trnd_popd',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: tlghtClipped.float(),
  description: 'Trend_Nightlights',
  folder: 'GEE data',
  fileNamePrefix: 'trnd_lght',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
Export.image.toDrive({
  image: tnppClipped.float(),
  description: 'Trend_NPP',
  folder: 'GEE data',
  fileNamePrefix: 'trnd_npp',
  region: country.geometry().bounds().getInfo(),
  scale: 1000,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});


