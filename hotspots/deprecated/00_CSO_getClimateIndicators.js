// Climate Security Observatory: Obtain climate derived indicators
// Alliance Bioversity-CIAT
// H. Achicanoy, 2021
// https://code.earthengine.google.com/1e7c8e56741f3fa1774fbae535ce6a62

// Load datasets
var countries = ee.FeatureCollection('USDOS/LSIB/2017'); // Countries
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY'); // Precipitation
var cwtdef = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filter(ee.Filter.date('1981-01-01','2020-12-01')).select('def'); // Climate water deficit
var tmax   = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filter(ee.Filter.date('1981-01-01','2020-12-01')).select('tmmx'); // Maximum temperature
var aet    = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filter(ee.Filter.date('1981-01-01','2020-12-01')).select('aet'); // Actual evapotranspiration

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
