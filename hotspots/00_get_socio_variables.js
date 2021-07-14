// African Crisis Observatory: Obtain land cover-population density indicators
// Alliance Bioversity-CIAT
// H. Achicanoy, 2021

// Load datasets
var countries = ee.FeatureCollection('USDOS/LSIB/2017'); // Countries
var srtm = ee.Image('CGIAR/SRTM90_V4').select('elevation'); // Elevation model 90 m
var landcover = ee.ImageCollection('MODIS/006/MCD12Q1').select('LC_Type1'); // Land cover 500 m
var ppdensity = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density").select('population_density'); // Population density
var hlth_accs = ee.Image('Oxford/MAP/accessibility_to_healthcare_2019'); // Accessibility to healthcare services 2019
var frct_srfc = ee.Image('Oxford/MAP/friction_surface_2019'); // Friction surface 2019
var lght_nght = ee.ImageCollection('NOAA/DMSP-OLS/NIGHTTIME_LIGHTS'); // Lightnigths

// Countries list
var cList = ['Kenya','Senegal','Uganda','Nigeria','Sudan','Mali','Zimbabwe']; // Countries list
// One specific country
var country = countries.filter(ee.Filter.inList('COUNTRY_NA',['Sudan']));
// Map.addLayer(country, {}, 'Country shp');

// Accessibility to healthcare 2019
var access = hlth_accs.select('accessibility');
var accessClipped = access.clip(country);

// Friction surface 2019
var frct = frct_srfc.select('friction');
var frctClipped = frct.clip(country);

//////////////////////////////////////////////
// Median
//////////////////////////////////////////////

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
var slght = lght_nght.reduce(ee.Reducer.stdDev());
var vlght = slght.divide(alght);
var vlghtClipped = vlght.clip(country);

//////////////////////////////////////////////
// Trend
//////////////////////////////////////////////

// Adds a band containing image date as years since 1981.
function createTimeBand(img) {
  var year = ee.Date(img.get('system:time_start')).get('year').subtract(1981);
  return ee.Image(year).byte().addBands(img);
}

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

//////////////////////////////////////////////
// Download processed files
//////////////////////////////////////////////

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

// CV: vLCvrClipped, vPopDClipped, vlghtClipped
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

// Trend: tPopDClipped, tlghtClipped
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
