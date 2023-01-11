var countries = ee.FeatureCollection('USDOS/LSIB/2017'); // Countries
var ppdensity = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density").select('population_density'); // Population density
var chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY'); // Precipitation
var cwtdef = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filter(ee.Filter.date('1981-01-01','2020-12-01')).select('def'); // Climate water deficit
var tmax   = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filter(ee.Filter.date('1981-01-01','2020-12-01')).select('tmmx'); // Maximum temperature
var aet    = ee.ImageCollection("IDAHO_EPSCOR/TERRACLIMATE").filter(ee.Filter.date('1981-01-01','2020-12-01')).select('aet'); // Actual evapotranspiration


// Countries list
var cList = ['Kenya','Senegal','Uganda','Nigeria','Sudan','Mali','Zimbabwe', "Ethiopia", "Guatemala", "Philippines"]; // Countries list
// One specific country
var country = countries.filter(ee.Filter.inList('COUNTRY_NA',['Kenya']));

//population density avg
var avgPopD = ppdensity.mean();
var avgPopDClipped = avgPopD.clip(country);

//precipitation p90

var yearlyRainfall = function(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');
  var filtered  = chirps.filter(ee.Filter.date(startDate, endDate));
  var total     = filtered.reduce(ee.Reducer.sum()).set({year: year, 'system:time_start':startDate});
  return total;
};

var years = ee.List.sequence(1981, 2020);
var yRain = ee.ImageCollection(years.map(yearlyRainfall));
var p90Rain = yRain.reduce(ee.Reducer.percentile([90]));
var p90RainClipped = p90Rain.clip(country);



var yearlyTmax = function(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');
  var filtered  = tmax.filter(ee.Filter.date(startDate, endDate));
  var median    = filtered.reduce(ee.Reducer.median()).set({year: year, 'system:time_start':startDate});
  return median;
};

var years = ee.List.sequence(1981, 2020);
var yTmax = ee.ImageCollection(years.map(yearlyTmax));
var p90Tmax = yTmax.reduce(ee.Reducer.percentile([90]));
var p90TmaxClipped = p90Tmax.clip(country);



var yearlyAET = function(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');
  var filtered  = aet.filter(ee.Filter.date(startDate, endDate));
  var median    = filtered.reduce(ee.Reducer.mean()).set({year: year, 'system:time_start':startDate});
  return median;
};

var years = ee.List.sequence(1981, 2020);
var yAET = ee.ImageCollection(years.map(yearlyAET));
var avgAET = yAET.mean();
var avgAETClipped = avgAET.clip(country);


var yearlyCWDf = function(year) {
  var startDate = ee.Date.fromYMD(year, 1, 1);
  var endDate   = startDate.advance(1, 'year');
  var filtered  = cwtdef.filter(ee.Filter.date(startDate, endDate));
  var median    = filtered.reduce(ee.Reducer.mean()).set({year: year, 'system:time_start':startDate});
  return median;
};

var years = ee.List.sequence(1981, 2020);
var yCWDf = ee.ImageCollection(years.map(yearlyCWDf));
var avgCWDf = yCWDf.mean();
var avgCWDfClipped = avgCWDf.clip(country);





// Access variables
// Image to GDrive: accessibility
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