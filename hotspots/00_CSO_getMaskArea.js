// Climate Security Observatory: Obtain mask area
// Alliance Bioversity-CIAT
// H. Achicanoy, 2021
// https://code.earthengine.google.com/03fa0d82a56d411f103e1305c0cd15c6

// Load datasets
var countries = ee.FeatureCollection('USDOS/LSIB/2017'); // Countries
var wcover    = ee.ImageCollection("ESA/WorldCover/v100").first().select('Map'); // ESA WorldCover 10m v100
var aux       = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density").first().select('population_density');

// Countries list
var cList = ['Kenya','Senegal','Uganda','Nigeria','Sudan','Mali','Zimbabwe','Guatemala', 'Philippines']; // Countries list
// One specific country
var country = countries.filter(ee.Filter.inList('COUNTRY_NA',['Philippines']));
// Map.addLayer(country, {}, 'Country shp');

Map.centerObject(country, 6)
// WorldCover map
var prj = aux.projection();

//var wcov2 =wcover.eq(20).and(wcover.eq())
// var wcov = wcover.eq(20).add(wcover.eq(30).add(wcover.eq(40).add(wcover.eq(50))))
var wcov = wcover.eq(30).add(wcover.eq(40).add(wcover.eq(50)));

var wcoverClipped = wcov.clip(country);
print("wcov",wcoverClipped)
Map.addLayer(wcoverClipped.updateMask(wcoverClipped.eq(1)),{},"wcov_print")
Map.centerObject(country,5)

// var palette = ['red', 'white', 'blue'];
// Map.addLayer(wcoverClipped, {palette: palette, min: 0, max: 1, bands: ['Map']}, 'Mask area');

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
