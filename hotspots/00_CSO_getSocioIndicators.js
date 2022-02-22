// Climate Security Observatory: Obtain land cover-population density indicators
// Alliance Bioversity-CIAT
// H. Achicanoy, 2021
// https://code.earthengine.google.com/d53224432e22aa7fb5c9fd4cf915e0aa

// Load datasets
var countries = ee.FeatureCollection('USDOS/LSIB/2017'); // Countries
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
var wcover    = ee.ImageCollection("ESA/WorldCover/v100").first().select('Map'); // ESA WorldCover 10m v100
var aux       = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density").first().select('population_density');

// Countries list
var cList = ['Kenya','Senegal','Uganda','Nigeria','Sudan','Mali','Zimbabwe']; // Countries list
// One specific country
var country = countries.filter(ee.Filter.inList('COUNTRY_NA',['Kenya']));
// Map.addLayer(country, {}, 'Country shp');

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

// WorldCover map
var prj = aux.projection();

var wcoverMode = wcover
    // Force the next reprojection to aggregate instead of resampling.
    .reduceResolution({
      reducer: ee.Reducer.mode(),
      maxPixels: 10000
    })
    // Request the data at the scale and projection of the MODIS image.
    .reproject({
      crs: prj
    });
var wcoverModeClipped = wcoverMode.clip(country);

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

/////////////////////////
// Net primary product
/////////////////////////

var nppTS = npp.map(createTimeBand);
var tnpp = nppTS.reduce(ee.Reducer.sensSlope());
var tnppClipped = tnpp.clip(country).select('slope');

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
