var options = {
    center: [41.6032, -73.0877],
    zoom: 9,
    zoomControl: false
};

var map = L.map('map', options);

// add zoom control to top right window position
L.control.zoom({
    position: 'topright'
}).addTo(map);

// Get basemap URL from Leaflet Providers
var basemap_url = 'https://cartodb-basemaps-{s}.global.ssl.fastly.net/rastertiles/dark_all/{z}/{x}/{y}.png';

// Get basemap attributes from Leaflet Providers
var basemap_attributes = {
    attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a> &copy; <a href="http://cartodb.com/attributions">CartoDB</a>',
    subdomains: 'abcd',
    maxZoom: 19
};
// requests some map tiles
var tiles = L.tileLayer(basemap_url, basemap_attributes).addTo(map);;

//define common styles for all circle markers
var commonStyles = {
    weight: 1,
    stroke: 1,
    fillOpacity: .8
};

function getRadius (area){
    var radius = Math.sqrt(area/Math.PI);
    return radius * 1.7;
}

// function getSize (d){
//     if (d > )
// }

function getColor (d){
        if  (d == -1){
            return "#018571"
        }if  (d == 0){
            return "#f5f5f5"
        }if (d == 1){
            return "#dfc27d"
        }
}

//function to define style
function style(feature){
    return{
        color: getColor(feature.properties.CWYrDiff20_YrMaxColdDiff),
        fillColor: getColor(feature.properties.CWYrDiff20_YrMaxColdDiff),
        radius: getRadius(feature.properties.CWYrDiff20_YrMaxDiff)
    };
}

//mouse on mouse off
function highlightFeature(e) {
    var layer = e.target;
    layer.setStyle({
        weight: 5,
        color: '#666',
        dashArray: '',
        fillOpacity: 0.7
    });
}

function resetHighlight(e) {
    geojson.resetStyle(e.target);
}

var geojson =
    L.geoJson(sites, {
        pointToLayer: function (feature, latlng) {
            return L.circleMarker(latlng,commonStyles);
        },
        style: style,
        onEachFeature: function(feature,layer) {
                layer.on({
                    mouseover: highlightFeature,
                    mouseout: resetHighlight
                });
                layer.bindTooltip(
                    '<b>' + layer.feature.properties.STA_SEQ + '</b>' +
                    '<br>' + layer.feature.properties.Station_Name +'<br>' +
                    '<b>' + layer.feature.properties.CWYrDiff20_SampYr1+ ': </b>' +
                    layer.feature.properties.CWYrDiff20_CWFishYr1 + ' FishPer100M' + '<br>' +
                    '<b>' + layer.feature.properties.CWYrDiff20_SampYr2 + ': </b>' +
                    layer.feature.properties.CWYrDiff20_CWFishYr2 + ' FishPer100M'+'<br>'+
                    '<b> Diff: </b>' +
                    layer.feature.properties.CWYrDiff20_YrMaxDiff + ' FishPer100M')
        }
    }).addTo(map);


