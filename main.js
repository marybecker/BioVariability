var options = {
    center: [41.6, -72.6],
    zoom: 10,
    zoomControl: false
};

var map = L.map('map', options);

// add zoom control to top right window position
L.control.zoom({
    position: 'topleft'
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
var tiles = L.tileLayer(basemap_url, basemap_attributes).addTo(map);

//define common styles for all circle markers
var commonStyles = {
    weight: 1,
    stroke: 1,
    fillOpacity: .8
};

var breaks = [-452,-60,-20,0,20,60,260];

function getRadius (area){
    var r = get_bin([area],bins=breaks)[0];
    var R = {0:12,1:10,3:8,4:6,5:4,6:2};
    return R[r];
}


 function get_bin(data,bins=[-134,-37,-13,0,16,51,263]){
    var B = [];
    for(var i = 0; i<data.length; i++){
        for(var b=1; b<bins.length; b++){
            if(data[i]>=bins[b-1] && data[i]<bins[b]){
                B.push(b);
            }
         }
     }
    return B;
 }

var color = [   [-1,'Cold to Not Cold','#018571'],
                [0,'Stable','#f5f5f5'],
                [1,'Not Cold to Cold','#dfc27d']];

function getColor (d){
        if  (d == color[0][0]){
            return color[0][2]
        }if  (d == color[1][0]){
            return color[1][2]
        }if (d == color[2][0]){
            return color[2][2]
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

function drawLegend (color){
    var legend = L.control({position: 'topright'});
    legend.onAdd = function(){
        var div = L.DomUtil.create('div','legend');
        div.innerHTML = "<h3>" +'CW Fish Diff'+ "</h3>";
        for (var i = 0; i < color.length; i++) {
            var c = getColor(color[i][0], color);
            div.innerHTML +=
                '<span style="background:' + c + '"></span> ' +
                '<label>'+(color[i][1]).toLocaleString()+ '</label>';
        }

        return div;
    };
    legend.addTo(map);
}

drawLegend(color);

