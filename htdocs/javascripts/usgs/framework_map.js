/**
 * Namespace: Framework_Map
 *
 * Framework_Map is a JavaScript library to provide a set of functions to build
 *  a Framework Web Site.
 *
 * version 2.14
 * April 10, 2024
*/

/*
###############################################################################
# Copyright (c) Oregon Water Science Center
# 
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
###############################################################################
*/
var grid_projection    = null;
var latlong_projection = null;
var grid_coordinates   = null;
var grid_graphic       = null;
var xy_multiplier      = 1.0;
var xy_units           = null;
var z_multiplier       = 1.0;
var z_units            = null;
var graph_x_axis       = null;
var graph_y_axis       = null;
      
var cell_geometry      = true;

var markerLayer;
var profileLayer;
var longlats           = [];

var myZoomFlag         = false;

var sliderValue        = 50;

// Retrieve framework information
//
var frameworkFile      = "framework_parameters.js";
        
var myInfo             = getInfo(frameworkFile);
        
var aboutFiles         = myInfo.aboutFiles;
console.log(myInfo);
console.log(aboutFiles);
        
var color_file         = myInfo.color_file;
        
var removecButtonUrl   = 'css/usgs/images/removeLocation.png';
var xsecButtonUrl      = 'css/usgs/images/xsec.png';

// Prepare when the DOM is ready 
//
$(document).ready(function() 
  {
  // Loading message
  //
  message = "Processing framework information ";
  openModal(message);

  console.log("Building map");

  //updateNavBar();

  jQuery("#lat").attr("placeholder", "Latitude").val("");   
  jQuery("#long").attr("placeholder", "Longitude").val("");   
  
  // Insert accordion text
  //
  jQuery.each(aboutFiles, function(keyItem, keyFile) {
  
      var InfoText = loadText(keyFile);
  
      jQuery("#" + keyItem).html(InfoText);
  
  });

  closeModal();

  // Set bounds based upon framework information
  //
  map_bounds = setBounds(myInfo);
		
  // Create the map object
  //
  map = new L.map('map', {zoomControl: false});

  map.fitBounds([
                 [map_bounds.lat_max, map_bounds.long_max], 
                 [map_bounds.lat_min, map_bounds.long_min]
                ]);

  map.setMaxBounds([
                 [map_bounds.lat_max, map_bounds.long_max], 
                 [map_bounds.lat_min, map_bounds.long_min]
                ]);

  // Disable controls for map view
  //
  //$(".leaflet-control-zoom").css("visibility", "visible");
  map.scrollWheelZoom.disable();    

  // Add home button
  //
  var zoomHome = L.Control.zoomHome();
  zoomHome.addTo(map);
	  
  // Add base map
  //
  map.addLayer(ESRItopoBasemap);
	  
  // Create the miniMap
  //
  miniMap = new L.Control.MiniMap(ESRItopoMinimap, { toggleDisplay: true }).addTo(map);

  // Add markerLayer and cross-section line
  //
  markerLayer  = new L.LayerGroup();
  profileLayer = new L.LayerGroup();
        
  // Add project boundary
  //
  //$.getJSON("ozark_line.geojson", function(data) {
  //  var boundaryLayer = L.geoJson(data, {
  //                                       style: function (feature) {
  //                                           return { color: "#f00", weight: 2, opacity: 1, fillOpacity: 0 };
  //                                      }
  //  });
  //  boundaryLayer.addTo(map).bringToBack();
  //  });
 
  // Add scale to map
  //
  //L.control.scale().addTo(map);
      
  // Add surficalGeology to base map
  //
  var imageUrl = 'grids/geotest_rgb.png';
  var imageBounds  = [
      [48.4194482, -121.8449579],
      [44.2608524, -115.3669151]
  ];
    
  // Doesn't preserve the cpras project surfical geology grids
  var imageUrl = 'gis/tt.tif';
  var imageBounds  = [
      [48.6095526509999942, -121.8449585989999946],
      [44.0062563269999956, -115.4507671469999934]
  ];
    
  var imageUrl = 'gis/geomap.tif';
  var imageBounds  = [
      [48.419427602, -121.844958599],
      [44.261075400, -115.367388015]
  ];
    
  var imageUrl = 'gis/test.tif';
  var imageBounds  = [
      [43.5438563270000003, -113.0552454900000043],
      [48.1471526509999990, -106.1331671469999947]
  ];
    
  var imageUrl = 'gis/notNeeded/t.wgs.png';
  var imageBounds  = [
      [48.4194482, -121.8449579],
      [44.2608524, -115.3669151]
  ];
    
  var imageUrl = 'gis/geomap.png';
  var imageBounds  = [
      [48.411817, -121.702958],
      [44.261075, -115.612214]
  ];

    
  var imageTxt = 'For Surface Geology see<a href="https://www.usgs.gov/publications/three-dimensional-model-geologic-framework-columbia-plateau-regional-aquifer-system">Burns and others (2010)</a>';


   // Surfical Geology overlay
   //
   var surficalGeology = L.imageOverlay(imageUrl,
                                        imageBounds,
                                        'Three-dimensional model of the geologic framework for the Columbia Plateau Regional Aquifer System, Idaho, Oregon, and Washington'
                                       ).addTo(map).bringToFront();        

   // Set initial opacity to 0.5 (Optional)
   //
   var opacityValue = sliderValue / 100;
   surficalGeology.setOpacity(opacityValue);	  
      
  // Slider for surfical geology
  //
  var controlGeologyOpacity = new L.Control.OpacitySlider(surficalGeology,
                                                              {
                                                                  sliderImageUrl: "css/usgs/images/opacity-slider.png",
                                                                  backgroundColor: "rgba(229, 227, 223, 0.9)",
                                                                  opacity: sliderValue,
                                                                  position: 'topright'
                                                              });
  controlGeologyOpacity.addTo(map);
  $('.opacity-control').prop('title','Move to change opacity of the surface geology on map');
  $('.opacity-control').prepend('<span>Surface Geology</span>');

  // Add zoom to your location
  //
      var myLocate = L.control.locate({
          drawCircle: false,
          drawMarker: false,
          clickBehavior: { outOfView: 'stop' },
          strings: { title: "Move and zoom to your location" }
      }).addTo(map);
//  jQuery('.leaflet-control-locate').on('click',(e)=>{
//      e.preventDefault();
//      e.stopPropagation();
//      console.log('Move and zoom to your location');
//      myLocate.stopFollowing;
//      var myAddress = { latlng: { lng: e.latlng.lng,  lat: e.latlng.lat } };
//
//      onMapClick(markerLayer, rasterLayer, myAddress);
//  });

  // Remove locations tool
  //
  var removeLocationTool = L.easyButton({
      id: 'removeLocationTool',
      states: [{
          icon: 'fa',
          title: 'Click to remove existing locations on map',
          onClick: function(control) {
              resetCrossSection(markerLayer, profileLayer);
          }
      }]
  }).addTo(map);
  $('#removeLocationTool').addClass('removeLocationTool');

  // Cross section tool
  //
  var xsecTool = L.easyButton({
      id: 'xsecTool',
      states: [{
          icon: 'fa',
          title: 'Click to make cross section between locations on map',
          onClick: function(control) {
              buildXsec(markerLayer, profileLayer);
          }
      }]
  }).addTo(map);
  $('#xsecTool').addClass('xsecTool');

  // Map bounds for geocoding tool
  //
  const corner1 = L.latLng(map_bounds.lat_max, map_bounds.long_ma);
  const corner2 = L.latLng(map_bounds.lat_min, map_bounds.long_mi);
  const bounds = L.latLngBounds(corner1, corner2);

  // Create the geocoding control and add it to the map
  //
  var searchControl = new L.esri.Controls.Geosearch({ zoomToResult: false, searchBounds: bounds }).addTo(map);
  jQuery('.geocoder-control').prop('title', "Enter address, intersection, or latitude/longitude");
  $(".geocoder-control").on("click", (e) => {
      e.preventDefault();
      e.stopPropagation();
  });

  searchControl.on('results', function(data) {
      console.log("geocoding results ",data);
  			
      if('latlng' in data)
        {
          var myAddress = { latlng: { lng: data.latlng.lng,  lat: data.latlng.lat } };

          //onMapClick(markerLayer, rasterLayer, myAddress);
          onMapClick(markerLayer, myAddress);
        }
  });
      
  // Base maps button from baseMaps.js
  //
  var basemapControl = L.control({position: 'topright'});
  basemapControl.onAdd = function (map) {
      var div = L.DomUtil.create('div', 'baseMaps NoJump');
      div.innerHTML = baseMapContent.join(" ");
      div.firstChild.onmousedown = div.firstChild.ondblclick = L.DomEvent.stopPropagation;
      return div;
  };
  basemapControl.addTo(map);      
  $('.baseMaps').prop('title','Click to change background on map');
      
  // Set current basemap active in dropdown menu
  //
  jQuery('#ESRItopoBasemap').addClass('active');
      
  // Clicked base map
  //
  $(".baseMaps").on("click", (e) => {
      e.preventDefault();
      e.stopPropagation();
  			
      // Clicked base map
      //
      jQuery('#basemapMenu li').click(function()
        {
         // Dont do anything if the current active layer is clicked
         //
         if(!map.hasLayer(window[$(this).prop("id")]))
           {			
            // Remove currently active basemap
            //
            jQuery("#basemapMenu li.active").each(function () 
               {
                console.log("removing ",$(this).prop("id"));

                map.removeLayer(window[$(this).prop("id")]);
                jQuery(this).removeClass("active");
               });
   			
            // Make new selection active and add to map
            //
            console.log("Active ",$(this).prop("id"));
            jQuery(this).addClass('active');
            map.addLayer(window[$(this).prop("id")]);
            miniMap.changeLayer(minimapObj[$(this).prop("id")]);
           }
   
         // Close after clicked
         //
         $("#basemapButton").dropdown("toggle");
      });
  });
 
  // Zoom message
  //
  $("#map").on("mouseover", function () {
      if(!myZoomFlag)
        {
         myZoomFlag = true
         message = "Use Shift-Left Mouse Drag: Select a region by pressing the Shift key and dragging the left mouse button"
         openModal(message);
         fadeModal(1000);
        }
     });
    
  // Show initial map zoom level
  //
  jQuery(".mapZoom" ).html( "<b>Zoom Level: </b>" + map.getZoom());

  // Refresh on extent change
  //
  map.on('zoomend dragend', function(evt) {
      var mapZoomLevel =  map.getZoom();
      if(mapZoomLevel > 12)
      {
          surficalGeology.setOpacity(0.0);
      }
      else
      {
          surficalGeology.setOpacity(sliderValue);
      }
      jQuery( ".mapZoom" ).html( "<b>Zoom Level: </b>" + map.getZoom());
      //jQuery( ".latlng" ).html(evt.latlng.lng.toFixed(3) + ", " + evt.latlng.lat.toFixed(3));
  });

  // Show mouse coordinates
  //
  map.on('mousemove', function(evt) {
	  jQuery( ".latlng" ).html(evt.latlng.lng.toFixed(3) + ", " + evt.latlng.lat.toFixed(3));
  });

  // Enable cell log on click
  //
  map.on('click', function(evt) {
     onMapClick(markerLayer, evt);
   });

  });

function userCoords() 
  {
    // Check coordinates
    //
    var user_longitude = jQuery("#long").prop("value");
    var user_latitude  = jQuery("#lat").prop("value");
  
    // if(user_latitude.length != 6 || user_latitude.indexOf('.') >= 0 || jQuery.isNumeric(user_latitude) === false)
    if(jQuery.isNumeric(user_latitude) === false && jQuery.isNumeric(user_longitude) === false)
      {
       jQuery("#lat").attr("placeholder", "Latitude").val("");   
       jQuery("#long").attr("placeholder", "Longitude").val("");   
       var message = "Please enter a six-digit number for latitude and longitude DMS (341098 degrees minutes seconds no spaces)";
       openModal(message);
       fadeModal(3000);
       return;
      }
  
    if(jQuery.isNumeric(user_latitude) === false)
      {
       jQuery("#lat").attr("placeholder", "Latitude").val("");   
       var message = "Please enter a six-digit number for latitude DMS (341098 degrees minutes seconds no spaces)";
       openModal(message);
       fadeModal(3000);
       return;
      }
    if(user_latitude.substr(0,1) == '-')
      {
       user_latitude = user_latitude.substr(1);
      }
  
    // if(user_longitude.length < 6 || user_longitude.length > 7 || user_longitude.indexOf('.') >= 0 || jQuery.isNumeric(user_longitude) === false)
    if(jQuery.isNumeric(user_longitude) === false)
      {
       jQuery("#long").attr("placeholder", "Longitude").val("");   
       var message = "Please enter a six-digit number for longitude DMS (902210 degrees minutes seconds no spaces)";
       openModal(message);
       fadeModal(3000);
       return;
      }
    if(user_longitude.substr(0,1) == '-')
      {
       user_longitude = user_longitude.substr(1);
      }

    // Translate map click location to model grid coordinates
    //
    //console.log("Coords " + user_longitude + " " + user_latitude);
    // var coords = DMS2dec(user_longitude, user_latitude);

    if(user_latitude.indexOf('.') == 0 || user_longitude.indexOf('.') == 0)
        {
            var coords = DMS2dec(user_longitude, user_latitude);
        }
    else
        {
            var coords = [];
            coords[0] = user_latitude;
            coords[1] = user_longitude;
        }

    //console.log("Coords " + coords[1] + " " + coords[0]);
    var lat    = coords[0];
    var long   = coords[1];
    if (long > 0){
        long = long * -1;
    }

    var coordinates            = User2User(
                                           { "x": long, "y": lat },
                                           latlong_projection,
                                           grid_projection
                                          );


    // Point inside model grid
    //
    if(isPointInPoly(framework, {x: long, y: lat}) > 0)
      {                                              
        // Place marker on map
        //
        marker = L.circle([lat, long], 
                          2, 
                          {
                           color: 'red',
                           fillColor: '#f03',
                           zIndexOffset: 999,
                           fillOpacity: 0.9
                          });
  
        // Add marker
        //
        markerLayer.addLayer(marker);
        map.addLayer(markerLayer);
 
        // Request cell log
        //
        var url = "well_log.html?";
        url    += "longitude=" + long;
        url    += "&latitude=" + lat;
        url    += "&x_coordinate=" + coordinates[0];
        url    += "&y_coordinate=" + coordinates[1];

        window.open(url);
      }
                                
    // Point outside model grid
    //
    else
      { 
       var message = "Point is located outside the boundaries of the Model grid";
       openModal(message);
       fadeModal(3000);
       return;
      }
}

function onMapClick(markerLayer, evt) 
  {
   if(typeof evt !== "undefined")
     {
      // Translate map click location to model grid coordinates
      //
      var long                   = evt.latlng.lng;
      var lat                    = evt.latlng.lat;
  
      var coordinates            = User2User(
                                             { "x": long, "y": lat },
                                             latlong_projection,
                                             grid_projection
                                            );
  
      // Point inside model grid
      //
      if(isPointInPoly(framework, {x: long, y: lat}) > 0)
        {                                              
          // Place marker on map
          //
          marker = L.circle([lat, long], 
                            2, 
                            {
                             color: 'red',
                             fillColor: '#f03',
                             zIndexOffset: 999,
                             fillOpacity: 0.9
                            });
    
          // Add marker
          //
          markerLayer.addLayer(marker);
          map.addLayer(markerLayer);
   
          // Request cell log
          //
          var url = "framework_cell_log.html?";
          url    += "longitude=" + long;
          url    += "&latitude=" + lat;
          url    += "&x_coordinate=" + coordinates[0];
          url    += "&y_coordinate=" + coordinates[1];
          //url    += "&color_file=" + color_file;
  
          //url    += "?" + (new Date()).getTime();
          //alert("URL " + url);
          window.open(url);
        }
                                  
      // Point outside model grid
      //
      else
        { 
         var message = "Point is located outside the boundaries of the Model grid";
         openModal(message);
         fadeModal(5000);
         return;
        }
    }
}

function updateLegend(v, currentMap, updateType) {

    //jQuery("#overlayMenu.sw").append('<li role="presentation" id="' + curSiteTypeInfo.overlayLayerName + '" class="' + curSiteTypeInfo.overlayLayerName  + '"><a role="menuitem" tabindex="-1"><div name="overlayLayers" ><img src="' + curSiteTypeInfo.singleMarkerURL + '"/><span>' + curSiteTypeInfo.legendLayerName + '</span></div></li>');}

}

function buildXsec(markerLayer, profileLayer)
  { 
   // Need two or more markers
   //
   if(markerLayer.getLayers().length < 2) 
     { 
      var message = "Need two or more locations in order to build a cross section";
      openModal(message);
      fadeModal(5000);
      return;
     }

   // Clear line segments
   //
   if(profileLayer.getLayers().length > 0)
      {
        map.removeLayer(profileLayer);
        profileLayer.clearLayers();
      }

   // Loop through markers
   //
   var i           = 0;
   var markers     = [];
   var points_txt  = [];
   var longlat_txt = [];

   markerLayer.eachLayer(function (layer)
                         {
                           var layer_id = layer._leaflet_id;
                           var color = "red";
                           i++;
                           if(i == 1 || i >=  markerLayer.getLayers().length )
                             {
                               color = "blue";
                             }

                           layer.setStyle({
                                           color: color,
                                           fillColor: color
                                          });
                           
                           //alert("markerLayer long " + layer.getLatLng().toString());
                           //alert("markerLayer long " + layer.getLatLng().lng);

                           var long = layer.getLatLng().lng;
                           var lat  = layer.getLatLng().lat;
                           markers.push([long, lat]);

                           var coordinates = User2User(
                                                       { "x": long, "y": lat },
                                                       latlong_projection,
                                                       grid_projection
                                                      );
                                    
                           points_txt.push( coordinates[0] + "," + coordinates[1] );
                           longlat_txt.push( long + "," + lat );
                         });

   // Add cross section line to map
   //
   var geojsonFeature = {
            "type": "Feature",
            "properties": {
                "ID": 1
            },
            "geometry": {
                "type": "LineString",
                "coordinates": markers
            }};

   var profile = L.geoJson(geojsonFeature, 
                            {
                              style: {
                                      color: "red",
                                      weight: 4
                                     }
                            });

   profileLayer.addLayer(profile);
   map.addLayer(profileLayer);
 
   // Request cell log
   //
   var url = "framework_xsec.html?";
   url    += "longlats=" + longlat_txt.join(" ");
   url    += "&points=" + points_txt.join(" ");
   //url    += "&color_file=" + color_file;
   //url    += "?" + (new Date()).getTime();
   //alert("URL " + url);
   window.open(url);
  }

function resetCrossSection(markerLayer, profileLayer)
  { 
   // Clear markers and line segments
   //
   if(markerLayer.getLayers().length > 0)
      {
       map.removeLayer(markerLayer);
       markerLayer.clearLayers();

       if(profileLayer.getLayers().length > 0)
          {
            map.removeLayer(profileLayer);
            profileLayer.clearLayers();
          }

       var message = "Removing all locations, please begin again. ";
       openModal(message);
       fadeModal(2000);
       return;
      }
  }
 
// Determine if a point is within Model grid
//
function isPointInPoly(poly, pt){
	for(var c = false, i = -1, l = poly.length, j = l - 1; ++i < l; j = i)
		((poly[i].y <= pt.y && pt.y < poly[j].y) || (poly[j].y <= pt.y && pt.y < poly[i].y))
		&& (pt.x < (poly[j].x - poly[i].x) * (pt.y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x)
		&& (c = !c);
	return c;
}


//      var surficalGeology = new L.tileLayer('geol_tiles/{z}/{x}/{y}.png', {
//          attribution: 'Map data &copy; ???',
//  }).addTo(map).bringToFront();
