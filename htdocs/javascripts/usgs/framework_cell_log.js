/**
 * Namespace: Framework_Cell_Log
 *
 * Framework_Cell_Log is a JavaScript library to build a single column of 
 *  framework information from the subsurface geologic layers.
 *
 * Special layout for MERAS project (addition line in explanation table for
 *  link to correlation web page CURRENTLY DISABLED).
 *
 * version 2.07
 * September 5, 2024
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
var frameworkFile      = "framework_parameters.js";

var map;
var configuration      = '{ "type": "map", "layers": { "group": { "name": "well log", "type": "marker", "details": "lat: -360, lng: -180", "latlng": [ -360, -180 ] } } }';

var map_options        = '{ "controls": { "zoom": false, "layers": { "autoZIndex": false, "sortLayers": false },';
map_options           += ' "attribution": { "prefix": "<a href= "https://github.com/gherardovarando/leaflet-map-builder"> leaflet-map-builder</a>" } },';
map_options           += ' "tooltip": { "marker": false, "polyline": false, "polygon": false, "circle": false, "rectangle": false }, ';
map_options           += ' "popup": { "marker": true, "polyline": false, "polygon": false, "circle": false, "rectangle": false     } }';

// Prepare when the DOM is ready 
//
$(document).ready(function() 
  {
   // Loading message
   //
   message = "Preparing cell log information";
   openModal(message);
        
    // Retrieve framework information
    //
    var myInfo       = getInfo(frameworkFile);
        
    var aboutFiles   = myInfo.wellLogFiles;
  
    // Insert accordion text
   //
    jQuery.each(aboutFiles, function(keyItem, keyFile) {
  
        var InfoText = loadText(keyFile);
  
        jQuery("#" + keyItem).html(InfoText);
  
    });

    // Parse url
    //
    longitude        = jQuery.url.param("longitude");
    latitude         = jQuery.url.param("latitude");
    x_coordinate     = jQuery.url.param("x_coordinate");
    y_coordinate     = jQuery.url.param("y_coordinate");
    //tiffs            = jQuery.url.param("tiffs");
    //color_file       = jQuery.url.param("color_file");

    var script_http  = myInfo.script_http;
    var rasters      = myInfo.rasters;
    var color_file   = myInfo.color_file;

    var northwest    = [ myInfo.northwest_x, myInfo.northwest_y ].join(",");
    var northeast    = [ myInfo.northeast_x, myInfo.northeast_y ].join(",");
    var southwest    = [ myInfo.southwest_x, myInfo.southwest_y ].join(",");
    var southeast    = [ myInfo.southeast_x, myInfo.southeast_y ].join(",");

    coordinates      = {
                        "longitude"    : longitude,
                        "latitude"     : latitude,
                        "x_coordinate" : x_coordinate,
                        "y_coordinate" : y_coordinate
                       };

    // Build cell log
    //
    createLog(coordinates, script_http, rasters, color_file);

  });
 
// Create Log
//
function createLog (coordinates, script_http, rasters, color_file) 
  {
    // Request for cell log information
    //
    var request_type = "GET";
    script_http      = "/cgi-bin/columbia/framework_well_log.pl";
    script_http      = "/cgi-bin/frameworkService/framework_location.py";
    var data_http    = "";
        data_http   += "longitude=" + coordinates.longitude;
        data_http   += "&latitude=" + coordinates.latitude;
        data_http   += "&x_coordinate=" + coordinates.x_coordinate;
        data_http   += "&y_coordinate=" + coordinates.y_coordinate;
        data_http   += "&color=" + color_file;
        data_http   += "&rasters=" + rasters;

    var dataType    = "json";

    webRequest(request_type, script_http, data_http, dataType, BuildCellGeometry);

  }
   
function showTooltip(name, x, y, contents) 
  {
    //alert("Tooltip " + contents + " at (" + x + ", " + y + ")");
    jQuery('<div id="tooltip">' + contents + '</div>').css( {
              top: y + 5,
              left: x + 5
          }).appendTo("#" + name).fadeIn(200);
  }

function locationMap(latitude, longitude)
  {
   //console.log("locationMap " + longitude + " " + latitude);

   var map = L.map('locationMap').setView([latitude, longitude], 8);
   map.scrollWheelZoom.disable();    

   L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}',
    {
      maxZoom: 17,
      minZoom: 9
    }).addTo(map);

   var color  = "#f06c00";
   var radius = 5;

   var circle = L.circleMarker([latitude, longitude], 
                               { 
                                radius: radius,
                                color: color,
                                fillColor: color,
                                fillOpacity: 0.15
   }).addTo(map);
      
  // Add surficalGeology to base map
  //
  var imageUrl = 'grids/geotest_rgb.png';
  var imageUrl = 'gis/geomap.tif';
  var imageTxt = 'For Surface Geology see<a href="https://www.usgs.gov/publications/three-dimensional-model-geologic-framework-columbia-plateau-regional-aquifer-system">Burns and others (2010)</a>';

  var imageBounds  = [
      [48.4194482, -121.8449579],
      [44.2608524, -115.3669151]
  ];

   // Surfical Geology overlay
   //
   var surficalGeology = L.imageOverlay(imageUrl,
                                        imageBounds
                                       ).addTo(map).bringToFront();        

   // Set initial opacity to 0.25 (Optional)
   //
   var opacityValue = 0.3;
   surficalGeology.setOpacity(opacityValue);	  

  }

function buildExplanation(cells)
  {
    console.log('buildExplanation');

    // Build explanation
    //
    var i = 0;
    var legend_html = [];
    legend_html.push('<table id="legend" class="cell_table">');
    legend_html.push(' <caption>' + "Explanation" + '</caption>');

    legend_html.push(' <thead class="fw-bold text-center">');
    legend_html.push(' <tr>');
    legend_html.push('  <th>Top<sup>1</sup></th>');
    legend_html.push('  <th>Bottom<sup>1</sup></th>');
    legend_html.push('  <th>Elevation</th>');
    legend_html.push('  <th>Thickness</th>');
    legend_html.push('  <th>Geologic Unit</th>');
    legend_html.push(' </tr>');
    legend_html.push(' </thead>');

    legend_html.push(' <tbody class="fw-medium">');

    // Set color specification array
    //
    console.log(cells);

    var top = 0;
    for(i = 0; i < cells.length; i++)
    {
      var unit        = cells[i].unit;
      var top_elev    = cells[i].top;
      var bot_elev    = cells[i].bottom;
      var top         = cells[i].depth;
      var thickness   = cells[i].thickness;
      var bot         = top + thickness;
      var color       = cells[i].color;
      var explanation = cells[i].explanation;
      var id          = unit
      var label       = id;

      console.log(`Unit ${unit} top ${top} bottom ${bot} color ${color} explanation ${explanation}`);

      // Build explanation
      //
      legend_html.push(` <tr id="${id}" bgcolor="${color}">`);

      legend_html.push(' <td>');
      legend_html.push(`  <div class="text-end"> ${Math.abs(top).toFixed(0)} </div>`);
      legend_html.push(' </td>');

      if(i < cells.length - 1)
      {
        legend_html.push(' <td>');
        legend_html.push(`  <div class="text-end"> ${Math.abs(bot).toFixed(0)} </div>`);
        legend_html.push(' </td>');
      }
      else
      {
        legend_html.push(' <td>');
        legend_html.push('  <div class="text-center">not determined</div>');
        legend_html.push(' </td>');
      }

      legend_html.push(' <td>');
      legend_html.push('  <div class="text-end">' + top_elev.toFixed(0) + '</div>');
      legend_html.push(' </td>');

      if(i < cells.length - 1)
      {
        legend_html.push(' <td>');
        legend_html.push('  <div class="text-end">' + Math.abs(thickness).toFixed(0) + '</div>');
        legend_html.push(' </td>');
      }
      else
      {
        legend_html.push(' <td>');
        legend_html.push('  <div class="text-center">not determined</div>');
        legend_html.push(' </td>');
      }

      legend_html.push(' <td style="background-color:#FFFFFF">');
      legend_html.push('  <div id="label_' + label + '" class="legendLabel">' + explanation + '</div>');
      legend_html.push(' </td>');

      legend_html.push(' </tr>');
    }

    legend_html.push(' <tr><td class="text-start" colspan="5"><sup>1</sup>' + "Values are depth below land surface, in feet" + '</td></tr>');
    //legend_html.push(' <tr><td id="cell_bottom" colspan="5"><a href="aq_names.html">Geologic correlation</a></td></tr>');
    legend_html.push(' </tbody>');
    legend_html.push('</table>');

    jQuery("#cell_table").html(legend_html.join("\n"));

  }

function BuildCellGeometry(json_data)
  { 
   console.log("BuildCellGeometry");
   console.log(json_data);

   // No subsurface
   //
   if(json_data.status != "success") 
     {
      var message = json_data.warning;
      if(typeof json_data.error !== "undefined") {message = json_data.error;}
      if(typeof json_data.warning !== "undefined") {message = json_data.warning;}

      openModal(message);
      fadeModal(3000);
      return;
     }

   // Check for returning warning or error
   //
   var message = json_data.warning;
   if(message) 
     {
      openModal(message);
      fadeModal(3000);
      return;
     }

   // Close modal dialog
   //
   closeModal();

   // General information
   //
   var cells             = json_data.cell_log;
   var longitude         = json_data.longitude;
   var latitude          = json_data.latitude;
   var easting           = json_data.easting;
   var northing          = json_data.northing;
   var row               = json_data.row;
   var col               = json_data.column;
   var rows              = json_data.nrows;
   var columns           = json_data.ncols;
   var time_units        = json_data.time_units;
   var length_units      = json_data.length_units;
   var cell_width        = json_data.cell_width;
   var x_axis_min        = json_data.x_axis_min;
   var x_axis_max        = json_data.x_axis_max;
   var y_axis_min        = json_data.y_axis_min;
   var y_axis_max        = json_data.y_axis_max;
   var cell_count        = json_data.cell_count;

   var LegendHash        = new Object();
   var ToolTipHash       = new Object();

   // No cell records
   //
   if(cell_count <= 0)
     {
      var warning  = "<p><b>No cell geometry for this row " + row + " and column " + col + "</b></p>";
      warning     += "<p><b>    All cells inactive</b></p>";
      openModal(message);
      return;
     }

   else
     {
       buildExplanation(cells);
     }

   // Page title
   //
   // var title = "Subsurface Information at Longitude " + parseFloat(longitude).toFixed(6) + " Latitude " + parseFloat(latitude).toFixed(6);
   var title = "Subsurface Information at Longitude " + parseFloat(longitude) + " Latitude " + parseFloat(latitude);
   jQuery(document).prop("title", title);
   jQuery("#page_title").html(title);

    // Cell log
    //
    plotCellLog(cells);

    // Location map
    //
    locationMap(latitude, longitude);
  }
