/**
 * Namespace: Framework_Tools
 *
 * Framework_Tools is a JavaScript library to read framework information such as
 * projection and boundary information.
 *
 * version 2.02
 * July 14, 2023
*/

/*
###############################################################################
# Copyright (c) Leonard Orzol Oregon Water Science Center
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

// Prevent jumping to top of page when clicking a href
//
jQuery('.noJump a').click(function(event){
   event.preventDefault();
});
   
// Request Framework information from json file
//
function getInfo(myFile) 
  {
    //alert("Grabbing framework parameters from " + myFile);

    var myInfo      = {};

    $.support.cors  = true;
    jQuery.ajax( 
                { url: myFile + "?_="+(new Date()).valueOf(),
                  dataType: "json",
                  async: false
                })
      .done(function(data)
            {
              myInfo  = parseInfo(data);
            })
      .fail(function() { 
          var message = "Error reading framework configuration file " + myFile;
          console.log(message);
          message_dialog("Error", message);
        });

    return myInfo;

  }

function parseInfo(json) 
  {
    return json;               
  }

function setBounds(myInfo) 
  {
    // Set projection from Grid coordinates to LatLong coordinates
    // 
    grid_projection            = myInfo.framework_projection;
    latlong_projection         = myInfo.latlong_projection;
                                                                 
    // Northwest corner
    // 
    grid_northwest             = myInfo.northwest_x + "," + myInfo.northwest_y;
    grid_x_origin              = myInfo.northwest_x;
    grid_y_origin              = myInfo.northwest_y;

    var latlong_coordinate     = proj4(grid_projection, latlong_projection, [grid_x_origin, grid_y_origin]);

    var long_nw                = latlong_coordinate[0];
    var lat_nw                 = latlong_coordinate[1];
                                                                 
    var framework_min_x        = grid_x_origin;
    var framework_max_x        = grid_x_origin;
    var framework_min_y        = grid_y_origin;
    var framework_max_y        = grid_y_origin;
                                                                 
    // Northeast corner
    // 
    grid_northeast             = myInfo.northeast_x + "," + myInfo.northeast_y;
    grid_x_northeast           = myInfo.northeast_x;
    grid_y_northeast           = myInfo.northeast_y;

    var latlong_coordinate     = proj4(grid_projection, latlong_projection, [grid_x_northeast, grid_y_northeast]);

    var long_ne                = latlong_coordinate[0];
    var lat_ne                 = latlong_coordinate[1];
                                                                 
    if(grid_x_northeast > framework_max_x) { framework_max_x = grid_x_northeast; }
    if(grid_x_northeast < framework_min_x) { framework_min_x = grid_x_northeast; }
    if(grid_y_northeast > framework_max_y) { framework_max_y = grid_y_northeast; }
    if(grid_y_northeast < framework_min_y) { framework_min_y = grid_y_northeast; }
                                                                 
    // Southeast corner
    // 
    grid_southeast             = myInfo.southeast_x + "," + myInfo.southeast_y;
    grid_x_southeast           = myInfo.southeast_x;
    grid_y_southeast           = myInfo.southeast_y;

    var latlong_coordinate     = proj4(grid_projection, latlong_projection, [grid_x_southeast, grid_y_southeast]);

    var long_se                = latlong_coordinate[0];
    var lat_se                 = latlong_coordinate[1];
                                                                 
    if(grid_x_southeast > framework_max_x) { framework_max_x = grid_x_southeast; }
    if(grid_x_southeast < framework_min_x) { framework_min_x = grid_x_southeast; }
    if(grid_y_southeast > framework_max_y) { framework_max_y = grid_y_southeast; }
    if(grid_y_southeast < framework_min_y) { framework_min_y = grid_y_southeast; }
                                                                 
    // Southwest corner
    // 
    grid_southwest             = myInfo.southwest_x + "," + myInfo.southwest_y;
    grid_x_southwest           = myInfo.southwest_x;
    grid_y_southwest           = myInfo.southwest_y;

    var latlong_coordinate     = proj4(grid_projection, latlong_projection, [grid_x_southwest, grid_y_southwest]);

    var long_sw                = latlong_coordinate[0];
    var lat_sw                 = latlong_coordinate[1];
                                                                 
    if(grid_x_southwest > framework_max_x) { framework_max_x = grid_x_southwest; }
    if(grid_x_southwest < framework_min_x) { framework_min_x = grid_x_southwest; }
    if(grid_y_southwest > framework_max_y) { framework_max_y = grid_y_southwest; }
    if(grid_y_southwest < framework_min_y) { framework_min_y = grid_y_southwest; }
                                                  
    // Framework boundary
    // 
    framework = [
                  {x: long_nw, y: lat_nw}, 
                  {x: long_ne, y: lat_ne}, 
                  {x: long_se, y: lat_se}, 
                  {x: long_sw, y: lat_sw}, 
                  {x: long_nw, y: lat_nw}
                 ]; 
                                                                 
    // Framework grid center
    // 
    var framework_x            = ( framework_min_x + framework_max_x ) * 0.50;
    var framework_y            = ( framework_min_y + framework_max_y ) * 0.50;

      //console.log(framework_x)
      //console.log(framework_y)

    var latlong_coordinate     = proj4(grid_projection, latlong_projection, [framework_x, framework_y]);

      //console.log(latlong_coordinate)

    var framework_ctr          = {
                                  long: latlong_coordinate[0], 
                                  lat : latlong_coordinate[1], 
                                 };

    var p_min                  = proj4(grid_projection, latlong_projection, [framework_min_x, framework_min_y]);
    var p_max                  = proj4(grid_projection, latlong_projection, [framework_max_x, framework_max_y]);

    framework_bnds             = {
                                  long_min: p_min[0], lat_min: p_min[1], 
                                  long_max: p_max[0], lat_max: p_max[1] 
                                 };

    return framework_bnds;               
  }

// Load text
//
function loadText(file_name) 
  {
    var myInfo = "";

    // Check file name
    //
    if(file_name.length < 1)
      {
        var message = "No file specified";
        openModal(message);
        fadeModal(4000);
        return;
      }

    // Load file
    //
    jQuery.ajax( 
                { url: file_name + "?_="+(new Date()).valueOf(),
                  dataType: "text",
                  async: false
                })
      .done(function(data)
            {
              myInfo = data;
            })
      .fail(function() 
            { 
              var message = "No file specified";
              openModal(message);
              fadeModal(4000);
              return;
            });

    return myInfo;
  }
