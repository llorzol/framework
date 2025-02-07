/**
 * Namespace: Framework_Xsec
 *
 * Framework_Xsec is a JavaScript library to build a profile of columns of 
 *  framework information from the subsurface geologic layers.
 *
 * version 2.01
 * July 19, 2023
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
var frameworkFile  = "framework_parameters.js";

var null_color     = "#D8D8D8";

// Prepare when the DOM is ready 
//
$(document).ready(function() 
  {
   // Retrieve framework information
   //
   var myInfo       = getInfo(frameworkFile);

   var aboutFiles   = myInfo.xsecFiles;

   // Insert accordion text
   //
   jQuery.each(aboutFiles, function(keyItem, keyFile) {
  
        var InfoText = loadText(keyFile);
  
        jQuery("#" + keyItem).html(InfoText);
  
    });
        
    // Parse url
    //
    var longlats     = jQuery.url.param("longlats");
    var points       = jQuery.url.param("points");

    // Retrieve model and modflow information
    //
    var myInfo       = getInfo(frameworkFile);

    var script_http  = myInfo.script_cgi_bin;
    var rasters      = myInfo.rasters;
    var color_file   = myInfo.color_file;

    coordinates      = {
                        "longlats"    : longlats,
                        "points"      : points
                       };

    // Build cross section
    //
    createXsec(coordinates, script_http, rasters, color_file);

  });
 
// Create cross section
//
function createXsec (coordinates, script_http, rasters, color_file) 
  {

    // Request for cell log information
    //
    var request_type = "GET";
    var script_http  = [script_http, "framework_xsec.pl?"].join("/");
    var script_http  = '/cgi-bin/frameworkService/framework_xsec.py';
    var data_http    = "";
        data_http   += "longlats=" + coordinates.longlats;
        data_http   += "&points=" + coordinates.points;
        data_http   += "&rasters=" + rasters.join(" ");
        data_http   += "&color=" + color_file;

    var dataType    = "json";

    webRequest(request_type, script_http, data_http, dataType, BuildCellXsec);

  }

function YtickFormatter(v, axis)
  { 
   vv = v ; 
   return vv;
  }

function roundNumber(num, dec) {
	var result = Math.round(num*Math.pow(10,dec))/Math.pow(10,dec);
	return result;
}

function Y2tickFormatter(v, axis)
  { 
   vv = v;
   if(y_axis_interval > 2)
     {
      vv = roundNumber(vv, 0)
     }
   else
     {
      vv = roundNumber(vv, 2)
     }

   return vv;
  }

function BuildCellXsec(json_data)
  { 
   // No subsurface
   //
   if(json_data.status != "success") 
     {
      var message = json_data.warning;
      if(typeof json_data.error !== "undefined") {message = json_data.error;}
      if(typeof json_data.warning !== "undefined") {message = json_data.warning;}
      //alert(" Retrieved data " + message);

      openModal(message);
      fadeModal(3000);
      return;
     }

   // Check for returning warning or error
   //
   var message       = json_data.warning;
   if(typeof message !== "undefined") 
     {
      openModal(message);
      fadeModal(3000);
      return;
     }
   
   // General information
   //
   rows              = json_data.nrows;
   columns           = json_data.ncols;
   layers            = json_data.nlays;
   cell_width        = json_data.cell_width;
   x_axis_min        = json_data.x_axis_min;
   x_axis_max        = json_data.x_axis_max;
   y_axis_min        = json_data.y_axis_min;
   y_axis_max        = json_data.y_axis_max;
   cell_count        = json_data.cell_count;
   start_row         = json_data.start_row;
   start_col         = json_data.start_col;
   end_row           = json_data.end_row;
   end_col           = json_data.end_col;

   var LegendHash    = new Object();
   var ToolTipHash   = new Object();
   var data          = [];
   var data_sets     = [];
   var color_scheme  = [];

   var updateLegendTimeout = null;
   var latestPosition      = null;
   var previousPoint       = null;
   var placeholder_offset  = null;
   var placeholder_width   = null;
                                                                         
   // Data
   //
   data_sets     = [];
                                              
   // Rocks
   //
   //alert("Rocks " + json_data.cell_fields);
   var color_scheme = new Array();
               
   jQuery.each(json_data.rocks, 
               function(index) {
                           var array       = json_data.rocks[index].array;
                           var unit        = json_data.rocks[index].unit;
                           var color       = json_data.rocks[index].color;
                           var explanation = json_data.rocks[index].explanation;
                           var id          = color.replace("#","") 
                           var label       = id;
                           
                           color_scheme.push(color);
                                                                         
                           data_sets.push(
                                       { 
                                        label: label,
                                        lines: { show :      true, 
                                                 lineWidth : 0,
                                                 fill :      1.0
                                               },
                                        id   : id,
                                        data : array
                                       });
                                   
                           if(typeof LegendHash[label] === "undefined")
                             {
                              LegendHash[label] = 1;
                             }
                          });
   
   // Hard-code color indices to prevent them from shifting as
   // elements are turned on/off
   var i = 0;
   jQuery.each(data_sets, function(key, val) {
       val.color = i;
       //val.color = val.id;
       ++i;
   });

   // Title
   //
   var longlats = jQuery.url.param("longlats");
   longitude    = longlats.split(" ")[0].split(",")[0];
   latitude     = longlats.split(" ")[0].split(",")[1];

   var title     = "Geologic Framework Information starting at Longitude, Latitude";
   title        += " (" + parseFloat(longitude).toFixed(4) + ", " + parseFloat(latitude).toFixed(4) + ")";
   longitude    = longlats.split(" ")[1].split(",")[0];
   latitude     = longlats.split(" ")[1].split(",")[1];

   title        += " to (" + parseFloat(longitude).toFixed(4) + ", " + parseFloat(latitude).toFixed(4) + ")";
   jQuery(document).prop("title", title);
   jQuery("#page_title").html(title);
   
   // Data
   //
   data = data_sets;
    
   // Change title
   // 
   jQuery('#cellGraph').show();
   
   // Build explanation
   // 
   var i = 0;
   var legend_html = [];
   legend_html.push('<table id="legend" class="legendTableborder="1">');
   legend_html.push(' <tbody>');
   legend_html.push(' <tr><td id="explanation" colspan="4">' + "Geologic units" + '</tr></td>');
   jQuery.each(json_data.explanations, 
               function(index) {
                                var unit        = json_data.explanations[index].unit;
                                var color       = json_data.explanations[index].color; 
                                var explanation = json_data.explanations[index].explanation;
                                var label       = color.replace("#","");

                                if(typeof LegendHash[label] !== "undefined")
                                  {
                                   ToolTipHash[label] = explanation;

                                   legend_html.push(' <tr>');

                                   legend_html.push(' <td>');
                                   legend_html.push('  <div id="value_' + label + '" class="legendValue">&nbsp;&nbsp;&nbsp;</div>');
                                   legend_html.push(' </td>');

                                   legend_html.push(' <td class="checkBoxes">');
                                   legend_html.push('  <input type="checkbox" name="checkbox_' + label + '" value="on" id="' + label + '" checked="checked" />');
                                   legend_html.push(' </td>');

                                   legend_html.push(' <td>');
                                   legend_html.push('   <div id="legend_' + label + '" class="legendColorBox" style="border: solid ' + color + '">&nbsp;</div>');
                                   legend_html.push(' </td>');

                                   legend_html.push(' <td>');
                                   legend_html.push('  <div id="label_' + label + '" class="legendLabel">' + explanation + '</div>');
                                   legend_html.push(' </td>');

                                   legend_html.push(' </tr>');
                                   
                                   LegendHash[label] = 1;
                                  }

   });
   legend_html.push(' </tbody>');
   legend_html.push('</table>');
   jQuery("#cell_legend").html(legend_html.join(""));
                                                                     
   // Graph options
   //
   options = {
              legend: { 
                       show: false 
                      },
              lines:  { show: true, fill: 1.0 },
              grid:   { 
                       aboveData: true,
                       borderWidth: 2.0,
                       hoverable: true, 
                       clickable: true,
                       hoverFill: '#444',
                       hoverRadius: 5
                      },
              crosshair: { mode: "x" },
              yaxis:  { 
                       label: "Elevation",
                       labelWidth : 50,
                       ticks: 5, 
                       max: y_axis_max,
                       min: y_axis_min
                      },
              xaxis:  {
                       ticks: 5, 
                       min: x_axis_min,
                       max: x_axis_max
                      },
              tooltip: false,
              tooltipOpts: {
                            content: "%s for %x was %y"
                           },
              selection: { mode: "xy" },
              colors: color_scheme
             };
                                                                     
                                                                     
   // Graph options
   //
   reset_options = {
                  legend: { 
                           show: false 
                          },
                  lines: { show: true, fill: 1.0 },
                  grid: { 
                         aboveData: true,
                         borderWidth: 2.0,
                         hoverable: true, 
                         clickable: true,
                         hoverFill: '#444',
                         hoverRadius: 5
                        },
                  crosshair: { mode: "x" },
                  yaxis: { 
                          label: "Elevation",
                          labelWidth : 50,
                          ticks: 5, 
                          max: y_axis_max,
                          min: y_axis_min
                         },
                  xaxis: {
                          min: x_axis_min,
                          max: x_axis_max
                         },
                  tooltip: false,
                  tooltipOpts: {
                                content: "%s for %x was %y"
                               },
                  selection: { mode: "xy" },
                  colors: color_scheme
                 };
                                                                     
   // Graphs 
   //
   plot     = jQuery.plot(jQuery("#xsec_graph"), data_sets, options);
                                           
   overview = jQuery.plot(jQuery("#xsec_overview"), data_sets, {
                                                                legend: { 
                                                                         show: false 
                                                                        },
                                                                lines: { show: true, fill: 1.0 },
                                                                grid: { aboveData: true, color: "#000000", borderWidth: 2.0 },
                                                                yaxis: { 
                                                                        ticks: 3, 
                                                                        label: "Elevation",
                                                                        labelWidth : 50,
                                                                        min: y_axis_min,
                                                                        max: y_axis_max
                                                                       },
                                                                xaxis: {
                                                                        min: x_axis_min,
                                                                        max: x_axis_max
                                                                       },
                                                                selection: { mode: "xy" },
                                                                colors: color_scheme
                                                               });
                                                        
   jQuery("#xsec_overview").dblclick(function () { 
       var axes = plot.getAxes();
       axes.xaxis.min = json_data.x_axis_min;
       axes.xaxis.max = json_data.x_axis_max;
       axes.yaxis.min = json_data.y_axis_min;
       axes.yaxis.max = json_data.y_axis_max;
       jQuery.plot(jQuery("#xsec_graph"), 
                   data_sets,
                   jQuery.extend(true, {}, options, {
                                                     xaxis: { min: axes.xaxis.min, max: axes.xaxis.max },
                                                     yaxis: { min: axes.yaxis.min, max: axes.yaxis.max }
                     }));
     });
                                                    
   // Connect plot and overview
   //
   jQuery("#xsec_graph").bind("plotselected", function (event, ranges) {
        // clamp the zooming to prevent eternal zoom
        var delta_x = ranges.xaxis.to - ranges.xaxis.from;
        if (delta_x < 0.00001) {
            ranges.xaxis.to = ranges.xaxis.from + 0.00001;
        }
        var delta_y = ranges.yaxis.to - ranges.yaxis.from;
        if ( delta_y < 0.00001 ) {
            ranges.yaxis.to = ranges.yaxis.from + 0.00001;
        }
                                               
        // do the zooming
        plot = jQuery.plot(jQuery("#xsec_graph"), 
                           data_sets,
                           jQuery.extend(true, {}, options, {
                                                             xaxis: { min: ranges.xaxis.from, max: ranges.xaxis.to },
                                                             yaxis: { min: ranges.yaxis.from, max: ranges.yaxis.to }
                                                            }));
                                               
        // don't fire event on the overview to prevent eternal loop
        overview.setSelection(ranges, true);
       });
       
   // Show sleected area on overview
   //
   jQuery("#xsec_overview").bind("plotselected", 
                                 function (event, ranges) {
                                                           plot.setSelection(ranges);
                                                          });
                
   // Tooltip
   //
   jQuery("#xsec_graph").bind("plothover", function (event, pos, item) {
    
        // Check to see if axes are enabled
        //
        if(typeof pos.x === "undefined") { return; }
        if(typeof pos.y === "undefined") { return; }
    
        // Update value in legend
        //
        latestPosition = pos;
        if (!updateLegendTimeout)
          updateLegendTimeout = setTimeout(updateLegend, 50);

        //jQuery("#x").text(pos.x.toFixed(2));
        //jQuery("#y").text(pos.y.toFixed(2));

       // Tooltip window enabled
       //
       if(jQuery("#enableTooltip").is(':checked'))
         {
           if(item) 
             {
               if(previousPoint != item.datapoint) 
                 {
                   previousPoint  = item.datapoint;
                   
                   jQuery("#tooltip").remove();
                        
                   var x = item.datapoint[0].toFixed(0);
                   var y = item.datapoint[1].toFixed(0);
                         
                   var contents = ToolTipHash[item.series.label] + " Distance " + x + " Altitude " + y;
                        
                   var label_offset  = jQuery("#xsec_graph").offset();
                   var labelX_offset = label_offset.left;
                   var labelY_offset = label_offset.top;

                   x = item.pageX - labelX_offset;
                   y = item.pageY - labelY_offset;
                   
                   showTooltip("xsec_graph", x, y, contents);
                 }
             }
           else 
             {
               jQuery("#tooltip").remove();
               previousPoint = null;            

               jQuery("#distance_legend").text("Distance " + pos.x.toFixed(0));
             }
          }

       // Tooltip window disabled
       //
       else 
          {
            jQuery("#tooltip").remove();

            jQuery("#distance_legend").text("Distance " + pos.x.toFixed(0));
          }
   });
    
    function updateLegend() {
              
        updateLegendTimeout = null;
        
        var pos  = latestPosition;
        
        // Check position with active axes
        //
        var axes = plot.getAxes();

        if(typeof axes === "undefined")
           return;
                 
        if(pos.x < axes.xaxis.min || pos.x > axes.xaxis.max)
          {
            //alert("X axis min " + axes.xaxis.min + " max " + axes.xaxis.max + " X pos " + pos.x);
           return;
          }

        if(typeof pos.y  !== "undefined") 
          { 
           if(pos.y < axes.yaxis.min || pos.y > axes.yaxis.max)
              return;
          }

        // Determine value against active axes
        //
        var i, j, dataset = plot.getData();
        for(i = 0; i < dataset.length; ++i) 
           {
            var series  = dataset[i];
            var label   = series.label;
            var id      = series.id;
            var axis_nu = series.yaxis.n;
            if(LegendHash[id] > 0)
              {
                // Find the nearest points, x-wise
                //
                for (j = 0; j < series.data.length; ++j)
                    if (series.data[j][0] > pos.x)
                        break;
                
                // Interpolate
                //
                var y_str = "--";
                var p1 = series.data[j - 1];
                var p2 = series.data[j];
                if(typeof p1 === "undefined" || typeof p2 === "undefined") 
                   y = " ";
                else if (p1[1] === null)
                   y = " ";
                else if (p2[1] === null)
                   y = " ";
                else
                  {
                   y = p1[1];
                   var x1_diff = Math.abs(pos.x - p1[0]);
                   var x2_diff = Math.abs(pos.x - p2[0]);
                   if(x2_diff < x1_diff)
                      y = p2[1];
                   //y = p1[1] + (p2[1] - p1[1]) * (pos.x - p1[0]) / (p2[0] - p1[0]);
                   y_str = y.toString();
                  }
    
                jQuery('#value_' + id).html(y_str);
              }
            else
              {
                jQuery('#value_' + id).html("--");
              }
           }
    }

   // Single click
   //
   jQuery("#xsec_graph").bind("plotclick", function (event, pos, item) {
       if (item) {
           jQuery("#clickdata").text("You clicked in " + item.series.label + ".");
           plot.highlight(item.series, item.datapoint);
       }
       else {
           jQuery("#clickdata").text("");
       }
   });


   // Legend
   //
   $(':checkbox').click(function() 
                                {
                                 //var selected   = jQuery(this).attr('name');
                                 var id         = jQuery(this).attr('id');
                                 var isChecked  = jQuery(this).attr('checked');
                                 var name       = jQuery(this).attr('name');
                                 var legendshow = LegendHash[id];
                                 //alert("Legend id " + id + " name " + name + " checked " + isChecked + " legend " + legendshow);
                                    
                                 // Remove data set
                                 //
                                 if(legendshow > 0)
                                   {
                                    LegendHash[id] = 0;
                                    jQuery('#label_' + id).css({ "opacity": 0.4 });
                                    jQuery('#legend_' + id).css({ "opacity": 0.4 });
                                   }
                                 else
                                   {
                                    LegendHash[id] = 1;
                                    jQuery('#label_' + id).css({ "opacity": 1.0 });
                                    jQuery('#legend_' + id).css({ "opacity": 1.0 });
                                   }
                                 
                                 data = [];
                                 for(i = 0; i < data_sets.length; i++)
                                    {
                                     var id = data_sets[i].id;
                                     if(LegendHash[id] > 0)
                                       {
                                        data.push(data_sets[i]);
                                       }
                                    }
        
                                 // Determine min/max of xaxis
                                 //
                                 var axes = plot.getAxes();
                             
                                 jQuery.plot(jQuery("#xsec_graph"), 
                                             data,
                                             jQuery.extend(true, {}, options, {
                                                                               xaxis: { min: axes.xaxis.min, max: axes.xaxis.max },
                                                                               yaxis: { min: axes.yaxis.min, max: axes.yaxis.max }
                                                                              }
                                  ));
                               });
  }
