/**
 * Namespace: D3_Lithology
 *
 * D3_Lithology is a JavaScript library to provide a set of functions to build
 *  cell log column in svg format.
 *
 * version 2.03
 * September 13, 2024
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

// Set globals
//
var svg;
var jsonData;
var lithologyData;

var y_box_min   = 50;
var y_box_max   = 400;
var y_axis      = y_box_max - y_box_min;
var x_box_min   = 75;
var x_box_max   = x_box_min + 500;
var y_range;

var text_size   = 9;

// No information
//
function noLog(svgContainer)
  { 
   console.log("noLog");

   // No log label
   //
   label_txt       = "No Cell Cross Section Information";
   var  label      = "translate("
   label          += [( x_box_max + x_box_min ) * 0.5, + (y_box_max + y_box_min ) * 0.5].join(", ");
   label          += ") rotate(-90)";

   var myText      = svgContainer.append("text")
                                 .attr("transform", label)
                                 .attr('class', 'y_axis_label')
                                 .text(label_txt);
  }

// Plot Cell Log column
//
function plotCellXsec(cellData)
  {
   console.log("plotCellXsec");
   console.log(cellData);

   // Fade modal dialog
   //
   fadeModal(1000);
	
   // SVG canvas
   //
   var svg = d3.select("#xsec_graph");
	
   // Draw bore hole
   //
   axisBox(
           svg, 
           x_box_min, 
           x_box_max, 
           y_box_min, 
           y_box_max,
           "none"
          );

   // No cell log information
   //
   if(!cellData)
     {
      noLog(svg);
  
      return false;
     }
   else
     {
      return false;
     }
   
   // Plot specs
   //
   var min_value                  = 0.0;
   var max_value                  = null;
   var land_surface               = null;
   var min_elevation              = null;
   var max_elevation              = null;
   var min_depth                  = 0.0;
   var max_depth                  = null;
   var land_surface               = null;

   // Loop through cell records
   //
    var tempData     = cellData.slice();
    
    while ( tempData.length > 0 ) {
      var cellRecord  = tempData.shift();
      console.log(cellRecord);

      var unit        = cellRecord.unit;
      var top         = cellRecord.top;
      var bottom      = cellRecord.bottom;
      var top_depth   = cellRecord.depth;
      var thickness   = cellRecord.thickness;
      var bot_depth   = top_depth + thickness;
      var color       = cellRecord.color;
      var explanation = cellRecord.explanation;

      if(!min_elevation) { min_elevation = bottom; }
      else if(bottom < min_elevation) { min_elevation = bottom; }
      if(!max_elevation) { max_elevation = top; }
      else if(top > max_elevation) { max_elevation = top; }
      
      if(!max_depth) { max_depth = bot_depth; }
      else if(bot_depth > max_depth) { max_depth = bot_depth; }
    }
  
   // Plot specs
   //
   min_value                  = 0.0;
   max_value                  = max_elevation - min_elevation;
   land_surface               = max_elevation;
   returnList                 = get_max_min( min_value, max_value);
   y_min                      = returnList[0];
   y_max                      = returnList[1];
   y_interval                 = returnList[2];
   y_range                    = y_max - y_min
   if(y_min < 0) { y_min = 0.0; }
    console.log(`Y max ${y_max} Y min ${y_min} interval ${y_interval}`);
	
   // Add tooltip
   //
   var tooltip = addToolTip();

   // Geologic framework
   //
   addFramework(
                svg,
                x_box_min, 
                x_box_max, 
                y_box_min, 
                y_box_max, 
                y_interval, 
                cellData,
                tooltip
               )


   // Label axes
   //
   leftAxis(
            svg, 
            x_box_min, 
            x_box_max, 
            y_box_min, 
            y_box_max, 
            y_min, 
            y_max, 
            y_interval, 
            "Depth Below Land Surface, in feet"
           );

   var elevation_max = land_surface;
   var elevation_min = elevation_max - y_max;
   rightElevationAxis(
                      svg, 
                      x_box_min, 
                      x_box_max, 
                      y_box_min, 
                      y_box_max, 
                      elevation_min, 
                      elevation_max, 
                      y_interval, 
                      "Elevation, in feet"
                     );
  }

function addFramework(
                      svgContainer,
                      x_box_min, 
                      x_box_max, 
                      y_box_min, 
                      y_box_max, 
                      y_interval, 
                      cellData,
                      tooltip
                     )
  { 
   console.log("addFramework");
   //console.log(cellData);

   // Set
   //
   var y_range     = y_max - y_min;
   var y_axis      = y_box_max - y_box_min;

   // Loop through lithology
   //
   var tempData     = cellData.slice();
    
   while ( tempData.length > 0 ) {

        var cellRecord  = tempData.shift();
        console.log(cellRecord);

        var bottomFlag  = true;

        var id          = cellRecord.unit;
        var unit        = cellRecord.unit;
        var color       = cellRecord.color;
        var description = cellRecord.explanation;

        var top_depth   = cellRecord.depth;
        var thickness   = cellRecord.thickness;
        if(!thickness) { bottomFlag = false; }
        var bot_depth   = top_depth + thickness

        var width       = x_box_max - x_box_min

        var y_top       = y_box_min + y_axis * (top_depth - y_min) / y_range
        var y_bot       = y_box_min + y_axis * (bot_depth - y_min) / y_range
        var thickness   = y_bot - y_top
        if(tempData.length < 1) { thickness = y_box_max - y_top; }
     
        bot_depth = bot_depth.toFixed(0);
        if(!bottomFlag) { bot_depth = '"? ?"'; }
        var toolTip     = `${description}: Depth from ${top_depth.toFixed(0)} to ${bot_depth} feet`;
        var data        = [ {x:x_box_min, tooltip: toolTip}];

        // Add color
        //
           var lithology   = svgContainer.append("g")
                                         .attr("class", "lithology")
                                         .attr('id', id)
                                      .data(data)
           var myRect      = lithology.append("rect")
                                      .attr('x', x_box_min)
                                      .attr('y', y_top)
                                      .attr('width', width)
                                      .attr('height', thickness)
                                      .attr('fill', color)
                                   .on("mousemove", function(event, d) {
                                         tooltip
                                           .style("left", event.pageX + "px")
                                           .style("top", event.pageY + "px")
                                           .style("display", "inline-block")
                                           .html(d.tooltip);
                                   })
                                   .on("mouseout", function(d){ tooltip.style("display", "none");});

        if(tempData.length < 1) { 
          var myText = lithology.append("text")
              .text('?-?-?')
              .attr('x', x_box_min + 0.2 * (x_box_max - x_box_min))
              .attr('y', y_top + thickness)
        }
   }
  }

function addLegend(svgContainer, lithologyData, lithologyDefs)
  { 
   console.log("addLegend");

   var tempData     = lithologyDefs.slice();
   console.log(tempData);
  
   var x_legend     = x_box_max + 100
   var y_legend     = y_box_min
   var legend_box   = 20
   var y_top        = y_box_min

   var protocol     = window.location.protocol; // Returns protocol only
   var host         = window.location.host;     // Returns host only
   var pathname     = window.location.pathname; // Returns path only
   var url          = window.location.href;     // Returns full URL
   var origin       = window.location.origin;   // Returns base URL
   var webPage      = (pathname.split('/'))[1];

   console.log("protocol " + protocol);
   console.log("host " + host);
   console.log("pathname " + pathname);
   console.log("url " + url);
   console.log("origin " + origin);
   console.log("webPage " + webPage);

   var defs         = svgContainer.append("defs")

   // Loop through lithology
   //
   var Legend       = [];
   var LegendList   = [];
    
   while ( tempData.length > 0 ) {

        var lithRecord  = tempData.shift();

        var lithology   = lithRecord.lithology;
        var symbol      = lithRecord.symbol;
        var lithCode    = lithRecord.lithology.replace(/\s+&\s+/g, '');

        // Build legend
        //
        if(LegendList.indexOf(lithology) < 0)
          {
           var id          = symbol
           var svg_file    = symbol + ".svg"
           //var link_http   = [protocol + '/', host, webPage, "lithology_patterns", svg_file].join("/");
           var link_http   = ["lithology_patterns", svg_file].join("/");
   
           var pattern     = defs.append("pattern")
                                 .attr('id', id)
                                 .attr('patternUnits', 'userSpaceOnUse')
                                 .attr('width', 100)
                                 .attr('height', 100)
   
           var myimage     = pattern.append('image')
                                 .attr('xlink:href', link_http)
                                 .attr('width', 100)
                                 .attr('height', 100)
                                 .attr('x', 0)
                                 .attr('y', 0)

           LegendList.push(lithology);
           Legend.push({ 
                        'lithCode': lithCode,
                        'symbol': symbol,
                        'description': lithology,
                        'image': id
                       })

           //lithologyDefs[lithology].pattern = id;
          }
   }

   // Loop through lithology
   //
   var tempData     = Legend;
  
   var x_legend     = x_box_max + 100
   var y_legend     = y_box_min
   var legend_box   = 20
   var y_top        = y_box_min

   var descriptions = svgContainer.append("g")
                                  .attr("class", "legend_descriptions")
    
    while ( tempData.length > 0 ) {

        var Record      = tempData.shift();

        var lithCode    = Record.lithCode;
        var symbol      = Record.symbol;
        var description = Record.description
        var id          = Record.image
        var url         = 'url(#' + id + ')'

        var myRect      = descriptions.append("rect")
                                      .attr('x', x_legend)
                                      .attr('y', y_top)
                                      .attr('width', legend_box)
                                      .attr('height', legend_box)
                                      .attr('fill', url)
                                      .attr('stroke', 'black')
                                      .attr('stroke-width', 1)

        var myText      = descriptions.append("text")
                                      .text(description)
                                      .attr('class', lithCode)
                                      .attr('x', x_legend + legend_box * 1.25)
                                      .attr('y', y_top + legend_box * 0.5)
                                      .on('mouseover', function(d, i) {
                                         var lithClass = d3.select(this).attr('class');
                                         d3.selectAll("#" + lithClass)
                                           .transition()
                                           .duration(100)
                                           .attr('strokeWidth', 10)
                                           .attr('stroke', 'yellow')
                                      })
                                      .on('mouseout', function(d, i) {
                                         var lithClass = d3.select(this).attr('class');
                                         d3.selectAll("#" + lithClass)
                                           .transition()
                                           .duration(100)
                                           .attr('strokeWidth', 1)
                                           .attr('stroke', 'black')
                                      })

        y_top          += legend_box * 1.5
   }
  
   console.log("done addLegend");
  }

// Min and max
//
function get_max_min( min_value, max_value)
  { 
   var factor         = 0.01; 
   var interval_shift = 0.67; 
   var range          = max_value - min_value; 
        
   var interval       = factor; 
   range              = range / 5.0; 
        
   // Determine interval 
   // 
   while (range > factor) 
     { 
      if(range <= (factor * 1)) 
        { 
   	 interval = factor * 1; 
        } 
      else if (range <= (factor * 2))
        { 
   	 interval = factor * 2; 
        } 
      else if (range <= (factor * 2.5))
        { 
   	 if(factor < 10.0) 
           { 
            interval = factor * 2; 
           } 
         else 
           { 
            interval = factor * 2.5; 
           } 
        } 
      else if (range <= (factor * 5))
        { 
         interval = factor * 5;
        } 
      else
        { 
         interval = factor * 10;
        } 

       factor = factor * 10; 
    } 

   // Maximum
   //
   factor = parseInt(max_value / interval); 
   value  = factor * interval; 
   if(max_value >= value ) 
     { 
      value += interval;
     } 
   else if(Math.abs(max_value - value) <= interval_shift * interval) 
     { 
      max_value = value + interval; 
     } 
   else 
     { 
      max_value = value; 
     } 

   // Minimum
   //
   factor = parseInt(min_value / interval); 
   value  = factor * interval; 
   if(min_value >= value ) 
     { 
      value = (factor - 1) * interval; 
     } 
   if(Math.abs(min_value - value) <= interval_shift * interval) 
     { 
      min_value = value - interval; 
     } 
   else 
     { 
      min_value = value; 
     } 
      
   return [min_value, max_value, interval];
  }

