/**
 * Namespace: Flot_Utils
 *
 * Flot_Utils is a JavaScript functions to help with plotting of Modflow
 * input and output arrays.
 *
 * version 1.16
 * January 4, 2015
*/

/*
###############################################################################
# Copyright (c) 2015, Leonard Orzol Oregon Water Science Center
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

function findPos(obj) {
    var curleft = 0, curtop = 0;
    if (obj.offsetParent) {
        do {
            curleft += obj.offsetLeft;
            curtop  += obj.offsetTop;
        } while (obj = obj.offsetParent);
        return { x: curleft, y: curtop };
    }
    return undefined;
}

function rgbToHex(r, g, b) {
    if (r > 255 || g > 255 || b > 255)
        throw "Invalid color component";
    return ((r << 16) | (g << 8) | b).toString(16);
}

function hexToR(h) {return parseInt((cutHex(h)).substring(0,2),16)}
function hexToG(h) {return parseInt((cutHex(h)).substring(2,4),16)}
function hexToB(h) {return parseInt((cutHex(h)).substring(4,6),16)}
function cutHex(h) {return (h.charAt(0)=="#") ? h.substring(1,7):h}

function YtickFormatter(v, axis)
  { 
   var vv = v ; 
   return vv;
  }

function roundNumber(num, dec) 
{
	var result = Math.round(num*Math.pow(10,dec))/Math.pow(10,dec);
	return result;
}

function Y2tickFormatter(v, axis)
  { 
   var vv = v;
   if(y_axis_interval > 2)
     {
       vv = roundNumber(vv, 0);
     }
   else
     {
       vv = roundNumber(vv, 2);
     }

   return vv;
  }

function ConfigureModflowInput(action, configure_file, array_file, lookup_file)
  { 
   //alert("Configure Modflow Input");
   
   var form_text  = "";
       form_text += "<span id=\"file_input\">Add Modflow array file<input id=\"array_file\" name=\"array_file\" type=\"file\" \/><\/span>";
       form_text += "<script>";
       form_text += "$(function () {";
       form_text += "    $('#file_input').mouseover(function() {";
       form_text += "        $(this).addClass('over');";
       form_text += "        }).mouseout(function() {";
       form_text += "                $(this).removeClass('over');";
       form_text += "        });";
       form_text += "        $('#file_input input:file').change(function() {";
       form_text += "                var filepath = $(this).val();";
       form_text += "                alert(\"you choose: \" + filepath);";
       form_text += "        });";
       form_text += "});";
       form_text += "<\/script>";
        
       cell_dialog(form_text,"Modflow Input Configuration");

   return;
  }
        
// Tool tip 
//
function showTooltip(div_id, x, y, x_offset, y_offset, contents) 
  {
    var x_position = x + x_offset;
    var y_position = y + y_offset;
    jQuery("#tooltip").remove();
    jQuery('<div id="tooltip" class="flot_tooltip">' + contents + '</div>').css( {
          top:  parseInt(y_position),
          left: parseInt(x_position)
          }).appendTo(div_id).fadeIn(200);
  }

function SetCellGeometry(cell_status)
  { 
   if(cell_status == "undefined")
     { 
      cell_status = true; 
      $("#logs").text("Cell Logs (on)"); 
     }
   else if(cell_status == true) 
     { 
      cell_status = false; 
      $("#logs").text("Cell Logs (off)"); 
     }
   else              
     { 
      cell_status = true; 
      $("#logs").text("Cell Logs (on)"); 
     }
   
   return cell_status;
  }

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
   if(max_value > value ) 
     { 
      value = (factor + 1) * interval; 
     } 
   if(Math.abs(max_value - value) <= interval_shift * interval) 
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
   if(min_value > value ) 
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
