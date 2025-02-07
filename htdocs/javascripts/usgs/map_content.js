/**
 * Namespace: Map_Content
 *
 * Map_Content is a JavaScript library to provide instructions and general information.
 *
 * version 1.06
 * August 24, 2017
*/

/*
###############################################################################
# Copyright (c) 2017, Leonard Orzol Oregon Water Science Center
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

// Control for basemap dropdown
//
jQuery('#basemapMenu li').click(function(){

   // Dont do anything if the current active layer is clicked
   //
   if(!map.hasLayer(window[$(this).prop("id")])) 
     {			
       // Remove currently active basemap
       //
       jQuery("#basemapMenu li.active").each(function () 
          {
            // Console.log("removing ",$(this).attr("id"));
            //
            map.removeLayer(window[$(this).prop("id")]);
            jQuery(this).removeClass("active");
          });
			
       // Make new selection active and add to map
       //
       jQuery(this).addClass('active');
       map.addLayer(window[$(this).prop("id")]);
     }
});

// Monitor About click
//-----------------------------------------------
jQuery('.aboutButtons').click(function() {

    var div_id   = jQuery(this).prop("id");
    var file     = aboutFiles[div_id];
    var InfoText = loadText(file);

    jQuery("#aboutContent").html(InfoText);
   
});

// Monitor About click
//-----------------------------------------------
jQuery('#aboutTabs').click(function() {

    for (var div_id in aboutFiles)
      {
        var file     = aboutFiles[div_id];
        var InfoText = loadText(file);

        jQuery("#" + div_id).html(InfoText);
      }
   
});

// Add hover function to choices of data
//
jQuery('#dataDrop').hover(function() {

    var GW    = new Object;
    var SW    = new Object;
    var gw_re = /groundwater/;
    var sw_re = /surface water/;

    for(var category in wateruse_categories)
      {
          if(gw_re.test(category))
            {
              GW[category] = 1;
            };
          if(sw_re.test(category))
            {
              SW[category] = 1;
            };

      }

    var html = [];
    html.push(" <li>");
    html.push('  <a href="#" data-toggle="collapse" data-target=".navbar-collapse.in" onclick="fullExtent2(); return false;">');
    html.push('   &nbsp;&nbsp;Zoom To Full Extent');
    html.push('  </a>');
    html.push(" </li>");
   
    jQuery("#dataList").html(html.join(""));
});

// Add hover function to keep menus open
//
//jQuery('ul.nav li.dropdown').hover(function() 
//  {
//    jQuery(this).find('.dropdown-menu').stop(true, true).delay(100).fadeIn();
//  }, function() {
//    jQuery(this).find('.dropdown-menu').stop(true, true).delay(100).fadeOut();
//  });

// Monitor Legend hover
//-----------------------------------------------
jQuery('#legendDrop2').hover(function() {

    if(isLegend === true)
      {
         isLegend = false;
         var files = ["legend"];
         while(files.length > 0)
           {
             var file     = files.shift();
             var InfoText = loadText(file + ".txt");
             alert("Hover legend " + InfoText);
             var json     = JSON.stringify(eval("(" + InfoText + ")"));
             alert("Hover legend " + json["Boundary of active cells layer 1"].description);

             var legend_txt = [
                               '<li role="presentation">',
                               '<a role="menuitem" tabindex="-1">',
                               InfoText,
                               '</li>'
                              ].join(" ");
     
             jQuery("#overlayLegend").append(legend_txt);
             //jQuery("#overlayLegend").append(legend_txt);
             //jQuery("#overlayMenu.sw").append('<li role="presentation" id="' + curSiteTypeInfo.overlayLayerName + '" class="' + curSiteTypeInfo.overlayLayerName  + '"><a role="menuitem" tabindex="-1"><div name="overlayLayers" ><img src="' + curSiteTypeInfo.singleMarkerURL + '"/><span>' + curSiteTypeInfo.legendLayerName + '</span></div></li>');
           }
      }
   
});

function capitalizeEachWord(str) {
    return str.replace(/\w\S*/g, function(txt) {
        return txt.charAt(0).toUpperCase() + txt.substr(1).toLowerCase();
    });
}
