/**
 * Namespace: Loading
 *
 * Loading is a JavaScript library to make a Ajax request and show a 
 *  loading modal screen.
 *
 * version 1.07
 * October 16, 2016
*/

/*
###############################################################################
# Copyright (c) 2016 Leonard Orzol <llorzol@usgs.gov>
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

function webRequest(request_type, script_http, data_http, dataType, callFunction)
  {
    ajaxindicatorstart('loading data.. please wait..');  

    var myData   = null;

    $.support.cors = true;

    myPromise = $.ajax({
                        type: request_type,
                        url: script_http,
                        data: data_http,
                        dataType: dataType
                      })
        .done(function(myData) { 
                                ajaxindicatorstop();
                                callFunction(myData); 
                               })
        .fail(function(jqXHR, textStatus) {
             ajaxindicatorstop();

             var message = "";

             if(jqXHR.status === 0)
               {
                 message = "Not connect.n Verify Network.";
               }
            else if (jqXHR.status == 404)
               {
                 message = "Requested page not found. [404]";
               }
            else if (jqXHR.status == 500)
               {
                 message = "Internal Server Error [500].";
               }
            else if (exception === 'parsererror')
               {
                 message = "Requested JSON parse failed.";
               }
            else if (exception === 'timeout')
               {
                 message = "Time out error.";
               }
            else if (exception === 'abort')
               {
                 message = "Ajax request aborted.";
               }
            else
               {
                 message = "Uncaught Error.n";
               }

            message    += "<br> please wait while the page is refreshed";
            openModal(message);
            return false;
           });

  }

function webRequest2(request_type, script_http, data_http, dataType, callFunction)
  {
    var myData = null;

    ajaxindicatorstart('loading data.. please wait..');  

    $.support.cors = true;

    $.ajax({
            type: request_type,
            url: script_http,
            data: data_http,
            dataType: dataType
           }) 
        .done(function(myData)
           {
             ajaxindicatorstop();
             callFunction(myData);
           })
        .fail(function(jqXHR, textStatus)
           {
             ajaxindicatorstop();

             var message = "";

             if(jqXHR.status === 0)
               {
                 message = "Not connect.n Verify Network.";
               }
            else if (jqXHR.status == 404)
               {
                 message = "Requested page not found. [404]";
               }
            else if (jqXHR.status == 500)
               {
                 message = "Internal Server Error [500].";
               }
            else if (exception === 'parsererror')
               {
                 message = "Requested JSON parse failed.";
               }
            else if (exception === 'timeout')
               {
                 message = "Time out error.";
               }
            else if (exception === 'abort')
               {
                 message = "Ajax request aborted.";
               }
            else
               {
                 message = "Uncaught Error.n";
               }

            message    += "<br> please wait while the page is refreshed";
            openModal(message);
            return false;
           });
  }

function webRequest3(request_type, script_http, data_http, dataType, callBack)
  {

      //var myData = null;

    $.support.cors = true;

    $.when(

    $.ajax({
            type: request_type,
            url: script_http,
            data: data_http,
            dataType: dataType,
            beforeSend: function() { ajaxindicatorstart('loading data.. please wait..'); }
           }) 
        .done(function(myData)
           {
              ajaxindicatorstop();
              //return myData;
           })
        .fail(function(jqXHR, textStatus)
           {
             var message = "";

             if(jqXHR.status === 0)
               {
                 message = "Not connect.n Verify Network.";
               }
            else if (jqXHR.status == 404)
               {
                 message = "Requested page not found. [404]";
               }
            else if (jqXHR.status == 500)
               {
                 message = "Internal Server Error [500].";
               }
            else if (exception === 'parsererror')
               {
                 message = "Requested JSON parse failed.";
               }
            else if (exception === 'timeout')
               {
                 message = "Time out error.";
               }
            else if (exception === 'abort')
               {
                 message = "Ajax request aborted.";
               }
            else
               {
                 message = "Uncaught Error.n";
               }

            message    += "<br> please wait while the page is refreshed";
            openModal(message);
            return false;
           }))

        .then( function( myData, textStatus, jqXHR ) { return myData; } ); 

  }

function ajaxindicatorstart(text)
  {
    if(jQuery('body').find('#resultLoading').prop('id') != 'resultLoading')
      {
        jQuery('body').append('<div id="resultLoading" style="display:none"><div><img src="images/ajax-loader.gif"><div>'+text+'</div></div><div class="bg"></div></div>');
      }

    jQuery('#resultLoading').css({
                                  'width':'100%',
                                  'height':'100%',
                                  'position':'fixed',
                                  'z-index':'10000000',
                                  'top':'0',
                                  'left':'0',
                                  'right':'0',
                                  'bottom':'0',
                                  'margin':'auto'
    });
  
    jQuery('#resultLoading .bg').css({
                                  'background':'#000000',
                                  'opacity':'0.7',
                                  'width':'100%',
                                  'height':'100%',
                                  'position':'absolute',
                                  'top':'0'
    });
  
    jQuery('#resultLoading>div:first').css({
                                  'width': '250px',
                                  'height':'75px',
                                  'text-align': 'center',
                                  'position': 'fixed',
                                  'top':'0',
                                  'left':'0',
                                  'right':'0',
                                  'bottom':'0',
                                  'margin':'auto',
                                  'font-size':'16px',
                                  'z-index':'10',
                                  'color':'#ffffff'
  
    });
  
    jQuery('#resultLoading .bg').height('100%');
    jQuery('#resultLoading').fadeIn(300);
    jQuery('body').css('cursor', 'wait');
  }

function ajaxindicatorstop()
  {
    jQuery('#resultLoading .bg').height('100%');
    jQuery('#resultLoading').fadeOut(300);
    jQuery('body').css('cursor', 'default');
  }
