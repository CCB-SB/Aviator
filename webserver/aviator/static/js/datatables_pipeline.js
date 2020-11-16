//
// Pipelining function for DataTables. To be used to the `ajax` option of DataTables
//
$.fn.dataTable.pipeline = function (opts) {
    // Configuration options
    var conf = $.extend({
        pages: 5,     // number of pages to cache
        url: '',      // script url
        data: null,   // function or object with parameters to send to the server
                      // matching how `ajax.data` works in DataTables
        method: 'GET' // Ajax HTTP method
    }, opts);

    // Private variables for storing the cache
    var cacheLower = -1;
    var cacheUpper = null;
    var cacheLastRequest = null;
    var cacheLastJson = null;

    return function (request, drawCallback, settings) {
        var ajax = false;
        var requestStart = request.start;
        var drawStart = request.start;
        var requestLength = request.length;
        var requestEnd = requestStart + requestLength;

        if (settings.clearCache) {
            // API requested that the cache be cleared
            ajax = true;
            settings.clearCache = false;
        }
        else if (cacheLower < 0 || requestStart < cacheLower || requestEnd > cacheUpper) {
            // outside cached data - need to make a request
            ajax = true;
        }
        else if (JSON.stringify(request.order) !== JSON.stringify(cacheLastRequest.order) ||
            JSON.stringify(request.columns) !== JSON.stringify(cacheLastRequest.columns) ||
            JSON.stringify(request.search) !== JSON.stringify(cacheLastRequest.search)
        ) {
            // properties changed (ordering, columns, searching)
            ajax = true;
        }

        // Store the request for checking next time around
        cacheLastRequest = $.extend(true, {}, request);

        if (ajax) {
            // Need data from the server
            if (requestStart < cacheLower) {
                requestStart = requestStart - (requestLength * (conf.pages - 1));

                if (requestStart < 0) {
                    requestStart = 0;
                }
            }

            cacheLower = requestStart;
            cacheUpper = requestStart + (requestLength * conf.pages);

            request.start = requestStart;
            request.length = requestLength * conf.pages;

            // Provide the same `data` options as DataTables.
            if ($.isFunction(conf.data)) {
                // As a function it is executed with the data object as an arg
                // for manipulation. If an object is returned, it is used as the
                // data object to submit
                var d = conf.data(request);
                if (d) {
                    $.extend(request, d);
                }
            }
            else if ($.isPlainObject(conf.data)) {
                // As an object, the data given extends the default
                $.extend(request, conf.data);
            }

            settings.jqXHR = $.ajax({
                "type": conf.method,
                "url": conf.url,
                "data": request,
                "dataType": "json",
                "cache": false,
                "success": function (json) {
                    cacheLastJson = $.extend(true, {}, json);

                    if (cacheLower != drawStart) {
                        json.data.splice(0, drawStart - cacheLower);
                    }
                    if (requestLength >= -1) {
                        json.data.splice(requestLength, json.data.length);
                    }

                    drawCallback(json);
                   updatePlot(json.website_states);
                }
            });
        }
        else {
            json = $.extend(true, {}, cacheLastJson);
            json.draw = request.draw; // Update the echo for each response
            json.data.splice(0, requestStart - cacheLower);
            json.data.splice(requestLength, json.data.length);

            drawCallback(json);
        }
    }
};

// Register an API method that will empty the pipelined data, forcing an Ajax
// fetch on the next draw (i.e. `table.clearPipeline().draw()`)
$.fn.dataTable.Api.register('clearPipeline()', function () {
    return this.iterator('table', function (settings) {
        settings.clearCache = true;
    });
});

      //last 30 days annotation
      var checked_days = 30;
      var days = [];
      for(var d = 0; d < checked_days; ++d) {
        //var currentTime = new Date(Date.now()).getDate();
        var currentTime = new Date();
        currentTime.setDate(currentTime.getDate()-d);
        days.push(currentTime.getFullYear() + "-" +  ((currentTime.getMonth() + 1) < 10 ? "0" : "") + (currentTime.getMonth() + 1) + "-" + (currentTime.getDate()  < 10 ? "0" : "") + currentTime.getDate());
      }
      function updatePlot(states) {
        if(table == null) {
          return;
        }

        //cdata = table.columns({filter:'applied'}).data()[1];
        //adata = table.columns({filter:'applied'}).data()[6];
        var zdata_online = [];
        var zdata_offline = [];
        var zdata_na = [];
        var tdata = [];
        var xdata = [];

        for (var i = 0; i < states.length; ++i) {
              var tmp_data_online = [];
              var tmp_data_offline = [];
              var tmp_data_na = [];
              var tmp_tdata = [];
              if (states[i]["websites__states"] == null) {
                  continue;
              }
              for (var j = 0; j < states[i]["websites__states"].length; ++j) {
                let value = states[i]["websites__states"][j];
                if (value == null) {
                  tmp_data_online.push(null)
                  tmp_data_na.push(1)
                  tmp_data_offline.push(null)
                }
                else if(value){
                  tmp_data_online.push(1)
                  tmp_data_na.push(null)
                  tmp_data_offline.push(null)
                } else {
                  tmp_data_online.push(null)
                  tmp_data_na.push(null)
                  tmp_data_offline.push(1)
                }
                tmp_tdata.push("PMID: "+states[i]["pubmed_id"]);
              }

            zdata_online.push(tmp_data_online);
            zdata_offline.push(tmp_data_offline);
            zdata_na.push(tmp_data_na);
            tdata.push(tmp_tdata);
        }

        var online_data = {
          z: zdata_online,
          //z: [[1, 1, null], [null, null, 1]],
          x: days,
          text: tdata,
          colorscale: [[0, '#1a9850'], [1, '#1a9850']],
          colorbar: {
            len: 0.3,
            y: 1,
            yanchor: 'top',
            tickvals: [0],
            ticktext: [''],
            title: 'Online'
          },
          type: 'heatmap'
        };

        var offline_data = {
          z: zdata_offline,
          //z: [[null, null, 1], [1, null, null]],
          x: days,
          text: tdata,
          colorscale: [[0, '#d73027'], [1, '#d73027']],
          colorbar: {
            len: 0.3,
            y: 0.5,
            yanchor: 'middle',
            tickvals: [0],
            ticktext: [''],
            title: 'Offline'
          },
          type: 'heatmap'
        };

        var na_data = {
          z: zdata_na,
          //z: [[null, 1, null], [1, 1, null]],
          x: days,
          text: tdata,
          colorscale: [[0, '#969696'], [1, '#969696']],
          colorbar: {
            len: 0.3,
            y: 0,
            yanchor: 'bottom',
            tickvals: [0],
            ticktext: [''],
            title: 'NA'
          },
          type: 'heatmap'
        };

        var data = [online_data, offline_data, na_data];

        var layout = {
          plot_bgcolor:"#454d55",
          paper_bgcolor:"#454d55",
          font: {
            color: '#dee2e6'
          }
        };
        Plotly.newPlot(document.getElementById('heatmap'), data, layout);
      }

