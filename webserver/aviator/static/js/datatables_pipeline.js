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
                    updateStats(json.statistics);
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

function updatePlot(website_states) {
    if(table == null) {
      return;
    }

    var states = website_states.states;
    var dates = website_states.dates;
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
          wkey = "websites__states"
          if("states" in states[i]) {
              wkey = "states";
          }
          if (states[i][wkey] == null) {
              continue;
          }
          for (var j = 0; j < states[i][wkey].length; ++j) {
            let value = states[i][wkey][j];
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
      x: dates,
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
      x: dates,
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
      x: dates,
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
      },
      hovermode: 'closest',
    };
    Plotly.newPlot(document.getElementById('heatmap'), data, layout);
}

Chart.defaults.global.defaultFontColor = "#dee2e6";
var ctx1 = document.getElementById('overview').getContext('2d');
var ov_chart = new Chart(ctx1, {
    type: 'bar',
    data: {
        labels: ['Publications', 'Websites'],
        datasets: [{
            label: 'Offline',
            data: [0, 0],
            backgroundColor: [
                'rgba(255, 99, 132, 0.2)',
                'rgba(255, 99, 132, 0.2)'
            ],
            borderColor: [
                'rgba(255, 99, 132, 1)',
                'rgba(255, 99, 132, 1)'
            ],
            borderWidth: 2
        },{
            label: 'Temporarily offline',
            data: [0, 0],
            backgroundColor: [
              'rgba(255, 206, 86, 0.2)',
              'rgba(255, 206, 86, 0.2)',
            ],
            borderColor: [
                'rgba(255, 206, 86, 1)',
                'rgba(255, 206, 86, 1)',
            ],
            borderWidth: 2
        },{
            label: 'Online',
            data: [0, 0],
            backgroundColor: [
                'rgba(75, 192, 192, 0.2)',
                'rgba(75, 192, 192, 0.2)',
            ],
            borderColor: [
                'rgba(75, 192, 192, 1)',
                'rgba(75, 192, 192, 1)',
            ],
            borderWidth: 2
        },{
            label: 'Publications',
            data: [0, 0],
            backgroundColor: [
                'rgba(54, 162, 235, 0.2)',
                'rgba(255, 206, 86, 0.2)',

            ],
            borderColor: [
              'rgba(54, 162, 235, 1)',
              'rgba(255, 206, 86, 1)',
            ],
            borderWidth: 2
        }]
    },
    options: {
    title: {
        display: true,
        text: 'Database'
    },
        scales: {
            yAxes: [{
                ticks: {
                    beginAtZero: true
                },
              stacked: true,
            }],
          xAxes: [{
                stacked: true
          }],
        }
    }
});

var ctx2 = document.getElementById('c1').getContext('2d');
var temporal_chart = new Chart(ctx2, {
    type: 'line',
    data: {
        labels: [],
        datasets: [{
            label: 'Online',
            data: [],
            backgroundColor: [
                'rgba(75, 192, 192, 0.2)'
            ],
            borderColor: [
                'rgba(75, 192, 192, 1)'
            ],
            borderWidth: 2
        }, {
            label: 'Temporarily offline',
            data: [],
            backgroundColor: [
                'rgba(255, 206, 86, 0.2)',
            ],
            borderColor: [
                'rgba(255, 206, 86, 1)',
            ],
            borderWidth: 2
        }, {
            label: 'Offline',
            data: [],
            backgroundColor: [
                'rgba(255, 99, 132, 0.2)'
            ],
            borderColor: [
                'rgba(255, 99, 132, 1)'
            ],
            borderWidth: 2
        }]
    },
    options: {
        title: {
            display: true,
            text: 'Online/Offline'
        },
        scales: {
            yAxes: [{
                ticks: {
                    beginAtZero: true
                }
            }]
        }
    }
});

var pubs_per_year = new Chart(document.getElementById('pubs_per_year_availability').getContext('2d'), {
        type: 'bar',
        data: {
            labels: [],
            datasets: [{
                label: 'Online',
                data: [],
                backgroundColor: 'rgba(75, 192, 192, 0.2)',
                borderColor: 'rgba(75, 192, 192, 1)',
                borderWidth: 2
            }, {
                label: 'Temporarily offline',
                data: [],
                backgroundColor: 'rgba(255, 206, 86, 0.2)',
                borderColor: 'rgba(255, 206, 86, 1)',
                borderWidth: 2
            }, {
                label: 'Offline',
                data: [],
                backgroundColor: 'rgba(255, 99, 132, 0.2)',
                borderColor: 'rgba(255, 99, 132, 1)',
                borderWidth: 2
            }]
        },
        options: {
            title: {
                display: true,
                text: 'Availability per publication year'
            },
            scales: {
                yAxes: [{
                    ticks: {
                        beginAtZero: true
                    },
                  stacked: true,
                }],
              xAxes: [{
                                  stacked: true
              }],
            }
        }
    });

var top10_journals_chart = new Chart(document.getElementById('top10_journals_availability').getContext('2d'), {
        type: 'bar',
        data: {
            labels: [],
            datasets: [{
                label: 'Online',
                data: [],
                backgroundColor: 'rgba(75, 192, 192, 0.2)',
                borderColor: 'rgba(75, 192, 192, 1)',
                borderWidth: 2
            }, {
                label: 'Temporarily offline',
                data: [],
                backgroundColor: 'rgba(255, 206, 86, 0.2)',
                borderColor: 'rgba(255, 206, 86, 1)',
                borderWidth: 2
            }, {
                label: 'Offline',
                data: [],
                backgroundColor: 'rgba(255, 99, 132, 0.2)',
                borderColor: 'rgba(255, 99, 132, 1)',
                borderWidth: 2
            }]
        },
        options: {
            title: {
                display: true,
                text: 'Availability in 10 most popular journals'
            },
            scales: {
                yAxes: [{
                    ticks: {
                        beginAtZero: true
                    },
                  stacked: true,
                }],
              xAxes: [{
                    stacked: true
              }],
            }
        }
    });

function updateStats(stats){
    ov_chart.data.datasets[0].data[1] = stats.offline_count;
    ov_chart.data.datasets[1].data[1] = stats.temp_offline_count;
    ov_chart.data.datasets[2].data[1] = stats.online_count;
    ov_chart.data.datasets[3].data[0] = stats.paper_count;
    ov_chart.update();

    temporal_chart.data.labels = JSON.parse(stats.stat1_names)
    temporal_chart.data.datasets[0].data = JSON.parse(stats.stat1_online)
    temporal_chart.data.datasets[1].data = JSON.parse(stats.stat1_tmp_offline)
    temporal_chart.data.datasets[2].data = JSON.parse(stats.stat1_offline)
    temporal_chart.update()

    pubs_per_year.data.labels = JSON.parse(stats.pubs_per_year_names);
    pubs_per_year.data.datasets[0].data = JSON.parse(stats.pubs_per_year_online)
    pubs_per_year.data.datasets[1].data = JSON.parse(stats.pubs_per_year_tmp_offline)
    pubs_per_year.data.datasets[2].data = JSON.parse(stats.pubs_per_year_offline)
    pubs_per_year.update()

    top10_journals_chart.data.labels = JSON.parse(stats.top10_journals_names);
    top10_journals_chart.data.datasets[0].data = JSON.parse(stats.top10_journals_online)
    top10_journals_chart.data.datasets[1].data = JSON.parse(stats.top10_journals_tmp_offline)
    top10_journals_chart.data.datasets[2].data = JSON.parse(stats.top10_journals_offline)
    top10_journals_chart.update()
}
