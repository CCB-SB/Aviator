$(document).ready(function () {
    var numeric_cols = [4];

    var render_format_num = {
        display: function (data, type, row, meta) {
            return table_display_helper.numberWithCommas(data);
        }
    };
    var abstract_counter = 0;
    var table = $('#table').DataTable({
      drawCallback: function( settings ) {
          //updatePlot();
      },
        ajax: $.fn.dataTable.pipeline({
            url:tbl_data_url,
            pages: 5 // number of pages to cache
        }),
      "columnDefs": [
        { 'visible': false, 'targets': [11] },
      ],
        deferRender: true,
        processing: true,
        serverSide: true,
        lengthChange: false,
        order: [[0, "asc"]],
        pagingType: "input",
        autoWidth: false,
        columns: [
            {
                data: "title"
            }, {
                "data": "status", render: function ( data ) {
          var str = "";
          if (data != null) {
            for (var i=0; i < data.length; i++) {
              if(data[i] != null) {
                state = data[i];
                if (i > 0) {
                  str += "&nbsp;";
                }
                str += "<div style='display:none'>" + (state == null ? 1 : (state ? 0 : 2)) + "</div><span class='" + (state == null ? "orange" : (state ? "green" : "red")) + "-circle'></span>";
              }
            }
          }
          return str;
        }
            }, {
                "data": "percentage", render: function ( data ) {
          var str = "";
          if (data != null) {
            for (var i=0; i < data.length; i++) {
              if(data[i] != null) {
                percentage = data[i];
                if (i > 0) {
                  str += " / ";
                }
                str += "<div style='display:none'>" + (percentage < 10 ? ("00" + percentage) : (percentage < 100 ? ("0" + percentage) : percentage)) + "</div>" + (percentage == -1 ? "No Data" : (percentage+"% Online"));
              }
            }
          }
          return str;
        }
            }, {
                data: "authors",
            }, {
                data: "year",
            }, {
                data: "journal"
            }, {
                data: "pubmed_id"
            }, {
                "data": "abstract", render: function ( data ) {
				  ++abstract_counter;
				  return "<div style='display:none' id='abstract" + abstract_counter + "'>" + data.replaceAll('""""""""', '"') + "</div><button class=\"btn btn-outline-light my-2 my-sm-0\" type=\"button\" data-toggle=\"modal\" onclick=\"abstract(document.getElementById('abstract" + abstract_counter + "').innerHTML)\" data-target=\"#abstractModal\">Abstract</button>";
				}
            }, {
                data: "original_url", render: function ( data ) {
				var str = "";
				if (data != null) {
				  var i = 0;
				  for (i=0; i < data.length; i++) {
					if(data[i] != null) {
					  str += (i > 0 ? ", " : "") + data[i];
					}
				  }
				}
				return str;
			  }
            }, { data: "derived_url", render: function ( data ) {
        var str = "";
        if (data != null) {
          var i = 0;
          for (i=0; i < data.length; i++) {
            if(data[i] != null) {
              str += (i > 0 ? ", " : "") + data[i];
            }
          }
        }
        return str;
      }
      },
      { "data": "contact_mail"},
      { "data": "user_kwds", render: function ( data ) {
        var str = "";
        if (data != null) {
          var i = 0;
          for (i=0; i < data.length; i++) {
            if(data[i] != null) {
              str += (i > 0 ? ", " : "") + data[i];
            }
          }
        }
        return str;
      }
      },
      { data: "website_pks", render: function ( data ) {
        var str = "";
        if (data != null) {
          var i = 0;
          for (i=0; i < data.length; i++) {
            if(data[i] != null) {
              str += "<a class='btn btn-outline-light' href='details/" + data[i] + "'>Website</a>"
            }
          }
          if(i == 0) {
            return "<a class='btn btn-outline-light' href='details/"+data+"'>Website</a>"
          }
        }
        return str;
      }
    }
        ],
        initComplete: function () {
            var api = this.api();
            api.columns().every(function (i) {
                var column = this;
                if (numeric_cols.indexOf(i) !== -1) {
                    var el = $('<input class="form-control" placeholder="From" type="number"></th>')
                        .appendTo($(column.footer()).empty())
                        .on('change', function () {
                            column
                                .search($(this).val() + ";" + $(this).next().val())
                                .draw();
                        });
                    $('<input class="form-control" placeholder="To" type="number"></th>')
                        .appendTo(el.parent())
                        .on('change', function () {
                            column
                                .search($(this).prev().val() + ";" + $(this).val())
                                .draw();
                        });
                } else {
                    $('<input class="form-control" type="text"></th>')
                        .appendTo($(column.footer()).empty())
                        .on('change', function () {
                            column
                                .search($(this).val())
                                .draw();
                        });
                }
            });

            $('#table_wrapper .col-sm-12:eq(0)').html(api.buttons().container());
            $('#table_filter').remove();

            $(function () {
                $('[data-toggle="tooltip"]').tooltip()
            })
        },
        lengthMenu: [
            [10, 25, 50],
            ['10 rows', '25 rows', '50 rows']
        ],
        buttons: [
            'pageLength', 'excel', 'csv', 'colvis'
        ]
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
      function updatePlot() {
        if(table == null) {
          return;
        }
        cdata = table.columns({filter:'applied'}).data()[1];
        adata = table.columns({filter:'applied'}).data()[6];
        var zdata_online = [];
        var zdata_offline = [];
        var zdata_na = [];
        var tdata = [];
        var xdata = [];


        for (var k = 0; k < cdata.length; ++k) {
          for (var i = 0; i < cdata[k].length; ++i) {
            if(cdata[k][i] in websites) {
              var tmp_data_online = [];
              var tmp_data_offline = [];
              var tmp_data_na = [];
              var tmp_tdata = [];
              for (var j = 0; j < websites[cdata[k][i]]["states"].length; ++j) {
                let value = websites[cdata[k][i]]["states"][j];
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
                tmp_tdata.push("PMID: "+adata[i]);
              }
            }

            zdata_online.push(tmp_data_online);
            zdata_offline.push(tmp_data_offline);
            zdata_na.push(tmp_data_na);
            tdata.push(tmp_tdata);
          }
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

})
;
