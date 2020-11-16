$(document).ready(function () {
    var numeric_cols = [4];

    var render_format_num = {
        display: function (data, type, row, meta) {
            return table_display_helper.numberWithCommas(data);
        }
    };
    async function apply_foreign_filter() {
        if (foreign_filter_exists && scolumn != null) {
          foreign_filter_exists = false;
          scolumn.search(search_string).draw();
        }
      }
    var scolumn = null;
     var foreign_filter_exists = true;
    var abstract_counter = 0;
    var table = $('#table').DataTable({
      fnDrawCallback: function( settings ) {
          //updatePlot();
          //TODO Fix
          setTimeout(function(){
            apply_foreign_filter();
          }, 0);

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
                    $('<input class="form-control" type="text" value="'+(search_column == i ? search_string : "")+'"></th>')
                        .appendTo($(column.footer()).empty())
                        .on('change', function () {
                            column
                                .search($(this).val())
                                .draw();
                        });
                }
                if(search_column == i) {
                    scolumn = column;
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

});
