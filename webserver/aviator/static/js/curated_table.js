
var author_column = null;
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
          setTimeout(function(){
            apply_foreign_filter();
          }, 0);

      },
        ajax: $.fn.dataTable.pipeline({
            url: tbl_data_url,
            pages: 5 // number of pages to cache
        }),
      "columnDefs": [
        { 'visible': false, 'targets': [9, 11, 12, 13, 14] },
        { 'width': 50, 'targets': [1, 4, 13] },
        { 'width': 100, 'targets': [2, 5, 6, 7, 14, 15] },
        { 'width': 250, 'targets': [0] },
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
                  return "<div style='display:none'>" + data + "</div><span class='" + (data === "T" ? "orange" : (data === "O" ? "green" : (data === "F" ? "red" : "grey"))) + "-circle'></span>";
                }
            }, {
                "data": "percentage", render: function ( data ) {
                  return "<div style='display:none'>" + (data < 10 ? ("00" + data) : (data < 100 ? ("0" + data) : data)) + "</div>" + (data === null ? "No Data" : (data+"% Online"));;
                }
            }, {
                "data": "authors", render: function ( data ) {
				  str = "";
				  for(var n=0; n < data.length; n++) {
                      str += "<a href=\"#\" onclick=\"filterAuthor('"+data[n]+"');return false;\">"+data[n]+"</a>";
				      if(n < data.length - 1) {
				          str += ", ";
                      }
                  }
				  return str;
				}
            }, {
                data: "year",
            }, {
                data: "journal"
            }, {
                data: "pubmed_id"
            }, {
                "data": "description", render: function ( data ) {
				  ++abstract_counter;
				  return "<div style='display:none' id='abstract" + abstract_counter + "'>" + data + "</div><button class=\"btn btn-outline-light my-2 my-sm-0\" type=\"button\" data-toggle=\"modal\" onclick=\"abstract(document.getElementById('abstract" + abstract_counter + "').innerHTML)\" data-target=\"#abstractModal\">Abstract</button>";
				}
            }, {
                data: "url", render: function ( data ) {
				return "<a href=\""+data+"\">"+data+"</a>";
			  }
            },
          { "data": "tag_tags", render: function ( data ) {
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
          { data: "website", render: function ( data ) {
            return "<a class='btn btn-outline-light' href='details/"+data+"'>Show Details</a>"
          }
        }
        ],
        initComplete: function () {
            var api = this.api();
            api.columns().every(function (i) {
                var column = this;
                if (numeric_cols.indexOf(i) !== -1) {
                    var el = $('<input id="cs_from_'+i+'" class="form-control" placeholder="From" type="number"></th>')
                        .appendTo($(column.footer()).empty())
                        .on('change', function () {
                            column
                                .search($(this).val() + ";" + $(this).next().val())
                                .draw();
                        });
                    $('<input id="cs_to_'+i+'" class="form-control" placeholder="To" type="number"></th>')
                        .appendTo(el.parent())
                        .on('change', function () {
                            column
                                .search($(this).prev().val() + ";" + $(this).val())
                                .draw();
                        });
                    $( "#cs_from_"+i ).autocomplete({
                      source: function(request, response) {
                        $.getJSON(autocomplete_url, createTableSearchData(i), response);
                      },
                      minLength: 1
                    });
                    $( "#cs_to_"+i ).autocomplete({
                      source: function(request, response) {
                        $.getJSON(autocomplete_url, createTableSearchData(i), response);
                      },
                      minLength: 1
                    });
                } else {
                    $('<input id="cs_'+i+'" class="form-control" type="text" value="'+(search_column == i ? search_string : "")+'"></th>')
                        .appendTo($(column.footer()).empty())
                        .on('change', function () {
                            column
                                .search($(this).val())
                                .draw();
                        });
                    $( "#cs_"+i ).autocomplete({
                      source: function(request, response) {
                        $.getJSON(autocomplete_url, createTableSearchData(i), response);
                      },
                      minLength: 1
                    });
                    if ($( "#cs_"+i ).autocomplete().data("ui-autocomplete") != undefined) {
                        $( "#cs_"+i ).autocomplete().data("ui-autocomplete")._renderItem = function (ul, item) {
                            var newText = String(item.value).replace(
                            new RegExp(this.term, "gi"),
                            "<span style='color: white;'><strong>$&</strong></span>");
                            return $("<li></li>")
                            .data("item.autocomplete", item)
                            .append("<div>" + newText + "</div>")
                            .appendTo(ul);
                        };
                    }
                }
                if(search_column == i) {
                    scolumn = column;
                }
                if(i == 3) {
                    author_column = column;
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
function email(address) {
    window.location.href = "mailto:"+address.replace("[at]", "@");
}

function createTableSearchData(column) {
    data = {};
    for(var i=0; i <= 14; ++i) {
        if(i==4 && $('#cs_from_'+i).val() != undefined) {
            data[i] = $('#cs_from_'+i).val()+";"+$('#cs_to_'+i).val();
        } else if($('#cs_'+i).val() != undefined) {
            data[i] = $('#cs_'+i).val();
        }
    }
    data["q"] = column;
    return data;
}

function filterAuthor(author_name) {
    document.getElementById("cs_3").value = author_name;
  author_column.search(author_name).draw();
}