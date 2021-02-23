var author_column = null;
var tags_column = null;
var numeric_cols = [2, 4];
$(document).ready(function () {

    let cookie_name = "curated_hidden_columns";
    let standard_hidden_columns = [];
    let save_days = 60;

    async function apply_foreign_filter() {
        if (foreign_filter_exists && scolumn != null) {
            foreign_filter_exists = false;
            scolumn.search(search_string).draw();
        }
    }

    var scolumn = null;
    var foreign_filter_exists = true;
    var description_counter = 0;
    var collapse_counter = 0;

    function createCellText(text, without_p = false) {
        return '<div class="module"><a class="read-more collapsed" ' +
            'data-toggle="collapse" href="#collapse' + (++collapse_counter) + '" role="button"></a><div class="collapse" ' +
            'id="collapse' + (collapse_counter) + '" aria-expanded="false">' +
            (without_p ? '' : '<p>') + text + (without_p ? '' : '</p>') + '</div></div>';
    }

    var table = $('#table').DataTable({
        fnDrawCallback: function (settings) {
            setTimeout(function () {
                apply_foreign_filter();
            }, 0);

        },
        ajax: $.fn.dataTable.pipeline({
            url: tbl_data_url,
            pages: 5 // number of pages to cache
        }),
        "columnDefs": [
            {'visible': false, 'targets': getHiddenColumns()},
            {'width': 50, 'targets': [1, 4]},
            {'width': 100, 'targets': [2, 5, 6, 7]},
            {'width': 250, 'targets': [0]},
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
                "data": "title", render: function (data) {
                    return createCellText(data);
                }
            }, {
                "data": "status", render: function (data) {
                    return createCellText("<div style='display:none'>" + data + "</div><span class='" + (data === "TEMP_OFFLINE" ? "orange" : (data === "ONLINE" ? "green" : (data === "OFFLINE" ? "red" : "grey"))) + "-circle'></span>");
                }
            }, {
                "data": "percentage", render: function (data) {
                    return createCellText("<div style='display:none'>" + (data < 10 ? ("00" + data) : (data < 100 ? ("0" + data) : data)) + "</div>" + (data === null ? "No Data" : (parseFloat(data.toFixed(2)) + "% Online")), true);
                }
            }, {
                "data": "authors", render: function (data) {
                    str = "";
                    for (var n = 0; n < data.length; n++) {
                        str += "<a href=\"#\" onclick=\"filterAuthor('" + data[n] + "');return false;\">" + data[n] + "</a>";
                        if (n < data.length - 1) {
                            str += ", ";
                        }
                    }
                    return createCellText(str);
                }
            }, {
                "data": "year", render: function (data) {
                    return createCellText(data);
                }
            }, {
                "data": "journal", render: function (data) {
                    return createCellText(data);
                }
            }, {
                "data": "pubmed_id", render: function (data) {
                    return createCellText(`<a target="_blank" rel="noopener" href="https://pubmed.ncbi.nlm.nih.gov/${data}/">${data}</a>`);
                }
            }, {
                "data": "description", render: function (data) {
                    ++description_counter;
                    return createCellText("<div style='display:none' id='description" + description_counter + "'>" + data + "</div><button class=\"btn btn-outline-light my-2 my-sm-0\" type=\"button\" data-toggle=\"modal\" onclick=\"showDescription(document.getElementById('description" + description_counter + "').innerHTML)\" data-target=\"#descriptionModal\">Description</button>", true);
                }
            }, {
                data: "url", render: function (data) {
                    return createCellText("<a href=\"" + data + "\">" + data + "</a>");
                }
            }, {
                "data": "tag_tags", render: function (data) {
                    str = "";
                    for (var n = 0; n < data.length; n++) {
                        str += (n > 0 ? ", " : "") + "<a href=\"#\" onclick=\"filterTags('" + data[n] + "');return false;\">" + data[n] + "</a>";
                    }
                    return createCellText(str);
                }
            },
            {
                data: "website_pk", render: function (data) {
                    return createCellText("<a class='btn btn-outline-light' href='details/" + data[0] + "'>Show Details</a>")
                }
            }
        ],
        initComplete: function () {
            var api = this.api();
            api.columns().every(function (i) {
                var column = this;
                // websites
                if (i == 10) {
                    $('<p></p>').appendTo($(column.footer()).empty())
                    return;
                }
                // status
                if (i == 1) {
                    $('<select class="custom-select form-control">' +
                        '<option value></option>' +
                        '<option value="O">Online</option>' +
                        '<option value="F">Offline</option>' +
                        '<option value="T">Temp. offline</option>').appendTo($(column.footer()).empty())
                        .on('change', function () {
                            column
                                .search($(this).val())
                                .draw();
                        });
                    return;
                }
                if (numeric_cols.indexOf(i) !== -1) {
                    var el = $('<input id="cs_from_' + i + '" class="form-control" placeholder="From" type="number"></th>')
                        .appendTo($(column.footer()).empty())
                        .on('change', function () {
                            column
                                .search($(this).val() + ";" + $(this).next().val())
                                .draw();
                        });
                    $('<input id="cs_to_' + i + '" class="form-control" placeholder="To" type="number"></th>')
                        .appendTo(el.parent())
                        .on('change', function () {
                            column
                                .search($(this).prev().val() + ";" + $(this).val())
                                .draw();
                        });
                } else {
                    var el = $('<input id="cs_' + i + '" class="form-control" type="text" value="' + (search_column == i ? search_string : "") + '"></th>')
                        .appendTo($(column.footer()).empty());

                    el.autocomplete({
                        source: function (request, response) {
                            $.getJSON(autocomplete_url, createTableSearchData(i), response);
                        },
                        minLength: 1
                    });
                    el.on('change', function () {
                        column
                            .search($(this).val())
                            .draw();
                    });
                    if (el.autocomplete().data("ui-autocomplete") != undefined) {
                        el.autocomplete().data("ui-autocomplete")._renderItem = function (ul, item) {
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
                if (search_column == i) {
                    scolumn = column;
                }
                if (i == 3) {
                    author_column = column;
                }
                if (i == 9) {
                    tags_column = column;
                }
            });
            $('#table_wrapper .col-sm-12:eq(0)').html(api.buttons().container());
            $('#table_filter').remove();

            $(function () {
                $('[data-toggle="tooltip"]').tooltip()
            })
        },
        lengthMenu: [
            [5, 10, 25, 50],
            ['5 rows', '10 rows', '25 rows', '50 rows']
        ],
        buttons: [
            'pageLength', {text: '<i class="fas fa-download"></i> CSV', action: showExportCSVModal},
            {
                extend: 'colvis',
                action: function (e, dt, node, config) {
                    $.fn.dataTable.ext.buttons.collection.action.call(this, e, dt, node, config);
                },
                text: '<i class="fas fa-eye-slash"></i> Column visibility'
            }
        ]
    });

    var hidden_columns = getHiddenColumns();
    table.on('buttons-action', function (e, buttonApi, dataTable, node, config) {
        if (config.columns != undefined) {
            let index = config.columns;
            if (hidden_columns.includes(index)) {
                var sindex = hidden_columns.indexOf(index);
                if (sindex > -1) {
                    hidden_columns.splice(sindex, 1);
                }
            } else {
                hidden_columns.push(index);
            }
            setCookie(cookie_name, hidden_columns.toString(), save_days);
        }
    });

    function getHiddenColumns() {
        if (getCookie(cookie_name) == null) {
            return standard_hidden_columns;
        }
        values = getCookie(cookie_name).split(",");
        ret = [];
        for (let i = 0; i < values.length; i++) {
            ret.push(parseInt(values[i]));
        }
        return ret;
    }

    function setCookie(name, value, days) {
        var expires = "";
        if (days) {
            var date = new Date();
            date.setTime(date.getTime() + (days * 24 * 60 * 60 * 1000));
            expires = "; expires=" + date.toUTCString();
        }
        document.cookie = name + "=" + (value || "") + expires + "; path=/";
    }

    function getCookie(name) {
        var nameEQ = name + "=";
        var ca = document.cookie.split(';');
        for (var i = 0; i < ca.length; i++) {
            var c = ca[i];
            while (c.charAt(0) == ' ') c = c.substring(1, c.length);
            if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length, c.length);
        }
        return null;
    }
});

function email(address) {
    window.location.href = "mailto:" + address.replace("[at]", "@");
}

function createTableSearchData(column) {
    var data = {};
    for (var i = 0; i <= 14; ++i) {
        if (numeric_cols.indexOf(i) !== -1 && $('#cs_from_' + i).val() != undefined) {
            data[i] = $('#cs_from_' + i).val() + ";" + $('#cs_to_' + i).val();
        } else if ($('#cs_' + i).val() != undefined) {
            data[i] = $('#cs_' + i).val();
        }
    }
    data["q"] = column;
    return data;
}

function filterAuthor(author_name) {
    document.getElementById("cs_3").value = author_name;
    author_column.search(author_name).draw();
}

function filterTags(tag_name) {
    document.getElementById("cs_9").value = tag_name;
    tags_column.search(tag_name).draw();
}
