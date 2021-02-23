var author_column = null;
var numeric_cols = [2, 4];

$(document).ready(function () {
    let cookie_name = "pubs_hidden_columns";
    let standard_hidden_columns = [9, 11, 12, 13, 14];
    let save_days = 60;

    async function apply_foreign_filter() {
        if (foreign_filter_exists && scolumn != null) {
            foreign_filter_exists = false;
            scolumn.search(search_string).draw();
        }
    }

    var scolumn = null;
    var foreign_filter_exists = true;
    var abstract_counter = 0;
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
            {'width': 50, 'targets': [1, 4, 13]},
            {'width': 100, 'targets': [2, 5, 6, 7, 14, 15]},
            {'width': 250, 'targets': [0]},
        ],
        deferRender: true,
        processing: true,
        serverSide: true,
        lengthChange: false,
        order: [[0, "asc"]],
        pagingType: "input",
        autoWidth: false,
        columns: [{
            "data": "title", render: function (data) {
                return createCellText(data);
            }
        }, {
            "data": "status", render: function (data) {
                var str = "";
                if (data != null) {
                    for (var i = 0; i < data.length; i++) {
                        if (data[i] != null) {
                            state = data[i];
                            if (i > 0) {
                                str += "&nbsp;";
                            }
                            str += "<div style='display:none'>" + state + "</div><span class='" + (state === "T" ? "orange" : (state === "O" ? "green" : (state === "F" ? "red" : "grey"))) + "-circle'></span>";
                        }
                    }
                }
                return createCellText(str);
            }
        }, {
            "data": "percentage", render: function (data) {
                var str = "";
                if (data != null) {
                    for (var i = 0; i < data.length; i++) {
                        if (data[i] != null) {
                            percentage = data[i];
                            if (i > 0) {
                                str += " / ";
                            }
                            str += "<div style='display:none'>" + (percentage < 10 ? ("00" + percentage) : (percentage < 100 ? ("0" + percentage) : percentage)) + "</div>" + (percentage === null ? "No Data" : (parseFloat(percentage.toFixed(2)) + "% Online"));
                        }
                    }
                }
                return createCellText(str, true);
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
            "data": "abstract", render: function (data) {
                ++abstract_counter;
                return createCellText("<div style='display:none' id='abstract" + abstract_counter + "'>" + data.replaceAll('""""""""', '"') + "</div><button class=\"btn btn-outline-light my-2 my-sm-0\" type=\"button\" data-toggle=\"modal\" onclick=\"abstract(document.getElementById('abstract" + abstract_counter + "').innerHTML)\" data-target=\"#abstractModal\">Abstract</button>", true);
            }
        }, {
            data: "original_url", render: function (data) {
                var str = "";
                if (data != null) {
                    var i = 0;
                    for (i = 0; i < data.length; i++) {
                        var url = data[i];
                        if (url != null) {
                            if (!url.startsWith("http")) {
                                url = `http://${url}`
                            }
                            str += (i > 0 ? ", " : "") + `<a target="_blank" rel="noopener" href="${url}">${data[i]}</a>`;
                        }
                    }
                }
                return createCellText(str);
            }
        }, {
            data: "derived_url", render: function (data) {
                var str = "";
                if (data != null) {
                    var i = 0;
                    for (i = 0; i < data.length; i++) {
                        if (data[i] != null) {
                            str += (i > 0 ? ", " : "") + `<a target="_blank" rel="noopener" href="${data[i]}">${data[i]}</a>`;
                        }
                    }
                }
                return createCellText(str);
            }
        },
            {
                data: "contact_mail", render: function (data) {
                    var str = "";
                    if (data != null) {
                        var i = 0;
                        for (i = 0; i < data.length; i++) {
                            if (data[i] != null) {
                                str += "<a href=\"#\" onclick=\"email(this.innerHTML);return false;\">" + data[i] + "</a>";
                            }
                        }
                        if (i == 0) {
                            return createCellText("<a href=\"#\" onclick=\"email(this.innerHTML);return false;\">" + data + "</a>");
                        }
                    }
                    return createCellText(str);
                }
            },
            {
                "data": "user_kwds", render: function (data) {
                    var str = "";
                    if (data != null) {
                        var i = 0;
                        for (i = 0; i < data.length; i++) {
                            if (data[i] != null) {
                                str += (i > 0 ? ", " : "") + data[i];
                            }
                        }
                    }
                    return createCellText(str);
                }
            },
            {
                "data": "scripts", render: function (data) {
                    var str = "";
                    if (data != null) {
                        var i = 0;
                        for (i = 0; i < data.length; i++) {
                            if (data[i] != null) {
                                str += (i > 0 ? ", " : "") + data[i];
                            }
                        }
                    }
                    return createCellText(str);
                }
            },
            {
                "data": "ssl", render: function (data) {
                    var str = "";
                    if (data != null) {
                        var i = 0;
                        for (i = 0; i < data.length; i++) {
                            if (data[i] != null) {
                                str += (i > 0 ? ", " : "") + (data[i] == null ? "NA" : (data[i] ? "Yes" : "No"));
                            }
                        }
                    }
                    return createCellText(str);
                }
            },
            {
                "data": "heap_size", render: function (data) {
                    var str = "";
                    if (data != null) {
                        var i = 0;
                        for (i = 0; i < data.length; i++) {
                            if (data[i] != null) {
                                str += (i > 0 ? ", " : "") + (parseInt(data[i]) < 1000000 ? "< 1 mb" : (parseInt(parseInt(data[i]) / 1000000) + " mb"));
                            }
                        }
                    } else {
                        str = "NA";
                    }
                    return createCellText(str);
                }
            },
            {
                data: "website_pks", render: function (data) {
                    var str = "";
                    if (data != null) {
                        var i = 0;
                        for (i = 0; i < data.length; i++) {
                            if (data[i] != null) {
                                str += "<a class='btn btn-outline-light' href='details/" + data[i] + "'>Show Details</a>"
                            }
                        }
                        if (i == 0) {
                            return createCellText("<a class='btn btn-outline-light' href='details/" + data + "'>Show Details</a>")
                        }
                    }
                    return createCellText(str);
                }
            }
        ],
        initComplete: function () {
            var api = this.api();
            api.columns().every(function (i) {
                var column = this;
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
                    // websites
                    if(i == 15){
                        $('<p></p>').appendTo($(column.footer()).empty())
                        return;
                    }
                    // status
                    if(i == 1){
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
                    // SSL
                    if(i == 13){
                        $('<select class="custom-select form-control">' +
                            '<option value></option>' +
                            '<option value="true">Yes</option>' +
                            '<option value="false">No</option>').appendTo($(column.footer()).empty())
                        .on('change', function () {
                            column
                                .search($(this).val())
                                .draw();
                        });
                        return;
                    }
                    var el = $('<input id="cs_' + i + '" class="form-control" type="text" value="' + (search_column == i ? search_string : "") + '"></th>')
                        .appendTo($(column.footer()).empty());

                    el.autocomplete({
                        source: function (request, response) {
                            $.getJSON(autocomplete_url, createTableSearchData(i), response);
                        },
                        minLength: 1,
                        change: function () {
                            column
                                .search($(this).val())
                                .draw();
                        }
                    });

                    el.on('keypress', function(e){
                        if(e.key == "Enter") {
                            column
                                .search($(this).val())
                                .draw();
                            el.autocomplete("close");
                        }
                    })

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
            console.log("getCookie: " + standard_hidden_columns);
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

function createTableSearchData(column) {
    data = {};
    for (var i = 0; i <= 15; ++i) {
        if (numeric_cols.indexOf(i) !== -1 && $('#cs_from_' + i).val() != undefined) {
            data[i] = $('#cs_from_' + i).val() + ";" + $('#cs_to_' + i).val();
        } else if ($('#cs_' + i).val() != undefined) {
            data[i] = $('#cs_' + i).val();
        }
    }
    data["q"] = column;
    return data;
}

function email(address) {
    window.location.href = "mailto:" + address.replace("[at]", "@");
}

function filterAuthor(author_name) {
    document.getElementById("cs_3").value = author_name;
    author_column.search(author_name).draw();
}
