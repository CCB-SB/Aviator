(function (table_display_helper, undefined) {
    table_display_helper.numberWithCommas = function(x) {
        if (typeof x == "number") {
            if(x == 0) return "0";
            if(Math.abs(x) < 0.001) return x.toExponential(2);
            var parts = parseFloat(x).toString().split(".");
            parts[0] = parts[0].replace(/\B(?=(\d{3})+(?!\d))/g, ",");
            return parts.join(".");
        } else if (x === null) {
            return "NA";
        } else {
            return x;
        }
    };

    table_display_helper.addNumericalRenderingForColumns = function(tbl) {
        for (var j = 0; j < tbl.columns.length; j++) {
            (function (i) {
                if (tbl.types[i]) {
                    tbl.columns[i] = {
                        title: tbl.columns[i].title,
                        render: {
                            "display": function (data, type, row, meta) {
                                return table_display_helper.numberWithCommas(data);
                            }
                        }
                    }
                }
            })(j);
        }
    }
}(window.table_display_helper = window.table_display_helper || {}));
