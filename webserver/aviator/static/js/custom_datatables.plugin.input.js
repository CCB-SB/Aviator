(function ($) {
    function calcDisableClasses(oSettings) {
        var start = oSettings._iDisplayStart;
        var length = oSettings._iDisplayLength;
        var visibleRecords = oSettings.fnRecordsDisplay();
        var all = length === -1;

        var page = all ? 0 : Math.ceil(start / length);
        var pages = all ? 1 : Math.ceil(visibleRecords / length);

        var disableFirstPrevClass = (page > 0 ? '' : oSettings.oClasses.sPageButtonDisabled);
        var disableNextLastClass = (page < pages - 1 ? '' : oSettings.oClasses.sPageButtonDisabled);

        return {
            'first': disableFirstPrevClass,
            'previous': disableFirstPrevClass,
            'next': disableNextLastClass,
            'last': disableNextLastClass
        };
    }

    function calcCurrentPage(oSettings) {
        return Math.ceil(oSettings._iDisplayStart / oSettings._iDisplayLength) + 1;
    }

    function calcPages(oSettings) {
        return Math.ceil(oSettings.fnRecordsDisplay() / oSettings._iDisplayLength);
    }

    var firstClassName = 'first';
    var previousClassName = 'previous';
    var nextClassName = 'next';
    var lastClassName = 'last';

    $.fn.dataTableExt.oPagination.input = {
        "fnInit": function (oSettings, nPaging, fnCallbackDraw) {

            var pagList = document.createElement("ul");
            var nFirst = document.createElement('li');
            var nPrevious = document.createElement('li');
            var nNext = document.createElement('li');
            var nLast = document.createElement('li');
            var nInputli = document.createElement('li');
            var oLang = oSettings.oLanguage.oPaginate;

            nFirst.innerHTML = '<a href="#" class="page-link">&larr; ' + oLang.sFirst + '</a>';
            nPrevious.innerHTML = '<a href="#" class="page-link">&larr; ' + oLang.sPrevious + '</a>';
            nNext.innerHTML = '<a href="#" class="page-link">' + oLang.sNext + ' &rarr; </a>';
            nLast.innerHTML = '<a href="#" class="page-link">' + oLang.sLast + ' &rarr; </a>';

            nFirst.className = "paginate_button page-item first disabled";
            nPrevious.className = "paginate_button page-item previous disabled";
            nNext.className = "paginate_button page-item next";
            nLast.className = "paginate_button page-item last";

            nInputli.className = "page-item";
            nInputli.innerHTML = "<span style='padding: 3px 12px;' class='page-link'>&nbsp;&nbsp;&nbsp;Page <input type='text' class='form-control pagination-input'> of 0</span>"

            if (oSettings.sTableId !== '') {
                nPaging.setAttribute('id', oSettings.sTableId + '_paginate');
                nFirst.setAttribute('id', oSettings.sTableId + '_first');
                nPrevious.setAttribute('id', oSettings.sTableId + '_previous');
                nNext.setAttribute('id', oSettings.sTableId + '_next');
                nLast.setAttribute('id', oSettings.sTableId + '_last');
                nInputli.setAttribute('id', oSettings.sTableId + '_inputli');
                pagList.setAttribute("id", oSettings.sTableId + "_pagList");
            }
            pagList.appendChild(nFirst);
            pagList.appendChild(nPrevious);
            pagList.appendChild(nInputli);
            pagList.appendChild(nNext);
            pagList.appendChild(nLast);
            nPaging.appendChild(pagList);

            $(pagList).addClass('pagination');

            $(nFirst).click(function (e) {
                e.preventDefault();
                var iCurrentPage = calcCurrentPage(oSettings);
                if (iCurrentPage !== 1) {
                    oSettings.oApi._fnPageChange(oSettings, 'first');
                    fnCallbackDraw(oSettings);
                }
            });

            $(nPrevious).click(function (e) {
                e.preventDefault();
                var iCurrentPage = calcCurrentPage(oSettings);
                if (iCurrentPage !== 1) {
                    oSettings.oApi._fnPageChange(oSettings, 'previous');
                    fnCallbackDraw(oSettings);
                }
            });

            $(nNext).click(function (e) {
                e.preventDefault();
                var iCurrentPage = calcCurrentPage(oSettings);
                if (iCurrentPage !== calcPages(oSettings)) {
                    oSettings.oApi._fnPageChange(oSettings, 'next');
                    fnCallbackDraw(oSettings);
                }
            });

            $(nLast).click(function (e) {
                e.preventDefault();
                var iCurrentPage = calcCurrentPage(oSettings);
                if (iCurrentPage !== calcPages(oSettings)) {
                    oSettings.oApi._fnPageChange(oSettings, 'last');
                    fnCallbackDraw(oSettings);
                }
            });

            $(nInputli).find('input').keyup(function (e) {
                // 38 = up arrow, 39 = right arrow
                if (e.which === 38 || e.which === 39) {
                    this.value++;
                }
                // 37 = left arrow, 40 = down arrow
                else if ((e.which === 37 || e.which === 40) && this.value > 1) {
                    this.value--;
                }

                if (this.value === '' || this.value.match(/[^0-9]/)) {
                    /* Nothing entered or non-numeric character */
                    this.value = this.value.replace(/[^\d]/g, ''); // don't even allow anything but digits
                    return;
                }

                var iNewStart = oSettings._iDisplayLength * (this.value - 1);
                if (iNewStart < 0) {
                    iNewStart = 0;
                }
                if (iNewStart >= oSettings.fnRecordsDisplay()) {
                    iNewStart = (Math.ceil((oSettings.fnRecordsDisplay() - 1) / oSettings._iDisplayLength) - 1) * oSettings._iDisplayLength;
                }

                oSettings._iDisplayStart = iNewStart;
                fnCallbackDraw(oSettings);
            });

            /* Take the brutal approach to cancelling text selection */
//        $('span', nPaging).bind('mousedown', function () { return false; });
//        $('span', nPaging).bind('selectstart', function () { return false; });

            // If we can't page anyway, might as well not show it
//        var iPages = Math.ceil((oSettings.fnRecordsDisplay()) / oSettings._iDisplayLength);
//        if (iPages <= 1) {
//            $(nPaging).hide();
//        }
        },
        "fnUpdate": function (oSettings, fnCallbackDraw) {
            if (!oSettings.aanFeatures.p) {
                return;
            }
            var iPages = calcPages(oSettings);
            var iCurrentPage = calcCurrentPage(oSettings);

            var an = oSettings.aanFeatures.p;
            if (iPages <= 1) { // hide paging when we can't page
                $(an).hide();
            } else {
                var disableClasses = calcDisableClasses(oSettings);
                $(an).find('.' + firstClassName)
                    .removeClass(oSettings.oClasses.sPageButtonDisabled)
                    .addClass(disableClasses[firstClassName]);

                // Enable/Disable `prev` button.
                $(an).find('.' + previousClassName)
                    .removeClass(oSettings.oClasses.sPageButtonDisabled)
                    .addClass(disableClasses[previousClassName]);

                // Enable/Disable `next` button.
                $(an).find('.' + nextClassName)
                    .removeClass(oSettings.oClasses.sPageButtonDisabled)
                    .addClass(disableClasses[nextClassName]);

                // Enable/Disable `last` button.
                $(an).find('.' + lastClassName)
                    .removeClass(oSettings.oClasses.sPageButtonDisabled)
                    .addClass(disableClasses[lastClassName]);

                $(an).show();
                /* Loop over each instance of the pager */
                for (var i = 0, iLen = an.length; i < iLen; i++) {
                    $(an[i]).find('#' + oSettings.sTableId + '_inputli').children().contents().each(function () {
                        if (this.nodeType === 3)
                            this.nodeValue = this.nodeValue.replace(/of \d+/, 'of ' + iPages);
                    });
                    var inputs = an[i].getElementsByTagName('input');
                    inputs[0].value = iCurrentPage;
                }
            }
        }
    };
})(jQuery);
