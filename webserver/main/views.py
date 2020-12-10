from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.cache import cache_page
from .models import Website, WebsiteCall, Publication, Tag, CuratedWebsite
from collections import Counter
from django.core import serializers
from django.core.paginator import Paginator
from datetime import timedelta, date, datetime
import json
from django.contrib.postgres.aggregates import ArrayAgg, BoolOr, StringAgg
from django.db.models.functions import TruncDate, Cast
from django.db.models import Count
from django.db import models
from cache_memoize import cache_memoize
from django.conf import settings
from django_datatables_view.base_datatable_view import BaseDatatableView
from django.views.generic import TemplateView
from .stats import get_index_stats, get_all_statistics
from itertools import chain
import csv


def export_publications_csv(request):
    qs = Publication.objects.all().prefetch_related('websites').annotate(
        status=ArrayAgg('websites__status'), percentage=ArrayAgg('websites__percentage'),
        original_url=ArrayAgg('websites__original_url'),
        derived_url=ArrayAgg('websites__derived_url'), scripts=ArrayAgg('websites__script'),
        ssl=ArrayAgg('websites__certificate_secure'), heap_size=ArrayAgg('websites__last_heap_size'),
        website_pks=ArrayAgg('websites__pk'))
    columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id',
               'abstract', 'original_url', 'derived_url', 'contact_mail', 'user_kwds',
               'scripts', 'ssl', 'heap_size', 'website_pks']
    numeric_cols = {4}
    email_col = {10}
    return prepare_csv_export(qs, columns, request, numeric_cols, email_col)


def export_curated_csv(request):
    columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id', 'description', 'url', 'tag_tags', 'website']
    qs = CuratedWebsite.objects.all().prefetch_related('tags').annotate(tag_tags=ArrayAgg('tags__name'))
    numeric_cols = {4}
    email_col = {}
    return prepare_csv_export(qs, columns, request, numeric_cols, email_col)


def prepare_csv_export(qs, columns, request, numeric_cols, email_col):
    response = HttpResponse(content_type='text/csv')
    writer = csv.writer(response, delimiter ='\t')
    header_line = []
    for k in range(0, len(columns)):
        if str(k) in request.GET:
            header_line.append(columns[k])
    writer.writerow(header_line)
    filter = {}
    for k in range(0, len(columns)):
        if str(k) in request.GET and len(columns) > k:
            filter_string = request.GET.get(str(k))
            if len(filter_string) > 0 and filter_string != ";":
                field_name = columns[k]
                if k in email_col:
                    filter_string = filter_string.replace("[at]", "@")
                if k in numeric_cols:
                    min_str, max_str = filter_string.split(';')
                    if min_str and max_str:
                        filter["{}__range".format(field_name)] = (float(min_str), float(max_str))
                    elif min_str:
                        filter["{}__gte".format(field_name)] = float(min_str)
                    elif max_str:
                        filter["{}__lte".format(field_name)] = float(max_str)
                else:
                    filter["{}__icontains".format(field_name)] = filter_string
    if filter:
        qs = qs.filter(**filter)
    csv_content = []
    for k in range(0, len(columns)):
        if str(k) in request.GET:
            field_name = columns[k]
            values = qs.values_list(field_name)
            if k in email_col:
                new_values = []
                for i in range(0, len(values)):
                    new_value = []
                    for j in range(0, len(values[i][0])):
                        new_value.append(str(values[i][0][j]).replace("@", "[at]"))
                    new_values.append(str(new_value))
                csv_content.append(new_values)
            else:
                csv_content.append(values)
    writer.writerows(list(map(list, zip(*csv_content))))
    response['Content-Disposition'] = 'attachment; filename="export.csv"'
    return response


def index(request):
    context = get_index_stats()
    return render(request, 'index.html', context)


def overview(request):
    context = {'search_column': 0, 'search_string': ''}
    if request.method == 'POST':
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
    return render(request, 'overview.html', context)


@cache_memoize(settings.CACHE_TIMEOUT)
def get_publication_datatable_info():
    context = {}
    context["websites"] = json.dumps({x['pk']: x for x in list(
        Website.objects.all().values('pk', 'status', 'original_url', 'derived_url', 'percentage',
                                     'states'))})
    return context


def publications(request):
    context = {'search_column': -1, 'search_string': ''}
    if request.method == 'POST':
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
    return render(request, 'publications.html', context)

def curated(request):
    context = {'search_column': -1, 'search_string': ''}
    if request.method == 'POST':
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
    return render(request, 'curated.html', context)

def details(request, pk):
    context = {}
    website = get_object_or_404(Website, pk=pk)
    context['calls'] = website.calls.all().order_by('datetime')
    context['website'] = website
    return render(request, 'details.html', context)


def publication(request, pk):
    context = {}
    paper = get_object_or_404(Publication, pk=pk)
    context['paper'] = paper
    context['websites'] = paper.websites.all()
    return render(request, 'publication.html', context)


def author(request):
    context = {}
    context['websites'] = Website.objects.all()
    context['tags'] = Tag.objects.all()
    return render(request, 'author.html', context)

def api(request):
    context = {}
    return render(request, 'api.html', context)

# cache for 6 hours
@cache_page(60 * 60 * 6)
def websiteData(request):
    return JsonResponse({"data": list(
        Website.objects.all().values('original_url', 'derived_url', 'status', 'created_at',
                                     'updated_at', 'pk', 'papers'))})


# cache for 6 hours
@cache_page(60 * 60 * 6)
def paperData(request):
    data_papers = list(
        Publication.objects.all().values('pk', 'title', 'url', 'authors', 'abstract', 'year',
                                         'journal', 'pubmed_id', 'contact_mail',
                                         'user_kwds').annotate(websites=ArrayAgg('websites')))
    return JsonResponse({"data": data_papers})


def statistics(request):
    context = get_all_statistics(Publication.objects.all())
    return render(request, 'statistics.html', context)


def autocomplete(request):
    def remove_duplicates(x):
        return list(dict.fromkeys(x))
    if 'q' in request.GET:
        qs = Publication.objects.all().prefetch_related('websites').annotate(
            status=ArrayAgg('websites__status'), percentage=ArrayAgg('websites__percentage'),
            original_url=ArrayAgg('websites__original_url'),
            derived_url=ArrayAgg('websites__derived_url'), scripts=ArrayAgg('websites__script'),
            ssl=ArrayAgg('websites__certificate_secure'), heap_size=ArrayAgg('websites__last_heap_size'),
            website_pks=ArrayAgg('websites__pk'))
        columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id',
                   'abstract', 'original_url', 'derived_url', 'contact_mail', 'user_kwds',
                   'scripts', 'ssl', 'heap_size', 'website_pks']
        numeric_cols = {4}
        email_col = {10}
        filter = {}
        for k in range(0, len(columns)):
            if str(k) in request.GET and len(columns) > k:
                filter_string = request.GET.get(str(k))
                if len(filter_string) > 0 and filter_string != ";":
                    field_name = columns[k]
                    if k in email_col:
                        filter_string = filter_string.replace("[at]", "@")
                    if k in numeric_cols:
                        min_str, max_str = filter_string.split(';')
                        if min_str and max_str:
                            filter["{}__range".format(field_name)] = (float(min_str), float(max_str))
                        elif min_str:
                            filter["{}__gte".format(field_name)] = float(min_str)
                        elif max_str:
                            filter["{}__lte".format(field_name)] = float(max_str)
                    else:
                        filter["{}__icontains".format(field_name)] = filter_string
        if filter:
            qs = qs.filter(**filter)
        listed_cols = {1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15}
        ignore_cols = {4, 15}
        field_name = columns[int(request.GET.get('q'))]
        if int(request.GET.get('q')) in ignore_cols:
            return JsonResponse(list(), safe=False)
        elif int(request.GET.get('q')) in listed_cols:
            search = request.GET.get(request.GET.get('q')).lower()
            seen = {}
            sList = []
            values = qs.values(field_name)
            for hm in values:
                for item in hm[field_name]:
                    if item in seen: continue
                    if search not in item.lower(): continue
                    seen[item] = 1
                    sList.append(item)
            sList.sort()
            sList = sList[:5]
            if int(request.GET.get('q')) in email_col:
                for i in range(len(sList)):
                    sList[i] = sList[i].replace("@", "[at]")
            return JsonResponse(sList, safe=False)
        else:
            return JsonResponse(remove_duplicates(list(hm[field_name] for hm in qs.order_by(field_name).values(field_name)))[:5], safe=False)
    return JsonResponse(list(), safe=False)

def curated_autocomplete(request):
    def remove_duplicates(x):
        return list(dict.fromkeys(x))
    if 'q' in request.GET:
        columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id', 'description', 'url', 'tag_tags', 'website']
        qs = CuratedWebsite.objects.all().prefetch_related('tags').annotate(tag_tags=ArrayAgg('tags__name'))
        numeric_cols = {4}
        email_col = {}
        filter = {}
        for k in range(0, len(columns)):
            if str(k) in request.GET and len(columns) > k:
                filter_string = request.GET.get(str(k))
                if len(filter_string) > 0 and filter_string != ";":
                    field_name = columns[k]
                    if k in email_col:
                        filter_string = filter_string.replace("[at]", "@")
                    if k in numeric_cols:
                        min_str, max_str = filter_string.split(';')
                        if min_str and max_str:
                            filter["{}__range".format(field_name)] = (float(min_str), float(max_str))
                        elif min_str:
                            filter["{}__gte".format(field_name)] = float(min_str)
                        elif max_str:
                            filter["{}__lte".format(field_name)] = float(max_str)
                    else:
                        filter["{}__icontains".format(field_name)] = filter_string
        if filter:
            qs = qs.filter(**filter)
        listed_cols = {3, 9}
        ignore_cols = {10}
        field_name = columns[int(request.GET.get('q'))]
        if int(request.GET.get('q')) in ignore_cols:
            return JsonResponse(list(), safe=False)
        elif int(request.GET.get('q')) in listed_cols:
            search = request.GET.get(request.GET.get('q')).lower()
            seen = {}
            sList = []
            values = qs.values(field_name)
            for hm in values:
                for item in hm[field_name]:
                    if item in seen: continue
                    if search not in item.lower(): continue
                    seen[item] = 1
                    sList.append(item)
            sList.sort()
            sList = sList[:5]
            if int(request.GET.get('q')) in email_col:
                for i in range(len(sList)):
                    sList[i] = sList[i].replace("@", "[at]")
            return JsonResponse(sList, safe=False)
        else:
            return JsonResponse(remove_duplicates(list(hm[field_name] for hm in qs.order_by(field_name).values(field_name)))[:5], safe=False)
    return JsonResponse(list(), safe=False)

class Table(BaseDatatableView):
    model = Publication
    columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id',
               'abstract', 'original_url', 'derived_url', 'contact_mail', 'user_kwds',
                'scripts', 'ssl', 'heap_size', 'website_pks']
    max_display_length = 500

    escape_values = False

    def get_initial_queryset(self):
        return Publication.objects.all().prefetch_related('websites').annotate(
            status=ArrayAgg('websites__status'), percentage=ArrayAgg('websites__percentage'),
            original_url=ArrayAgg('websites__original_url'),
            derived_url=ArrayAgg('websites__derived_url'), scripts=ArrayAgg('websites__script'),
            ssl=ArrayAgg('websites__certificate_secure'), heap_size=ArrayAgg('websites__last_heap_size'),
            website_pks=ArrayAgg('websites__pk'))

    def render_column(self, row, column):
        # We want to render user as a custom column
        if column == "contact_mail":
            return list(c.replace("@", "[at]") for c in row.contact_mail)
        return super(Table, self).render_column(row, column)

    def filter_queryset(self, qs):
        filter_keys = [k for k in self.request.GET.keys() if "[value]" in k and "columns" in k]
        numeric_cols = {4}
        filter = {}
        for k in filter_keys:
            if len(self.request.GET[k]) > 0:
                col = int(k.split('[')[1].split(']')[0])
                field_name = self.columns[col]
                if col in numeric_cols:
                    min_str, max_str = self.request.GET[k].split(';')
                    if min_str and max_str:
                        filter["{}__range".format(field_name)] = (float(min_str), float(max_str))
                    elif min_str:
                        filter["{}__gte".format(field_name)] = float(min_str)
                    elif max_str:
                        filter["{}__lte".format(field_name)] = float(max_str)
                else:
                    filter["{}__icontains".format(field_name)] = self.request.GET[k]
        if filter:
            qs = qs.filter(**filter)
        return qs

    def prepare_results(self, qs):
        data = []
        for item in qs:
            data.append({column: self.render_column(item, column) for column in self.get_columns()})
        return data

    def extract_datatables_column_data(self):
        """ Helper method to extract columns data from request as passed by Datatables 1.10+
        """
        request_dict = self._querydict
        col_data = []
        if not self.pre_camel_case_notation:
            counter = 0
            data_name_key = 'columns[{0}][name]'.format(counter)
            while data_name_key in request_dict:
                searchable = True if request_dict.get(
                    'columns[{0}][searchable]'.format(counter)) == 'true' else False
                orderable = True if request_dict.get(
                    'columns[{0}][orderable]'.format(counter)) == 'true' else False

                col_data.append({'name': request_dict.get(data_name_key),
                                 'data': request_dict.get('columns[{0}][data]'.format(counter)),
                                 'searchable': searchable,
                                 'orderable': orderable,
                                 'search.value': request_dict.get(
                                     'columns[{0}][search][value]'.format(counter)),
                                 'search.regex': request_dict.get(
                                     'columns[{0}][search][regex]'.format(counter)),
                                 })
                counter += 1
                data_name_key = 'columns[{0}][name]'.format(counter)
        return col_data

    def ordering(self, qs):
        """ Get parameters from the request and prepare order by clause
        """
        # Number of columns that are used in sorting
        sorting_cols = 0
        if self.pre_camel_case_notation:
            try:
                sorting_cols = int(self._querydict.get('iSortingCols', 0))
            except ValueError:
                sorting_cols = 0
        else:
            sort_key = 'order[{0}][column]'.format(sorting_cols)
            while sort_key in self._querydict:
                sorting_cols += 1
                sort_key = 'order[{0}][column]'.format(sorting_cols)

        order = []
        order_columns = self.get_order_columns()

        for i in range(sorting_cols):
            # sorting column
            sort_dir = 'asc'
            try:
                if self.pre_camel_case_notation:
                    sort_col = int(self._querydict.get('iSortCol_{0}'.format(i)))
                    # sorting order
                    sort_dir = self._querydict.get('sSortDir_{0}'.format(i))
                else:
                    sort_col = int(self._querydict.get('order[{0}][column]'.format(i)))
                    # sorting order
                    sort_dir = self._querydict.get('order[{0}][dir]'.format(i))
            except ValueError:
                sort_col = 0

            sdir = '-' if sort_dir == 'desc' else ''
            sortcol = order_columns[sort_col]

            if isinstance(sortcol, list):
                for sc in sortcol:
                    if isinstance(sc, list):
                        for o in sc:
                            order.append('{0}{1}'.format(sdir, o.replace('.', '__')))
                    else:
                        order.append('{0}{1}'.format(sdir, sc.replace('.', '__')))
            else:
                if isinstance(sortcol, list):
                    for o in sortcol:
                        order.append('{0}{1}'.format(sdir, o.replace('.', '__')))
                else:
                    order.append('{0}{1}'.format(sdir, sortcol.replace('.', '__')))
        if order:
            return qs.order_by(*order)
        return qs

    def get_context_data(self, *args, **kwargs):
        try:
            self.initialize(*args, **kwargs)

            # prepare columns data (for DataTables 1.10+)
            self.columns_data = self.extract_datatables_column_data()

            # determine the response type based on the 'data' field passed from JavaScript
            # https://datatables.net/reference/option/columns.data
            # col['data'] can be an integer (return list) or string (return dictionary)
            # we only check for the first column definition here as there is no way to return list and dictionary
            # at once
            self.is_data_list = True
            if self.columns_data:
                self.is_data_list = False
                try:
                    int(self.columns_data[0]['data'])
                    self.is_data_list = True
                except ValueError:
                    pass

            # prepare list of columns to be returned
            self._columns = self.get_columns()

            # prepare initial queryset
            qs = self.get_initial_queryset()

            # store the total number of records (before filtering)
            total_records = qs.count()

            # apply filters
            qs = self.filter_queryset(qs)

            # number of records after filtering
            total_display_records = qs.count()

            stats = get_all_statistics(qs)

            # apply ordering
            qs = self.ordering(qs)

            website_states = list(qs.values('pubmed_id', 'websites__states'))
            latest_date = WebsiteCall.objects.latest("datetime").datetime

            state_dates = [(latest_date-timedelta(days=day_delta)).date().strftime("%Y-%m-%d") for day_delta in range(settings.TEMPORAL_INFO_DAYS, -1, -1)]

            # apply pagintion
            qs = self.paging(qs)

            # prepare output data
            if self.pre_camel_case_notation:
                aaData = self.prepare_results(qs)

                ret = {'sEcho': int(self._querydict.get('sEcho', 0)),
                       'iTotalRecords': total_records,
                       'iTotalDisplayRecords': total_display_records,
                       'aaData': aaData
                       }
            else:
                data = self.prepare_results(qs)

                ret = {'draw': int(self._querydict.get('draw', 0)),
                       'recordsTotal': total_records,
                       'recordsFiltered': total_display_records,
                       'data': data
                       }
            ret['website_states'] = {
                           "states": website_states,
                           "dates": state_dates
                       }
            ret['statistics'] = stats
            return ret
        except Exception as e:
            return self.handle_exception(e)

class CuratedTable(Table):
    model = CuratedWebsite
    columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id', 'description', 'url', 'tag_tags', 'website_pk']
    max_display_length = 500

    escape_values = False

    def get_initial_queryset(self):
        return CuratedWebsite.objects.all().prefetch_related('tags').annotate(tag_tags=ArrayAgg('tags__name'), website_pk=ArrayAgg('website'))

    def get_context_data(self, *args, **kwargs):
        try:
            self.initialize(*args, **kwargs)

            # prepare columns data (for DataTables 1.10+)
            self.columns_data = self.extract_datatables_column_data()

            # determine the response type based on the 'data' field passed from JavaScript
            # https://datatables.net/reference/option/columns.data
            # col['data'] can be an integer (return list) or string (return dictionary)
            # we only check for the first column definition here as there is no way to return list and dictionary
            # at once
            self.is_data_list = True
            if self.columns_data:
                self.is_data_list = False
                try:
                    int(self.columns_data[0]['data'])
                    self.is_data_list = True
                except ValueError:
                    pass

            # prepare list of columns to be returned
            self._columns = self.get_columns()

            # prepare initial queryset
            qs = self.get_initial_queryset()

            # store the total number of records (before filtering)
            total_records = qs.count()

            # apply filters
            qs = self.filter_queryset(qs)

            # number of records after filtering
            total_display_records = qs.count()

            stats = get_all_statistics(qs, True)

            # apply ordering
            qs = self.ordering(qs)

            website_states = list(qs.values('pubmed_id', 'states'))
            latest_date = WebsiteCall.objects.latest("datetime").datetime
            if CuratedWebsite.objects.all().count() > 0:
                dates = CuratedWebsite.objects.all()[0].dates
                latest_date = dates[len(dates) - 1]

            state_dates = [(latest_date-timedelta(days=day_delta)).date().strftime("%Y-%m-%d") for day_delta in range(settings.TEMPORAL_INFO_DAYS, -1, -1)]

            # apply pagintion
            qs = self.paging(qs)

            # prepare output data
            if self.pre_camel_case_notation:
                aaData = self.prepare_results(qs)

                ret = {'sEcho': int(self._querydict.get('sEcho', 0)),
                       'iTotalRecords': total_records,
                       'iTotalDisplayRecords': total_display_records,
                       'aaData': aaData
                       }
            else:
                data = self.prepare_results(qs)

                ret = {'draw': int(self._querydict.get('draw', 0)),
                       'recordsTotal': total_records,
                       'recordsFiltered': total_display_records,
                       'data': data
                       }
            ret['website_states'] = {
                           "states": website_states,
                           "dates": state_dates
                       }
            ret['statistics'] = stats
            return ret
        except Exception as e:
            return self.handle_exception(e)
