from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.cache import cache_page
from .models import Website, WebsiteCall, Publication
from collections import Counter
from django.core import serializers
from django.core.paginator import Paginator
from datetime import timedelta, date, datetime
import json
from django.contrib.postgres.aggregates import ArrayAgg, BoolOr
from django.db.models.functions import TruncDate, Cast
from django.db.models import Count
from django.db import models
from cache_memoize import cache_memoize
from django.conf import settings
from django_datatables_view.base_datatable_view import BaseDatatableView
from django.views.generic import TemplateView
from .stats import get_index_stats, get_all_statistics

def index(request):
    context = get_index_stats()
    return render(request, 'index.html', context)

def overview(request):
    context = {'search_column':0, 'search_string':''}
    if request.method == 'POST':
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
    return render(request, 'overview.html', context)

@cache_memoize(settings.CACHE_TIMEOUT)
def get_publication_datatable_info():
    context = {}
    context["websites"] = json.dumps({x['pk']: x for x in list(
        Website.objects.all().values('pk', 'status', 'original_url', 'derived_url', 'percentage', 'states'))})
    return context

def publications(request):
    context = {}
    context.update(get_publication_datatable_info())
    return render(request, 'publications.html', context)

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
    return render(request, 'author.html', context)

# cache for 6 hours
@cache_page(60 * 60 * 6)
def websiteData(request):
    return JsonResponse({"data": list(Website.objects.all().values('original_url', 'derived_url', 'status', 'created_at', 'updated_at', 'pk', 'papers'))})

# cache for 6 hours
@cache_page(60 * 60 * 6)
def paperData(request):
    data_papers = list(Publication.objects.all().values('pk', 'title', 'url', 'authors', 'abstract', 'year', 'journal', 'pubmed_id', 'contact_mail', 'user_kwds').annotate(websites=ArrayAgg('websites')))
    return JsonResponse({"data": data_papers})

def statistics(request):
    context = get_all_statistics()
    return render(request, 'statistics.html', context)

class Table(BaseDatatableView):
    model = Publication
    #columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id', 'abstract', 'original_url', 'derived_url', 'contact_mail', 'user_kwds', 'states', 'websites']
    columns = ['title', 'authors', 'year', 'journal', 'pubmed_id', 'abstract', 'contact_mail', 'user_kwds', 'websites']

    max_display_length = 500

    escape_values = False

    def render_column(self, row, column):
        # We want to render user as a custom column
        #TODO
        #if column == 'status':
        #    return list(row.websites.values_list("status", flat=True))
        #elif column == 'percentage':
        #    return list(row.websites.values_list("percentage", flat=True))
        #elif column == 'states':
        #    return list(row.websites.values_list("states", flat=True))
        #elif column == 'original_url':
        #    return list(row.websites.values_list("original_url", flat=True))
        #elif column == 'derived_url':
        #    return list(row.websites.values_list("derived_url", flat=True))
        #elif column == 'websites':
        #    return list(row.websites.values_list("pk", flat=True))
        #else:
        #    return super(Table, self).render_column(row, column)

        if column == 'websites':
            return list(row.websites.values_list("pk", flat=True))
        else:
            return super(Table, self).render_column(row, column)

    def filter_queryset(self, qs):
        filter_keys = [k for k in self.request.GET.keys() if "[value]" in k and "columns" in k]
        numeric_cols = {2, 3, 4, 7, 8}
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
                searchable = True if request_dict.get('columns[{0}][searchable]'.format(counter)) == 'true' else False
                orderable = True if request_dict.get('columns[{0}][orderable]'.format(counter)) == 'true' else False

                col_data.append({'name': request_dict.get(data_name_key),
                                 'data': request_dict.get('columns[{0}][data]'.format(counter)),
                                 'searchable': searchable,
                                 'orderable': orderable,
                                 'search.value': request_dict.get('columns[{0}][search][value]'.format(counter)),
                                 'search.regex': request_dict.get('columns[{0}][search][regex]'.format(counter)),
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

