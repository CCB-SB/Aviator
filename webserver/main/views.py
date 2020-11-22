from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.cache import cache_page
from .models import Website, WebsiteCall, Publication
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
    context = {'search_column': -1, 'search_string':''}
    if request.method == 'POST':
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
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


def autocomplete(request):
    def remove_duplicates(x):
      return list(dict.fromkeys(x))
    if 'q' in request.GET:
        #Apply Filters
        qs = Publication.objects.all().prefetch_related('websites')
        filter = {}
        if '0' in request.GET:
            qs = qs.filter(title__icontains=request.GET.get('0'))
        if '1' in request.GET:
            qs = qs.filter(websites__status__icontains=request.GET.get('1'))
        if '2' in request.GET:
            qs = qs.filter(websites__percentage__icontains=request.GET.get('2'))
        if '3' in request.GET:
            qs = qs.filter(authors__icontains=request.GET.get('3'))
        if '4' in request.GET:
            min_str, max_str = request.GET['4'].split(';')
            if min_str and max_str:
                qs = qs.filter(authors__range=(float(min_str), float(max_str)))
            elif min_str:
                qs = qs.filter(authors__gte=float(min_str))
            elif max_str:
                qs = qs.filter(authors__lte=float(max_str))
        if '5' in request.GET:
            qs = qs.filter(journal__icontains=request.GET.get('5'))
        if '6' in request.GET:
            qs = qs.filter(pubmed_id__icontains=request.GET.get('6'))
        if '7' in request.GET:
            qs = qs.filter(abstract__icontains=request.GET.get('7'))
        #if '8' in request.GET:
        #    qs = qs.filter(websites__original_url__icontains=request.GET.get('8'))
        #if '9' in request.GET:
        #    qs = qs.filter(websites__derived_url__icontains=request.GET.get('9'))
        if '10' in request.GET:
            qs = qs.filter(contact_mail__icontains=request.GET.get('10'))
        if '11' in request.GET:
            qs = qs.filter(user_kwds__icontains=request.GET.get('11'))
        #if '12' in request.GET:
        #    qs = qs.filter(websites__pk__icontains=request.GET.get('12'))
        #Get autocomplete strings
        if request.GET.get('q') == '0':
            qs = qs.filter(title__icontains=request.GET.get('0')).order_by('title')
            return JsonResponse(remove_duplicates(list(hm['title'] for hm in qs.values('title')[:5])), safe=False)
        if request.GET.get('q') == '3':
            qs = qs.filter(authors__icontains=(request.GET.get('3')))
            search = request.GET.get('3').lower()
            seen = {}
            sList = []
            values = qs.values('authors')
            for hm in values:
                for item in hm['authors']:
                    if item in seen: continue
                    if not item.lower().startswith(search): continue
                    seen[item] = 1
                    sList.append(item)
            sList.sort();
            return JsonResponse(sList[:5], safe=False)
        if request.GET.get('q') == '5':
            qs = qs.filter(journal__icontains=request.GET.get('5')).order_by('journal')
            return JsonResponse(remove_duplicates(list(hm['journal'] for hm in qs.values('journal')))[:5], safe=False)
        if request.GET.get('q') == '6':
            qs = qs.filter(pubmed_id__istartswith=request.GET.get('6')).order_by('pubmed_id')
            return JsonResponse(remove_duplicates(list(hm['pubmed_id'] for hm in qs.values('pubmed_id')[:5])), safe=False)
        if request.GET.get('q') == '7':
            qs = qs.filter(abstract__icontains=request.GET.get('7')).order_by('abstract')
            return JsonResponse(remove_duplicates(list(hm['abstract'] for hm in qs.values('abstract')[:5])), safe=False)
        if request.GET.get('q') == '8':
            qs = qs.filter(websites__original_url__icontains=request.GET.get('8'))
            return JsonResponse(list(hm['websites__original_url'] for hm in qs.values('websites__original_url')[:5]), safe=False)
        if request.GET.get('q') == '9':
            qs = qs.filter(websites__derived_url__icontains=request.GET.get('9'))
            return JsonResponse(remove_duplicates(list(hm['websites__derived_url'] for hm in qs.values('websites__derived_url')[:5])), safe=False)
        if request.GET.get('q') == '10':
            qs = qs.filter(contact_mail__icontains=request.GET.get('10'))
            return JsonResponse(remove_duplicates(list(hm['contact_mail'].replace("@", "[at]") for hm in qs.values('contact_mail')[:5])), safe=False)
        if request.GET.get('q') == '11':
            qs = qs.filter(user_kwds__icontains=request.GET.get('11'))
            return JsonResponse(remove_duplicates(list(hm['user_kwds'] for hm in qs.values('user_kwds')[:5])), safe=False)
    return JsonResponse(list(), safe=False)

class Table(BaseDatatableView):
    model = Publication
    columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id', 'abstract', 'original_url', 'derived_url', 'contact_mail', 'user_kwds', 'scripts', 'ssl', 'website_pks']
    max_display_length = 500

    escape_values = False

    def get_initial_queryset(self):
        return Publication.objects.all().prefetch_related('websites').annotate(status=ArrayAgg('websites__status'), percentage=ArrayAgg('websites__percentage'), original_url=ArrayAgg('websites__original_url'), derived_url=ArrayAgg('websites__derived_url'), scripts=ArrayAgg('websites__script'), ssl=ArrayAgg('websites__certificate_secure'), website_pks=ArrayAgg('websites__pk'))

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

            # apply ordering
            qs = self.ordering(qs)

            website_states = list(qs.values('pubmed_id', 'websites__states'))

            # apply pagintion
            qs = self.paging(qs)

            # prepare output data
            if self.pre_camel_case_notation:
                aaData = self.prepare_results(qs)

                ret = {'sEcho': int(self._querydict.get('sEcho', 0)),
                       'iTotalRecords': total_records,
                       'iTotalDisplayRecords': total_display_records,
                       'website_states': website_states,
                       'aaData': aaData
                       }
            else:
                data = self.prepare_results(qs)

                ret = {'draw': int(self._querydict.get('draw', 0)),
                       'recordsTotal': total_records,
                       'recordsFiltered': total_display_records,
                       'website_states': website_states,
                       'data': data
                       }
            return ret
        except Exception as e:
            return self.handle_exception(e)

