from django.http import HttpResponse, JsonResponse
from django.shortcuts import get_object_or_404, render
from django.views.decorators.cache import cache_page
from .models import Website, WebsiteCall, Publication, Tag, CuratedWebsite, GlobalStatistics
from datetime import timedelta
import json
from django.contrib.postgres.aggregates import ArrayAgg, BoolOr
from django.db.models import Q, BooleanField
from cache_memoize import cache_memoize
from django.conf import settings
from django_datatables_view.base_datatable_view import BaseDatatableView
from .stats import get_index_stats, get_all_statistics
import csv
from .forms import CaptchaForm


def filter_queryset(qs, filter_col2v, colnames, email_field):
    numeric_fields = {"year"}
    bool_agg_fields = {"ssl": "websites__certificate_secure"}
    if hasattr(qs.model, "websites"):
        numeric_agg_field = {"percentage": "websites__percentage"}
    elif hasattr(qs.model, "website"):
        numeric_agg_field = {"percentage": "website__percentage"}
    else:
        numeric_agg_field = {}
    qs_filter = {}
    for col, v in filter_col2v.items():
        field_name = colnames[col]
        if field_name in email_field:
            v = v.replace("[at]", "@")
        if field_name in bool_agg_fields:
            agg_filter = {}
            agg_filter["{}__gte".format(bool_agg_fields[field_name])] = v == 'true'
            annot_addition = {
                "{}_filter".format(field_name): BoolOr(Q(**agg_filter),
                                                       output_field=BooleanField())
            }
            qs = qs.annotate(**annot_addition)
            qs_filter["{}_filter".format(field_name)] = True
        elif field_name in numeric_fields:
            min_str, max_str = v.split(';')
            if min_str and max_str:
                qs_filter["{}__range".format(field_name)] = (float(min_str), float(max_str))
            elif min_str:
                qs_filter["{}__gte".format(field_name)] = float(min_str)
            elif max_str:
                qs_filter["{}__lte".format(field_name)] = float(max_str)
        elif field_name in numeric_agg_field:
            min_str, max_str = v.split(';')
            if min_str or max_str:
                agg_filter = {}
                if min_str and max_str:
                    agg_filter["{}__range".format(numeric_agg_field[field_name])] = (
                        float(min_str), float(max_str))
                elif min_str:
                    agg_filter["{}__gte".format(numeric_agg_field[field_name])] = float(min_str)
                elif max_str:
                    agg_filter["{}__lte".format(numeric_agg_field[field_name])] = float(max_str)
                annot_addition = {
                    "{}_filter".format(field_name): BoolOr(Q(**agg_filter),
                                                           output_field=BooleanField())
                }
                qs = qs.annotate(**annot_addition)
                qs_filter["{}_filter".format(field_name)] = True
        else:
            qs_filter["{}__icontains".format(field_name)] = v
    if qs_filter:
        qs = qs.filter(**qs_filter)
    return qs


def prepare_csv_export(qs, columns, header, request, email_fields, ignore_fields, name):
    response = HttpResponse(content_type='text/csv')
    writer = csv.writer(response, delimiter=',')
    header_line = []
    shown = []
    if "columns" in request.POST:
        shown = request.POST["columns"].split(';')
    for k in range(0, len(columns)):
        if str(k) in shown and (columns[k] not in ignore_fields):
            header_line.append(header[k])
    writer.writerow(header_line)

    filter_col2v = {k: request.POST[str(k)] for k in range(0, len(columns)) if
                    str(k) in shown and (columns[k] not in ignore_fields) and request.POST[str(k)] not in {"", ";"}}
    qs = filter_queryset(qs, filter_col2v, columns, email_fields)

    csv_content = []
    for k in range(0, len(columns)):
        if str(k) in shown and (columns[k] not in ignore_fields):
            field_name = columns[k]
            values = qs.values_list(field_name)
            if field_name in email_fields:
                new_values = []
                for i in range(0, len(values)):
                    new_value = []
                    for j in range(0, len(values[i][0])):
                        new_value.append(str(values[i][0][j]).replace("@", "[at]"))
                    new_values.append(str(new_value))
                csv_content.append(new_values)
            else:
                csv_content.append([x[0] for x in values])
    writer.writerows(list(map(list, zip(*csv_content))))
    response['Content-Disposition'] = 'attachment; filename="aviator_' + name + '.csv"'
    return response


def index(request):
    context = get_index_stats()
    context['overall_calls'] = 0
    context['overall_size'] = 0
    for gs in GlobalStatistics.objects.all():
        context['overall_calls'] = "{:,}".format(gs.num_calls)
        context['overall_size'] = "{:,}".format(gs.data_size // 1000000000)
    return render(request, 'index.html', context)


def help(request):
    return render(request, 'help.html')


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
    if request.POST and 'search_column' not in request.POST:
        form = CaptchaForm(request.POST)
        if form.is_valid():
            qs = Publication.objects.all().prefetch_related('websites').annotate(
                status=ArrayAgg('websites__status'), percentage=ArrayAgg('websites__percentage'),
                original_url=ArrayAgg('websites__original_url'),
                derived_url=ArrayAgg('websites__derived_url'), scripts=ArrayAgg('websites__script'),
                ssl=ArrayAgg('websites__certificate_secure'),
                heap_size=ArrayAgg('websites__last_heap_size'),
                website_pks=ArrayAgg('websites__pk'))
            columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id',
                       'abstract', 'original_url', 'derived_url', 'contact_mail', 'user_kwds',
                       'scripts', 'ssl', 'heap_size', 'website_pks']
            header = ['Title', 'Status', 'Last 30 days', 'Authors', 'Year', 'Journal', 'PubMed',
                      'Abstract',
                      'Original URL', 'Derived URL', 'Contact Mail', 'Keywords',
                      'Programming Languages', 'SSL', 'RAM Usage', 'Websites']
            email_fields = {"contact_mail"}
            ignore_fields = {"website_pks"}
            return prepare_csv_export(qs, columns, header, request, email_fields,
                                      ignore_fields, "publications")
        else:
            return HttpResponse("Wrong Captcha")
    else:
        form = CaptchaForm()
    context = {'search_column': -1, 'search_string': '', 'form': form}
    if request.method == 'POST' and 'search_column' in request.POST:
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
    if request.method == 'GET' and 'search_column' in request.GET:
        context['search_column'] = request.GET['search_column']
        context['search_string'] = request.GET['search_string']
    return render(request, 'publications.html', context)


def curated(request):
    if request.POST and 'search_column' not in request.POST:
        form = CaptchaForm(request.POST)
        if form.is_valid():
            columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id',
                       'description', 'url', 'tag_tags', 'website']
            header = ['Tool', 'Status', 'Last 30 days', 'Authors', 'Year', 'Journal', 'PubMed',
                      'Description', 'URL', 'Keywords', 'Website']
            qs = CuratedWebsite.objects.all().prefetch_related('tags').annotate(
                tag_tags=ArrayAgg('tags__name'))
            email_fields = {}
            ignore_fields = {"website"}
            return prepare_csv_export(qs, columns, header, request, email_fields,
                                      ignore_fields, "curated")
        else:
            return HttpResponse("Wrong Captcha")
    else:
        form = CaptchaForm()
    context = {'search_column': -1, 'search_string': '', 'form': form}
    if request.method == 'POST' and 'search_column' in request.POST:
        context['search_column'] = request.POST['search_column']
        context['search_string'] = request.POST['search_string']
    if request.method == 'GET' and 'search_column' in request.GET:
        context['search_column'] = request.GET['search_column']
        context['search_string'] = request.GET['search_string']
    return render(request, 'curated.html', context)


def details(request, pk):
    context = {}
    website = get_object_or_404(Website, pk=pk)
    show_for_x_days = 30
    latest_time = website.calls.latest("datetime").datetime
    first_call = website.calls.earliest("datetime").datetime
    context['calls'] = website.calls.filter(datetime__date__gte=latest_time-timedelta(days=show_for_x_days)).order_by('-datetime')
    context['website'] = website
    context["first_call"] = first_call
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
    context = get_all_statistics(
        Publication.objects.all().prefetch_related('websites').annotate(
            website_pks=ArrayAgg('websites'),
            status=ArrayAgg('websites__status')))
    context['weekdays_online'] = [0,0,0,0,0,0,0,0]
    context['weekdays_offline'] = [0,0,0,0,0,0,0,0]
    context['weekdays_average'] = [0,0,0,0,0,0,0]
    context['recovery_rate'] = [0,0,0,0,0,0,0]
    for gs in GlobalStatistics.objects.all():
        context['weekdays_online'] = gs.weekdays_online
        context['weekdays_offline'] = gs.weekdays_offline
        for n in range(len(context['weekdays_average'])):
            context['weekdays_average'][n] = context['weekdays_online'][0] - context['weekdays_online'][n + 1]
        #recovery_rate
        start = 0
        index = 0
        state = 0
        for n in range(len(gs.recovery_rate)):
            if n == 0:
                if gs.recovery_rate[n] == 0:
                    start = 1
            if gs.recovery_rate[n] == 0:
                if state == 0:
                    index = n
                    state = 1
            else:
                state = 0

        context['recovery_rate'] = []
        context['recovery_rate_legend'] = []
        for n in range(start, index + 1):
            context['recovery_rate'].append(gs.recovery_rate[n])
            if n == 1:
                context['recovery_rate_legend'].append(f'{n} day')
            else:
                context['recovery_rate_legend'].append(f'{n} days')
        ##
    return render(request, 'statistics.html', context)


def autocomplete(request):
    if 'q' in request.GET:
        qs = Publication.objects.all().prefetch_related('websites').annotate(
            status=ArrayAgg('websites__status'), percentage=ArrayAgg('websites__percentage'),
            original_url=ArrayAgg('websites__original_url'),
            derived_url=ArrayAgg('websites__derived_url'), scripts=ArrayAgg('websites__script'),
            ssl=ArrayAgg('websites__certificate_secure'),
            heap_size=ArrayAgg('websites__last_heap_size'),
            website_pks=ArrayAgg('websites'))
        columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id',
                   'abstract', 'original_url', 'derived_url', 'contact_mail', 'user_kwds',
                   'scripts', 'ssl', 'heap_size', 'website_pks']
        email_fields = {"contact_mail"}

        filter_col2v = {k: request.GET[str(k)] for k in range(0, len(columns)) if
                        str(k) in request.GET and request.GET[str(k)] not in {"", ";"}}
        qs = filter_queryset(qs, filter_col2v, columns, email_fields)

        listed_cols = {1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15}
        ignore_cols = {2, 4, 15}
        return generate_autocomplete_list(columns, email_fields, ignore_cols, listed_cols, qs,
                                          request)
    return JsonResponse(list(), safe=False)


def curated_autocomplete(request):
    if 'q' in request.GET:
        columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id',
                   'description', 'url', 'tag_tags', 'website']
        qs = CuratedWebsite.objects.all().prefetch_related('tags').annotate(
            tag_tags=ArrayAgg('tags__name'))
        email_fields = {}
        filter_col2v = {k: request.GET[str(k)] for k in range(0, len(columns)) if
                        str(k) in request.GET and request.GET[str(k)] not in {"", ";"}}
        qs = filter_queryset(qs, filter_col2v, columns, email_fields)

        listed_cols = {3, 9}
        ignore_cols = {10}
        return generate_autocomplete_list(columns, email_fields, ignore_cols, listed_cols, qs,
                                          request)
    return JsonResponse(list(), safe=False)


def generate_autocomplete_list(columns, email_fields, ignore_cols, listed_cols, qs,
                               request):
    def remove_duplicates(x):
        return sorted(set(e.strip() for e in x))

    column_i = int(request.GET.get('q'))
    field_name = columns[column_i]
    if column_i in ignore_cols:
        return JsonResponse(list(), safe=False)
    elif column_i in listed_cols:
        search = request.GET.get(request.GET.get('q')).lower()
        if field_name in email_fields:
            search = search.replace("[at]", "@")
        seen = {}
        sList = []
        values = qs.values(field_name)
        for hm in values:
            for item in hm[field_name]:
                if str(item) in seen: continue
                if search not in str(item).lower(): continue
                seen[item] = 1
                sList.append(item)
        sList.sort()
        sList = sList[:5]
        if field_name in email_fields:
            for i in range(len(sList)):
                sList[i] = sList[i].replace("@", "[at]")
            sList = remove_duplicates(sList)
        return JsonResponse(sList, safe=False)
    else:
        return JsonResponse(remove_duplicates(
            list(hm[field_name] for hm in qs.order_by(field_name).values(field_name)))[:5],
                            safe=False)


class Table(BaseDatatableView):
    model = Publication
    columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id',
               'abstract', 'original_url', 'derived_url', 'contact_mail', 'user_kwds',
               'scripts', 'ssl', 'heap_size', 'website_pks']
    max_display_length = 500

    escape_values = False

    def get_initial_queryset(self, without_annots=False):
        if without_annots:
            return Publication.objects.all()
        else:
            return Publication.objects.all().prefetch_related('websites').annotate(
                status=ArrayAgg('websites__status'), percentage=ArrayAgg('websites__percentage'),
                original_url=ArrayAgg('websites__original_url'),
                derived_url=ArrayAgg('websites__derived_url'), scripts=ArrayAgg('websites__script'),
                ssl=ArrayAgg('websites__certificate_secure'),
                heap_size=ArrayAgg('websites__last_heap_size'),
                website_pks=ArrayAgg('websites__pk'))

    def render_column(self, row, column):
        # We want to render user as a custom column
        if column == "contact_mail":
            return list(c.replace("@", "[at]") for c in row.contact_mail)
        return super(Table, self).render_column(row, column)

    def filter_queryset(self, qs):
        filter_keys = [k for k in self.request.GET.keys() if
                       "[value]" in k and "columns" in k and len(self.request.GET[k]) > 0]
        filter_col2v = {int(k.split('[')[1].split(']')[0]): self.request.GET[k] for k in
                        filter_keys}
        return filter_queryset(qs, filter_col2v, self.columns, email_field={'contact_mail'})

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
            total_records = self.get_initial_queryset(without_annots=True).count()

            # apply filters
            qs = self.filter_queryset(qs)

            # apply ordering
            qs = self.ordering(qs)

            stats = get_all_statistics(qs, curated=False)

            # number of records after filtering
            total_display_records = stats["paper_count"]

            # apply pagination
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
                "states": stats["website_states"],
                "dates": stats["state_dates"]
            }
            del stats["website_states"]
            del stats["state_dates"]
            ret['statistics'] = stats
            return ret
        except Exception as e:
            return self.handle_exception(e)


class CuratedTable(Table):
    model = CuratedWebsite
    columns = ['title', 'status', 'percentage', 'authors', 'year', 'journal', 'pubmed_id',
               'description', 'url', 'tag_tags', 'website_pk']
    max_display_length = 500
    escape_values = False

    def get_initial_queryset(self):
        return CuratedWebsite.objects.all().prefetch_related('tags').annotate(
            tag_tags=ArrayAgg('tags__name'), website_pk=ArrayAgg('website'))

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

            # apply ordering
            qs = self.ordering(qs)

            stats = get_all_statistics(qs, curated=True)

            # number of records after filtering
            total_display_records = stats["paper_count"]

            # apply pagination
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
                "states": stats["website_states"],
                "dates": stats["state_dates"]
            }
            del stats["website_states"]
            del stats["state_dates"]
            ret['statistics'] = stats
            return ret
        except Exception as e:
            return self.handle_exception(e)

def aviator_api(request):    
    if "input" in request.GET and request.GET["input"].isdigit():
        return HttpResponse(sum(int(e) for e in request.GET["input"]), content_type="text/plain")
    return HttpResponse()
