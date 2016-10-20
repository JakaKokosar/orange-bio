"""GenAPI"""
import json
import io
import gzip
import requests
import re


from math import ceil
from collections import OrderedDict
from genesis import Genesis, GenData
from Orange.data import ContinuousVariable, StringVariable, Domain, Table, DiscreteVariable
from .tools import transpose_table
from .. import dicty

CallBack = dicty.CallBack
transformValues = dicty.transformValues
averageAttributes = dicty.averageAttributes
example_tables = dicty.example_tables

view_model = ['experiment', 'growth', 'genotype', 'treatment', 'strain', 'time', 'replicate']
DEFAULT_EMAIL = 'anonymous@genialis.com'
DEFAULT_PASSWD = 'anonymous'
DEFAULT_URL = 'https://dictyexpress.research.bcm.edu'


class GenAPI(object):
    """
    Python module that leverages Genesis PyAPI (Python API for accsess to DictyExpress database).
    It supports connection to the server and data retrieval functionalities.
    """
    def __init__(self, email=DEFAULT_EMAIL, password=DEFAULT_PASSWD, url=DEFAULT_URL):

        self._gen = Genesis(email, password, url)
        self.email = email

    def fetch_etc_objects(self, callback=lambda: None):
        """ Function downloads all available :obj:`GenData` etc objects from DictyExpress database.

        Returns:
            :obj:`list`: of :obj:`GenData` objects

        """
        cbc = CallBack(1, callback)
        try:
            list_of_experiments = self._gen.api.data.get(case_ids__contains='5535115cfad58d5e03006217', status='done',
                                                         type__startswith='data:etc:')['objects']

        except requests.exceptions.ConnectionError as e:
            raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

        cbc()
        store_experiments = [GenData(exp, self._gen) for exp in list_of_experiments]
        cbc.end()
        return store_experiments

    def download_etc_data(self, gen_data_id, callback=lambda: None):
        """ Function downloads etc data of a chosen experiment from the server.

        Args:
            gen_data_id (str): id of :obj:`GenData` object to download.

        Returns:
             :obj:`dict`: data in json like format


        """
        cbc = CallBack(3, callback, callbacks=70)

        try:
            response = next(self._gen.download([gen_data_id], 'output.etcfile'))
        except requests.exceptions.ConnectionError as e:
            raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

        cbc()
        if not response.ok:
            response.raise_for_status()

        response_gzipped = io.BytesIO(response.content)
        cbc()
        response_content = io.TextIOWrapper(gzip.GzipFile(fileobj=response_gzipped), encoding="utf-8")
        cbc()

        try:
            json_object = json.load(response_content)
        except ValueError as e:
            raise ValueError('Downloaded data is not a valid JSON') from e

        cbc.end()
        return json_object

    def etc_to_table(self, etc_json, time_var=False, callback=lambda: None):
        """ Converts data from Json to :obj:`Orange.data.table`

        Args:
            etc_json (dict): Data in json like format
            time_var (bool): Create column of time points. Default is set to False.
        Returns:
            :obj:`Orange.data.Table`
        """
        cbc = CallBack(2, callback, callbacks=30)

        variables = []
        time_point = 1
        for time in etc_json['etc']['timePoints']:
            var = ContinuousVariable('TP ' + str(time_point))
            var.attributes['Time'] = str(time)
            variables.append(var)
            time_point += 1

        meta_attr = StringVariable.make('Gene')
        domain = Domain(variables, metas=[meta_attr])
        cbc()

        table = []
        for row in etc_json['etc']['genes']:
            gene_expression = [exp for exp in etc_json['etc']['genes'][row]]
            gene_expression.append(row)
            table.append(gene_expression)

        orange_table = Table(domain, table)

        if time_var:
            orange_table = transpose_table(orange_table)
            cbc()

        cbc.end()
        return orange_table

    def download_exprs_data(self, ids, typ, callback=lambda: None):
        """ Download data from the server.

        Args:
            ids (list): List of expression IDs to download

            typ (str): Expression type to download

        Returns:
             :obj:`dict`: Keys are genes, values are lists of time-point values

        """
        out = OrderedDict()
        cbc = CallBack(len(ids), callback, callbacks=80)
        try:
            for response in self._gen.download(ids, 'output.' + typ):
                if not response.ok:
                    response.raise_for_status()

                response_gzipped = io.BytesIO(response.content)
                response_content = io.TextIOWrapper(gzip.GzipFile(fileobj=response_gzipped), encoding="utf-8")
                cbc()

                for l in response_content.read().split('\n')[1:]:
                    if l:
                        gene, val = l.split('\t')
                        # is there OrderDict that can have default value?
                        if str(gene) in out:
                            out[str(gene)].append(str(val))
                        else:
                            out[str(gene)] = [str(val)]

        except requests.exceptions.ConnectionError as e:
            raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

        cbc.end()
        return out

    def get_filters(self):
        """ Sets filters. For now it supports filtering by projects.
        """
        # For now it supports filtering by projects:
        filters = dict()
        try:
            all_projects = self._gen.api.case.get()['objects']
        except requests.exceptions.ConnectionError as e:
            raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

        filters['Projects'] = [(project['name'], project['id']) for project in all_projects]
        filters['Projects'].sort(key=lambda x: x[0])
        return filters

    def get_annotations(self, project_id, offset, callback=lambda: None):
        """ Get annotations for each sample in selected project.

        Args:
            project_id (str):  id of a project that contains expressions

        Returns:
            :obj:`dict`: Dictionary of annotations for each sample in provided project. Key is an ID of sample and
            value is a list of its annotations.

        """
        if project_id:
            cbc = CallBack(3, callback, callbacks=100)

            if project_id == "all":
                project_id = None  # We want to download data from all projects

            try:
                exprs = self._gen.api.data.get(type__startswith='data:expression:',
                                               case_ids__contains=project_id,
                                               offset=offset,
                                               limit=100)
                if not exprs['objects']:
                    cbc.end()
                    raise ValueError('Sorry no data for selected experiment')

                """
                The problem is that the data objects do not include annotations.
                Annotations were attached to the reads objects that are on inputs,
                2 levels back: reads -> mapping -> expressions
                """

                cbc()
                alignment_ids = set(d['input']['alignment'] for d in exprs['objects'])
                alignments = self._gen.api.data('set')(';'.join(alignment_ids)).get()['objects']
                cbc()

                read_ids = set(d['input']['reads'] for d in alignments)
                reads = self._gen.api.data('set')(';'.join(read_ids)).get()['objects']
                cbc()

            except requests.exceptions.ConnectionError as e:
                cbc.end()
                raise requests.exceptions.ConnectionError('Server not accessible, check your connection.') from e

            metas = exprs['meta']
            exprs_map = OrderedDict((exp['id'], exp) for exp in exprs['objects'])
            alignments_map = {alignment['id']: alignment for alignment in alignments}
            reads_map = {read['id']: read for read in reads}

            cbc.end()

            def handle_annotations(exprs_map, alignments_map, reads_map):
                """ Returns a dict where key is ID of expression object and value is a list of its annotations.

                example:('5534fee7fad58d5de20061fb', ['Dd_Dp_genome_biology', 'K.a.',
                                                                            'wildtype', 'Prestalk', 'NC4', 16, 1]
                """

                annot_out = OrderedDict()

                for exprs_id, exprs in exprs_map.items():
                    aligment_id = exprs['input']['alignment']
                    reads_id = alignments_map[aligment_id]['input']['reads']
                    anno = reads_map[reads_id]['var']

                    exp_labels = []
                    for label in view_model:
                        if label == 'experiment':
                            exp_labels.append(anno[label] if label in anno else '?')
                        elif label == 'replicate':
                            exp_labels.append(anno['replicates'][label] if label in anno['replicates'] else '?')
                        else:
                            exp_labels.append(anno['sample'][label] if label in anno['sample'] else '?')
                    exp_labels.append(exprs['date_created'])
                    annot_out[exprs_id] = exp_labels

                return annot_out

            def handle_metas(metas):
                """ Support pagination.
                """
                out = {}
                next = metas['next']
                previous = metas['previous']
                offset = int(metas['offset'])
                limit = int(metas['limit'])
                total = int(metas['total_count'])
                out['next_offset'] = None
                out['prev_offset'] = None

                if next:
                    out['next_offset'] = int(re.search('offset=(.*)', next).group(1))
                if previous:
                    out['prev_offset'] = int(re.search('offset=(.*)', previous).group(1))

                if offset == 0:
                    out['currentPage'] = 1
                else:
                    out['currentPage'] = (offset // limit) + 1

                out['allPages'] = ceil(total / limit)
                return out

            return handle_annotations(exprs_map, alignments_map, reads_map), handle_metas(metas)
