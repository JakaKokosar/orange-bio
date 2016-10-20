"""ResolweAPI"""
import io
import gzip
import requests


from resdk import Resolwe
from collections import OrderedDict
from urllib import parse as urlparse
from .. import dicty

CallBack = dicty.CallBack


class ResolweAPI(object):

    def __init__(self, user, password, url):

        """
        :param user: Username of logged in user
        :param password: Password of logged in user
        :param url: Server URL
        :param buffer: Default None
        """

        self._res = Resolwe(user, password, url)
        self.aut = self._res.auth

    def fetch_etc_objects(self):
        """
        Gets all projects available for currently logged in user. Resolwe objects are stored in buffer.

        # gets all available experiments for currently logged in user
        res_objects = self._res.data.filter()

        if self.buffer is not None:
            # store res_objects to buffer
            self.to_buffer(self.user, res_objects)

        """
        raise  NotImplementedError

    def get_json(self, object_id):

        """
        :param object_id: ID of Resolwe object to get from buffer
        :return: Returns data in json like format
        :rtype : '_io.TextIOWrapper'

        return self.from_buffer(self.user, object_id).download()
        """
        raise NotImplementedError

    def etc_to_table(self, etc_json):

        """
        :param etc_json: Data in json format from Resolwe
        :return: it returns `Orange.data.Table`


        json_object = json.load(etc_json)

        variables = []
        tp = 1
        for time in json_object['etc']['timePoints']:
            var = ContinuousVariable("TP " + str(tp))
            var.attributes["Time"] = float(time)
            variables.append(var)
            tp += 1

        meta_attr = StringVariable.make("DDB")
        domain = Domain(variables, metas=[meta_attr])

        table = []
        for row in json_object['etc']['genes']:
            gene = []
            for i in [time for time in json_object['etc']['genes'][row]]:
                gene.append(i)
            gene.append(row)
            table.append(gene)

        orangeTable = Table(domain, table)

        return orangeTable

        """
        raise NotImplementedError

    def get_data(self):
        """

        Returns:

        """
        data = self._res.data.filter(type='data:expression:')
        out = OrderedDict()

        def parse_descriptor(exprs_sample):
            view = ['name', 'experiment_type', 'molecule', 'genotype', 'cell_type', 'organism', 'source', 'date']
            anno = []
            for label in view:
                if label == 'name':
                    anno.append(exprs_sample.name)
                elif label == 'date':
                    anno.append(exprs_sample.created)
                else:
                    if label in exprs_sample.descriptor['geo']:
                        anno.append(exprs_sample.descriptor['geo'][label])
                    else:
                        anno.append('?')

            return anno

        for data_obj in data:
            exprs_sample = self._res.sample.filter(data=data_obj.id)
            if exprs_sample:
                out[str(data_obj.id)] = parse_descriptor(exprs_sample[0])

        return out

    def download_exprs_data(self, fids, ftype, callback=lambda: None):
        """
        Args:
            fids (list): list of data ids to download.
            ftype (str): expression type.

        Returns:
            :obj:`dict` keys are genes, values are lists of time-point values

        """
        out = OrderedDict()

        def get_file_type(fid, ftype):
            try:
                anno = self._res.data.get(int(fid)).annotation
                return anno['output.' + ftype]['value']['file']
            except (KeyError, TypeError):
                    raise ValueError('Selected expression type is not available')

        cbc = CallBack(len(fids), callback, callbacks=80)

        for fid in fids:
            file_type = get_file_type(fid, ftype)
            file_url = urlparse.urljoin(self._res.url, 'data/{}/{}'.format(fid, file_type))

            response = requests.get(file_url, stream=True, auth=self._res.auth)
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
        cbc.end()
        return out
