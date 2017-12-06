from unittest import TestCase
from orangecontrib.bio.ncbi.taxonomy import common_taxids, name, Taxonomy


class TaxonomyTest(TestCase):

    def setUp(self):
        self.common_ids = common_taxids()
        self.organisms = [(name(tax_id), tax_id) for tax_id in self.common_ids]
        self.taxon = Taxonomy()

    def test_ids_count(self):
        self.assertGreater(len(self.taxon.taxids()), len(self.common_ids))

    def test_common_organisms(self):
        for id in self.common_ids:
            # create taxon from common organisms
            self.taxon.get_entry(id)
