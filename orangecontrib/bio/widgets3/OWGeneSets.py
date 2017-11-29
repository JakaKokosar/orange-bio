""" GeneSets """
import threading
import numpy as np
import operator


from requests.exceptions import ConnectionError
from functools import reduce
from collections import defaultdict

from AnyQt.QtWidgets import (
    QTreeView, QTableView, QTreeWidget, QTreeWidgetItem, QTreeWidgetItemIterator, QButtonGroup, QGridLayout,
    QStackedWidget, QHeaderView, QCheckBox, QItemDelegate, QCompleter
)
from AnyQt.QtCore import (
    Qt, QObject, pyqtSignal, QRunnable, QSize, pyqtSlot, QModelIndex, QStringListModel, QThread, QThreadPool,
    Slot, QSortFilterProxyModel
)
from AnyQt.QtGui import (
    QBrush, QColor, QFont, QStandardItemModel, QStandardItem
)

from Orange.widgets.gui import (
    vBox, comboBox, lineEdit, ProgressBar, rubber, button, widgetBox, LinkRole, LinkStyledItemDelegate,
    auto_commit, widgetLabel
)

from Orange.widgets.widget import OWWidget, Msg
from Orange.widgets.settings import Setting, ContextSetting, DomainContextHandler
from Orange.widgets.utils.itemmodels import DomainModel
from Orange.widgets.utils.signals import Output, Input

from Orange.data import ContinuousVariable, DiscreteVariable, StringVariable, Domain, Table

from orangecontrib.bio.utils import serverfiles
from orangecontrib.bio import gene, geneset, taxonomy, utils


CATEGORY, GENES, MATCHED, TERM = range(4)
DATA_HEADER_LABELS = ["Category", "Genes", "Matched", "Term"]
HIERARCHY_HEADER_LABELS = ["Category"]


class Signals(QObject):
    finished = pyqtSignal()
    error = pyqtSignal(Exception)
    result = pyqtSignal(object)
    progress = pyqtSignal()


class Worker(QRunnable):
    """ Worker thread
    """

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = Signals()

        if self.kwargs:
            if self.kwargs['progress_callback']:
                self.kwargs['progress_callback'] = self.signals.progress

    @pyqtSlot()
    def run(self):
        try:
            result = self.fn(*self.args, **self.kwargs)
        except Exception as e:
            self.signals.error.emit(e)
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()


def hierarchy_tree(tax_id, gene_sets):
    def tree():
        return defaultdict(tree)

    collection = tree()

    def collect(col, set_hierarchy):
        if set_hierarchy:
            collect(col[set_hierarchy[0]], set_hierarchy[1:])

    for hierarchy, t_id, _ in gene_sets:
        collect(collection[t_id], hierarchy)

    return tax_id, collection[tax_id]


def download_gene_sets(tax_id, gene_sets, progress_callback):

    # get only those sets that are not already downloaded
    for hierarchy, tax_id in [(hierarchy, tax_id) for hierarchy, tax_id, local in gene_sets if not local]:

        serverfiles.localpath_download(geneset.sfdomain, geneset.filename(hierarchy, tax_id),
                                       callback=progress_callback.emit)

    return tax_id, gene_sets


class OWGeneSets(OWWidget):
    name = "Gene Sets"
    description = ""
    icon = ""
    priority = 1
    want_main_area = True

    # settings
    selected_organism = Setting(0)
    auto_commit = Setting(True)
    auto_apply = Setting(True)

    class Inputs:
        genes = Input("Genes", Table)

    class Outputs:
        gene_list = Output("Gene List", Table)
        matched_genes = Output("Matched Genes", Table)

    class Information(OWWidget.Information):
        pass

    class Error(OWWidget.Error):
        cant_reach_host = Msg("Host orange.biolab.si is unreachable.")
        cant_load_organisms = Msg("No available organisms, please check your connection.")

    def __init__(self):
        super().__init__()

        # commit
        self.commit_button = None

        # progress bar
        self.progress_bar = None

        # data
        self.input_data = None
        self.input_genes = None
        self.input_genes_matched = None
        self.input_genes_unmatched = None
        self.input_genes_mapped = None
        self.organisms = list()
        self.input_info = None

        # filter
        self.lineEdit_filter = None
        self.search_pattern = ''
        self.organism_select = None

        # data model view
        self.data_view = None
        self.data_model = None

        # gene matcher NCBI
        self.gene_matcher = None

        # filter proxy model
        self.filter_proxy_model = None

        # hierarchy widget
        self.hierarchy_widget = None
        self.hierarchy_state = None

        # threads
        self.threadpool = QThreadPool()

        # gui
        self.setup_gui()

        try:
            self.available_organisms()
        except ConnectionError:
            self.Error.cant_load_organisms()

        self.set_data(self.input_genes)
        #self.on_organism_change()

    @Inputs.genes
    def set_data(self, data):
        if data:
            self.input_data = data
            self.input_genes = set(data.metas.flat)
            self.generate_gene_sets()
        else:
            self.on_organism_change()

        self.update_info_box()

    def update_info_box(self):
        info_string = ''
        if self.input_genes:
            info_string += '{} unique gene names on input.\n'.format(len(self.input_genes))
            if self.input_genes_matched:
                ratio = (len(self.input_genes_matched) / len(self.input_genes)) * 100
                info_string += '{} ({:.2f}%) gene names matched.\n'.format(len(self.input_genes_matched), ratio)
        else:
            info_string += 'No genes on input.\n'

        self.input_info.setText(info_string)

    def progress_advance(self):
        # GUI should be updated in main thread. That's why we are calling advance method here
        if self.progress_bar:
            self.progress_bar.advance()

    def on_gene_sets_download(self, result):
        # make sure this happens in the main thread.
        # Qt insists that widgets be created within the GUI(main) thread.
        assert threading.current_thread() == threading.main_thread()

        tax_id, sets = result
        self.set_hierarchy_model(self.hierarchy_widget, *hierarchy_tree(tax_id, sets))

        self.organism_select.setEnabled(True)  # re-enable combobox
        self.progress_bar.finish()
        self.setStatusMessage('')

    def set_hierarchy_model(self, model, tax_id, sets):
        # TODO: maybe optimize this code?
        for key, value in sets.items():
            item = QTreeWidgetItem(model, [key])
            item.setFlags(item.flags() & (Qt.ItemIsUserCheckable | ~Qt.ItemIsSelectable | Qt.ItemIsEnabled))
            item.setData(0, Qt.CheckStateRole, Qt.Checked)
            item.setExpanded(True)
            item.tax_id = tax_id
            item.hierarchy = key

            if value:
                item.setFlags(item.flags() | Qt.ItemIsTristate)
                self.set_hierarchy_model(item, tax_id, value)
            else:
                if item.parent():
                    item.hierarchy = ((item.parent().hierarchy, key), tax_id)

            if not item.childCount() and not item.parent():
                item.hierarchy = ((key,), tax_id)

    def on_organism_change(self):
        if not self.organisms:
            return

        # get selected organism
        tax_id, _ = self.organisms[self.selected_organism]
        # init matcher
        self.gene_matcher = gene.matcher([gene.GMNCBI(tax_id), gene.GMKEGG(tax_id),
                                          gene.GMGO(tax_id), gene.GMDicty(tax_id)])
        # .set_targets(gene.NCBIGeneInfo(tax_id).keys())

        self.Error.clear()
        # do not allow user to change organism when download task is running
        self.organism_select.setEnabled(False)
        # reset hierarchy widget state
        self.hierarchy_widget.clear()
        # clear data view
        self.init_item_model()

        # get all gene sets for selected organism
        gene_sets = geneset.list_all(org=tax_id)
        # init progress bar
        self.progress_bar = ProgressBar(self, iterations=len(gene_sets) * 100)
        # status message
        self.setStatusMessage('downloading sets')

        worker = Worker(download_gene_sets, tax_id, gene_sets, progress_callback=True)
        worker.signals.progress.connect(self.progress_advance)
        worker.signals.result.connect(self.on_gene_sets_download)
        worker.signals.finished.connect(self.generate_gene_sets)
        worker.signals.error.connect(self.handle_error)

        # move download process to worker thread
        self.threadpool.start(worker)

    def generate_gene_sets(self):
        self.setStatusMessage('collecting sets')
        worker = Worker(geneset.collections, *self.get_selected_hierarchies())
        worker.signals.result.connect(self.display_gene_sets)
        worker.signals.error.connect(self.handle_error)

        self.threadpool.start(worker)

    def handle_error(self, ex):
        self.progress_bar.finish()
        self.setStatusMessage('')
        if isinstance(ex, ConnectionError):
            self.organism_select.setEnabled(True)  # re-enable combobox
            self.Error.cant_reach_host()

    def display_gene_sets(self, result):
        assert threading.current_thread() == threading.main_thread()

        if self.input_genes:
            all_genes = reduce(operator.ior, (set(g.genes) for g in result), set())
            self.gene_matcher.set_targets(all_genes)
            self.input_genes_mapped = [(symbol, self.gene_matcher.umatch(symbol)) for symbol in self.input_genes]
            self.input_genes_matched = set(filter(lambda x: x[1] is not None, self.input_genes_mapped))
            self.input_genes_unmatched = set(filter(lambda x: x[1] is None, self.input_genes_mapped))

        self.init_item_model()
        self.progress_bar = ProgressBar(self, iterations=len(result))
        self.update_info_box()

        for gene_set in result:
            category_column = QStandardItem()
            name_column = QStandardItem()
            matched_column = QStandardItem()
            genes_column = QStandardItem()

            category_column.setData(", ".join(gene_set.hierarchy), Qt.DisplayRole)
            name_column.setData(gene_set.name, Qt.DisplayRole)
            name_column.setData(gene_set.link, Qt.ToolTipRole)
            name_column.setData(gene_set.link, LinkRole)
            name_column.setForeground(QColor(Qt.blue))

            if self.input_genes:
                matched_set = gene_set.genes & {input_name for input_name, _ in self.input_genes_matched}
                matched_column.setData(matched_set, Qt.UserRole)
                matched_column.setData(len(matched_set), Qt.DisplayRole)

            genes_column.setData(len(gene_set.genes), Qt.DisplayRole)
            genes_column.setData(gene_set.genes, Qt.UserRole)  # store genes to get then on output on selection

            row = [category_column, genes_column, matched_column, name_column]
            self.data_model.appendRow(row)

            self.progress_advance()

        # adjust column width
        for i in range(len(DATA_HEADER_LABELS)):
            self.data_view.resizeColumnToContents(i)

        self.progress_bar.finish()
        self.setStatusMessage('')

    def get_selected_hierarchies(self):
        """ return selected hierarchy
        """
        sets_to_display = list()
        iterator = QTreeWidgetItemIterator(self.hierarchy_widget, QTreeWidgetItemIterator.Checked)

        while iterator.value():
            # note: if hierarchy value is not a tuple, then this is just top level qTreeWidgetItem that
            #       holds subcategories. We don't want to display all sets from category
            if type(iterator.value().hierarchy) is not str:
                sets_to_display.append(iterator.value().hierarchy)
            iterator += 1

        return sets_to_display

    def available_organisms(self):
        sets = geneset.list_all()
        tax_ids = set(taxonomy.common_taxids() + [tid for _, tid, _ in sets])

        for tax_id in tax_ids:
            try:
                tax_name = taxonomy.name(tax_id)
            except taxonomy.UnknownSpeciesIdentifier:
                tax_name = None

            self.organisms.append((tax_id, tax_name))

        self.organisms.sort(key=lambda x: x[1])
        self.organism_select.addItems([name for _, name in self.organisms])

    def commit(self):
        selection_model = self.data_view.selectionModel()

        if selection_model:
            genes_from_set = selection_model.selectedRows(GENES)
            matched_genes = selection_model.selectedRows(MATCHED)

            if genes_from_set:
                genes = [model_index.data(Qt.UserRole) for model_index in genes_from_set]
                data = [[gene_name] for gene_name in list(set.union(*genes))]
                domain = Domain([], metas=[StringVariable('Genes')])
                table = Table(domain, data)

                #data = Table.from_numpy(domain, X=np.array([1, 2]), metas=[StringVariable('Gene')])
                #print(np.array(set_of_genes, dtype=str)[np.newaxis].T)
                self.Outputs.gene_list.send(table)

            if matched_genes and self.input_genes:
                genes = [model_index.data(Qt.UserRole) for model_index in matched_genes]
                output_genes = [gene_name for gene_name in list(set.union(*genes))]
                gene_attribute = self.input_data.domain.metas
                selected = [i for i, ex in enumerate(self.input_data) if ex[gene_attribute[0]] in output_genes]

                #domain = Domain([], metas=[StringVariable('Genes')])
                #table = Table(domain, data)

                self.Outputs.matched_genes.send(self.input_data[selected])

    def setup_gui(self):
        # control area
        info_box = vBox(self.controlArea, 'Input info')
        self.input_info = widgetLabel(info_box)

        organism_box = vBox(self.controlArea, 'Organisms')
        self.organism_select = comboBox(organism_box, self,
                                        'selected_organism',
                                        callback=self.on_organism_change)

        hierarchy_box = widgetBox(self.controlArea, "Entity Sets")
        self.hierarchy_widget = QTreeWidget(self)
        self.hierarchy_widget.setEditTriggers(QTreeView.NoEditTriggers)
        self.hierarchy_widget.setHeaderLabels(HIERARCHY_HEADER_LABELS)
        self.hierarchy_widget.itemClicked.connect(self.generate_gene_sets)
        hierarchy_box.layout().addWidget(self.hierarchy_widget)

        self.commit_button = auto_commit(self.controlArea, self, "auto_commit", "&Commit", box=False)

        #rubber(self.controlArea)

        # main area
        self.filter_proxy_model = QSortFilterProxyModel(self.data_view)
        self.filter_proxy_model.setFilterKeyColumn(3)

        self.data_view = QTreeView()
        self.data_view.setModel(self.filter_proxy_model)
        self.data_view.setAlternatingRowColors(True)
        self.data_view.setSortingEnabled(True)
        self.data_view.setSelectionMode(QTreeView.ExtendedSelection)
        self.data_view.setEditTriggers(QTreeView.NoEditTriggers)
        self.data_view.viewport().setMouseTracking(True)
        self.data_view.setItemDelegateForColumn(TERM, LinkStyledItemDelegate(self.data_view))

        self.data_view.selectionModel().selectionChanged.connect(self.commit)

        self.lineEdit_filter = lineEdit(self.mainArea, self, 'search_pattern', 'Filter gene sets:')
        self.lineEdit_filter.setPlaceholderText('search pattern ...')
        self.lineEdit_filter.textChanged.connect(self.filter_proxy_model.setFilterRegExp)

        self.mainArea.layout().addWidget(self.data_view)

    def init_item_model(self):
        self.data_model = QStandardItemModel()
        self.data_model.setSortRole(Qt.UserRole)
        self.data_model.setHorizontalHeaderLabels(DATA_HEADER_LABELS)
        self.filter_proxy_model.setSourceModel(self.data_model)

    def sizeHint(self):
        return QSize(1280, 960)


if __name__ == "__main__":
    from AnyQt.QtWidgets import QApplication
    app = QApplication([])
    ow = OWGeneSets()
    ow.show()
    app.exec_()
