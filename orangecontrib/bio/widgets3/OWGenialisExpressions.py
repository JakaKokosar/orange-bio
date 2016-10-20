import math
import os
import requests_cache


from PyQt4.QtGui import QTreeWidget, QTreeWidgetItem, QLabel, QFont, QFrame, QLineEdit, QHeaderView
from Orange.widgets.widget import OWWidget
from Orange.widgets import widget, gui, settings
from Orange.data import Table
from PyQt4.QtCore import Qt
from requests.exceptions import HTTPError, ConnectionError

from orangecontrib.bio import resolwe
from ..utils import environ


# Creates line separator
def h_line():
    line = QFrame()
    line.setFrameShape(QFrame.HLine)
    line.setFrameShadow(QFrame.Sunken)
    return line


LABELS = []

#  Support cache with requests_cache module
cache_path = os.path.join(environ.buffer_dir, "resolwe")

try:
    os.makedirs(cache_path)
except OSError:
    pass

cache_file = os.path.join(cache_path, 'GenialisExpressions_cache')
requests_cache.install_cache(cache_name=cache_file, backend='sqlite')


TIMEPOINT_COLUMN = 6
REPLICATE_COLUMN = 7


def tfloat(s):
    try:
        return float(s)
    except ValueError:
        return None


class OWGenialisExpressions(OWWidget):

    name = "GenialisExpressions"
    description = "Genialis (expressions) widget"
    icon = "../widgets/icons/GenCloud.svg"
    want_main_area = True
    priority = 100

    inputs = []
    outputs = [("Data", Table)]

    username = settings.Setting("")
    password = settings.Setting("")
    store_credentials = settings.Setting([('', ''), ('', '')])
    selectedServer = settings.Setting(0)
    expType = settings.Setting(0)
    log2 = settings.Setting(False)
    joinreplicates = settings.Setting(False)
    genesAsColumns = settings.Setting(False)
    setTimeVariable = settings.Setting(False)

    def __init__(self):
        super().__init__()

        self.servers = [
            ('https://dictyexpress.research.bcm.edu', 'dictyExpress', 'genesis'),
            ('http://localhost:8000/', 'Resolwe', 'resolwe')
        ]
        self.search = ""
        self.server = None
        self.items = []
        self.lastSelected = None  # store last selected customTreeItems
        self.expTypes = [('Expression RPKUM (polyA)', 'exp'), ('Expression RPKM (polyA)', 'rpkmpolya'),
                         ('Expression RPKM', 'rpkm'), ('Expression RPKUM', 'rpkum'),
                         ('Read counts (polyA)', 'rc'), ('Read counts (raw)', 'rc_raw')]

        self.labels = []

        self.selectedProject = ''

        # pagination support
        self.currentPage = 0
        self.allPages = 0
        self.next_offset = None
        self.prev_offset = None
        self.active_offset = 0

        self.controlArea.setMaximumWidth(350)
        self.controlArea.setMinimumWidth(350)

        """ Login Section """

        box = gui.widgetBox(self.controlArea, 'Server:')

        gui.comboBox(box, self, "selectedServer",
                     items=[title for url, title, api in self.servers], callback=self.updateServer)

        self.namefield = gui.lineEdit(box, self, "username", "Username:",
                                      labelWidth=100,
                                      orientation='horizontal',
                                      callback=self.AuthChanged)

        self.passfield = gui.lineEdit(box, self, "password", "Password:",
                                      labelWidth=100,
                                      orientation='horizontal',
                                      callback=self.AuthChanged)

        self.passfield.setEchoMode(QLineEdit.Password)

        self.controlArea.layout().addWidget(h_line())

        """ Options Section """

        box = gui.widgetBox(self.controlArea, 'Expression type:')

        gui.comboBox(box, self, "expType", items=[label for label, name in self.expTypes], callback=None)

        self.replicates_chechBox = gui.checkBox(self.controlArea, self, "joinreplicates",
                     "Average replicates (use median)")
        self.replicates_chechBox.setToolTip('Averages identical experiments by using medians as values.')

        self.transpose_checkBox = gui.checkBox(self.controlArea, self, "genesAsColumns", "Genes as columns")
        self.transpose_checkBox.setToolTip('Transpose genes from rows to columns')

        self.transform_checkBox = gui.checkBox(self.controlArea, self, "log2",
                 "Logarithmic (base 2) transformation")
        self.transform_checkBox.setToolTip('Returns log2(value+1) for each value')

        self.controlArea.layout().addWidget(h_line())


        box = gui.widgetBox(self.controlArea, "Filters:")

        self.filter_box = QTreeWidget(box)
        self.filter_box.itemClicked.connect(self.set_annotations)
        self.filter_box.setHeaderHidden(True)
        box.layout().addWidget(self.filter_box)

        self.controlArea.layout().addWidget(h_line())

        self.commit_button = gui.button(self.controlArea, self, "Commit", callback=self.commit)
        self.handle_commit_button(False)

        gui.rubber(self.controlArea)

        """ MAIN AREA """

        label = QLabel("Available experiments:")
        myFont = QFont()
        myFont.setBold(True)
        label.setFont(myFont)
        self.mainArea.layout().addWidget(label)

        self.mainArea.layout().addWidget(h_line())

        gui.lineEdit(self.mainArea, self, "search", "Search:",
                     labelWidth=100,
                     orientation='horizontal',
                     callbackOnType=True,
                     callback=self.search_update)

        """ Experiment Section """

        self.experimentsWidget = QTreeWidget(alternatingRowColors=True,
                                             rootIsDecorated=False,
                                             uniformRowHeights=True,
                                             sortingEnabled=True)

        self.experimentsWidget.setUniformRowHeights(True)

        self.experimentsWidget.setSelectionMode(QTreeWidget.ExtendedSelection)

        self.experimentsWidget.setItemDelegateForColumn(
            0, gui.IndicatorItemDelegate(self, role=Qt.DisplayRole))

        self.experimentsWidget.selectionModel().selectionChanged.connect(
            self.on_selection_changed)

        #for i in range(len(self.headerLabels)):
            #self.experimentsWidget.resizeColumnToContents(i)

        self.mainArea.layout().addWidget(self.experimentsWidget)

        box = gui.widgetBox(self.mainArea, orientation=Qt.Horizontal)

        self.previous = gui.button(box, self, "Previous", callback=self.previous_page)
        self.info = gui.label(box, self, 'Pages: ' + str(self.currentPage) + ' / ' + str(self.allPages))
        self.next = gui.button(box, self, "Next",  callback=self.next_page)
        self.next.setDisabled(True)
        self.previous.setDisabled(True)

        gui.rubber(box)

        self.AuthSet()
        if self.username and self.password:
            self.connect()

    def reset(self):
        self.items = []
        self.lastSelected = None
        self.experimentsWidget.clear()
        self.handle_commit_button(False)
        self.warning(1)
        self.error(1)

    def on_selection_changed(self):
        self.handle_commit_button(True)

    def handle_commit_button(self, handle):
        self.commit_button.setEnabled(handle)

    def set_annotations(self, item):
        with self.progressBar(100) as progress:
            self.reset()
            self.selectedProject = item.data(0, Qt.UserRole)  # This is id of selected project
            try:
                exprs_anno, metas = self.server.get_annotations(self.selectedProject,
                                                                0,
                                                                progress.advance)
            except ConnectionError as e:
                self.error(1, e.args[0])
            except ValueError as e:
                self.warning(1, e.args[0])

            else:
                if exprs_anno:
                    self.load_tree_items(exprs_anno)
                    self.handle_pages(metas)

    def turn_page(self):
        self.reset()
        with self.progressBar(100) as progress:

            exprs_anno, metas = self.server.get_annotations(self.selectedProject,
                                                                self.active_offset,
                                                                progress.advance)
            if exprs_anno:
                self.load_tree_items(exprs_anno)
                self.handle_pages(metas)

    def previous_page(self):
        self.active_offset = self.prev_offset
        self.turn_page()

    def next_page(self):
        self.active_offset = self.next_offset
        self.turn_page()

    def pagination_reset(self):
        self.currentPage = 0
        self.allPages = 0
        self.next_offset = None
        self.prev_offset = None
        self.active_offset = 0
        self.handle_next_previous()
        self.label_update()

    def handle_next_previous(self):

        self.next.setDisabled(not self.next_offset)
        self.previous.setDisabled(not self.prev_offset)
        if self.prev_offset == 0:
            self.previous.setDisabled(False)

    def handle_pages(self, metas):
        self.currentPage = metas['currentPage']
        self.allPages = metas['allPages']
        self.next_offset = metas['next_offset']
        self.prev_offset = metas['prev_offset']

        self.handle_next_previous()
        self.label_update()

    def label_update(self):
        self.info.setText('Pages: ' + str(self.currentPage) + ' / ' + str(self.allPages))

    def add_filters(self, parent, elements):
        self.filter_box.clear()
        for text, children in elements.items():
            root = self.add_parent(parent, 0, text)
            for title, proj_id in children:
                self.add_child(root, 0, title, proj_id)

    def add_parent(self, parent, column, title):
        item = QTreeWidgetItem(parent, [title])
        item.setData(column, Qt.UserRole, "all")
        return item

    def add_child(self, parent, column, title, proj_id):
        item = QTreeWidgetItem(parent, [title])
        item.setData(column, Qt.UserRole, proj_id)
        return item

    def AuthSet(self):
        self.passfield.setDisabled(not self.username)

    def AuthChanged(self):
        self.AuthSet()
        self.connect()

    def set_view(self, server_type):
        if server_type == 'genesis':
            self.labels = [' ', 'Experiment', 'Growth', 'Genotype', 'Treatment',
                           'Strain', 'Time', 'Replicate', 'Date created']
        else:
            self.labels = [' ', 'Sample', 'Experiment type', 'Molecule',
                           'Genotype', 'Cell type', 'Organism', 'Source', 'Date created']

        self.experimentsWidget.setHeaderLabels(self.labels)

    def search_update(self):
        parts = self.search.split()
        for item in self.items:
            item.setHidden(not all(s in item for s in parts))

    def connect(self):
        self.reset()
        self.server = None
        self.filter_box.clear()
        server_url = self.servers[self.selectedServer][0]
        server_type = self.servers[self.selectedServer][2]
        self.set_view(server_type)
        self.pagination_reset()

        self.store_credentials[self.selectedServer] = (self.username, self.password)
        login = self.store_credentials[self.selectedServer]
        if not login[0] or not login[1]:
            login = (self.username, self.password)

        self.namefield.setText(login[0])
        self.passfield.setText(login[1])

        try:
            self.server = resolwe.connect(login[0], login[1], server_url, server_type)

        except (resolwe.ResolweAuthException, ConnectionError) as e:
            self.error(1, e.args[0])

        else:
            if server_type == 'genesis':
                self.add_filters(self.filter_box.invisibleRootItem(), self.server.get_filters())
            else:
                 self.load_tree_items(self.server.get_data())

    def updateServer(self):
        self.username = self.store_credentials[self.selectedServer][0]
        self.password = self.store_credentials[self.selectedServer][1]
        self.connect()

    def load_tree_items(self, exprs_anno):
        self.experimentsWidget.clear()

        for exprs_id, anno in exprs_anno.items():
            self.items.append(CustomTreeItem(self.experimentsWidget, exprs_id, anno, self.labels))

        for i in range(len(self.labels)):
            self.experimentsWidget.resizeColumnToContents(i)

    def commit(self):
        self.warning(1)
        self.error(1)

        with self.progressBar(100) as progress:
            selected = [item for item in self.experimentsWidget.selectedItems()]

            if self.lastSelected:
                [last.setData(0, Qt.DisplayRole, "") for last in self.lastSelected]

            self.lastSelected = selected
            ids = [i.exprs_id for i in selected]
            anno = [i.annotations for i in selected]

            transform = None
            if self.log2:
                transform = lambda x: math.log(x + 1.0, 2)

            try:
                downloaded_exprs = self.server.download_exprs_data(ids,
                                                                   self.expTypes[self.expType][1], progress.advance)
            except (HTTPError, ConnectionError) as e:
                self.error(1, e.args[0])
            except ValueError as e:
                self.warning(1, e.args[0])
            else:
                data = resolwe.to_orange_table(downloaded_exprs, anno, self.labels[1:-1], transform,
                                               self.genesAsColumns, self.joinreplicates, progress.advance)

                [sel.setData(0, Qt.DisplayRole, " ") for sel in selected]

                self.send("Data", data)


class CustomTreeItem(QTreeWidgetItem):

    def __init__(self, parent, exprs_id, anno, labels):

        super(CustomTreeItem, self).__init__(parent)
        self.parent = parent

        self.exprs_id = exprs_id
        self.annotations = anno
        self.set_rows(self.annotations, labels)  # set rows in QtreeWidget

    def __contains__(self, text):
        return any(text.upper() in str(self.text(i)).upper() for i in range(self.columnCount()))

    def __lt__(self, other):
        col = self.parent.sortColumn()
        if col in [TIMEPOINT_COLUMN, REPLICATE_COLUMN]:
            left = tfloat(self.text(col))
            right = tfloat(other.text(col))
            if isinstance(left, float) and isinstance(right, float):
                return left < right

        return QTreeWidgetItem.__lt__(self, other)

    def set_rows(self, annot, labels):
        for index, label in enumerate(labels):
            if index > 0:
                self.setText(index, str(annot[index-1]))
