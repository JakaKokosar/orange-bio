""" update PPI """
import gzip


from server_update import *
from server_update.tests import test_BioGRID
from orangecontrib.bio import ppi


DOMAIN = 'PPI'
domain_path = sf_local.localpath(DOMAIN)
tmp_path = sf_local.localpath(DOMAIN, sf_temp)
create_folder(domain_path)

""" BioGrid """
ppi.BioGRID.download_data("https://thebiogrid.org/downloads/archives/Latest%20Release/BIOGRID-ALL-LATEST.tab2.zip")
# This download directory contains the most recent data release from the BioGRID. All files are named the same every
# month for consistency to allow for automated scripted downloads.

sfn = ppi.BioGRID.SERVER_FILE
filename = os.path.join(domain_path, sfn)

TITLE = 'BioGRID Protein interactions'
TAGS = ['protein interaction', 'BioGrid']
VERSION = ppi.BioGRID.VERSION

create_folder(tmp_path)
with gzip.GzipFile(os.path.join(tmp_path, sfn), "wb") as gz:
    gz.write(open(filename, "rb").read())

create_info_file(os.path.join(tmp_path, sfn), title=TITLE, tags=TAGS, version=VERSION,
                 compression='gz', uncompressed=file_size_bytes(filename))


helper = SyncHelper(DOMAIN, test_BioGRID.BioGRIDTest)
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
