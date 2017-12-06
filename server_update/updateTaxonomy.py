""" Taxonomy update """
import tarfile
import bz2


from server_update import *
from orangecontrib.bio.ncbi.taxonomy.utils import TaxonomyDB
from orangecontrib.bio.ncbi.taxonomy.config import DOMAIN, FILENAME_ALL, FILENAME, TAXDUMP_URL
from server_update.tests.test_Taxonomy import TaxonomyTest


TITLE = "NCBI Taxonomy"
TAGS = ["NCBI", "taxonomy", "organism names", "essential"]
helper = SyncHelper(DOMAIN, TaxonomyTest)

if not http_last_modified(TAXDUMP_URL) > sf_last_modified(DOMAIN, FILENAME_ALL):
    # remove update folder
    helper.remove_update_folder()
    sys.exit(up_to_date)

domain_path = sf_local.localpath(DOMAIN)
temp_path = os.path.join(domain_path, sf_temp)
create_folder(domain_path)
create_folder(temp_path)

taxdump_filename = os.path.join(domain_path, "taxdump.tar.gz")
db_filename = os.path.join(domain_path, FILENAME_ALL)

TaxonomyDB.download(domain_path)
TaxonomyDB.init_db(db_filename, tarfile.open(taxdump_filename))
create_info_file(db_filename)  # to run tests, we need .info file -> bio.taxonomy.pickled_cache

db_size = os.stat(db_filename).st_size  # store uncompressed database size

with bz2.BZ2File(os.path.join(temp_path, FILENAME_ALL), mode="w", compresslevel=9) as f:
    shutil.copyfileobj(open(db_filename, "rb"), f)
create_info_file(os.path.join(temp_path, FILENAME_ALL), title=TITLE, tags=TAGS, uncompressed=db_size, compression='bz2')


# smaller taxonomy database, use only hierarchy from common tax ids
db_common_filename = os.path.join(domain_path, FILENAME)
tax_db = TaxonomyDB(db_filename)
common_organisms = tax_db.commonly_used_organisms(db_common_filename)

db_size = os.stat(db_common_filename).st_size  # store uncompressed database size

with bz2.BZ2File(os.path.join(temp_path, FILENAME), mode="w", compresslevel=9) as f:
    shutil.copyfileobj(open(db_common_filename, "rb"), f)
create_info_file(os.path.join(temp_path, FILENAME), title=TITLE, tags=TAGS,
                 uncompressed=db_size, compression='bz2')

# sync files with remote server
helper.run_tests()
helper.sync_files()

helper.remove_update_folder()
