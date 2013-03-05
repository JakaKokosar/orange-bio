"""
Caching framework for cached kegg api calls.

"""
import os
import UserDict
import sqlite3
import cPickle as pickle

from contextlib import closing

from datetime import datetime, date, timedelta
from . import conf


class Store(object):
    def __init__(self):
        self.timestamp = 0

    def open(self):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


class Sqlite3Store(Store, UserDict.DictMixin):
    def __init__(self, filename):
        self.filename = filename
        self.con = sqlite3.connect(filename)
        self.con.execute("""
            CREATE TABLE IF NOT EXISTS cache
                (key TEXT UNIQUE,
                 value TEXT
                )
        """)
        self.con.execute("""
            CREATE INDEX IF NOT EXISTS cache_index
            ON cache (key)
        """)
        self.con.commit()

    def __getitem__(self, key):
        cur = self.con.execute("""
            SELECT value
            FROM cache
            WHERE key=?
        """, (key,))
        r = cur.fetchall()

        if not r:
            raise KeyError(key)
        else:
            return pickle.loads(str(r[0][0]))

    def __setitem__(self, key, value):
        value = pickle.dumps(value)
        self.con.execute("""
            INSERT OR REPLACE INTO cache
            VALUES (?, ?)
        """, (key, value))
        self.con.commit()

    def __delitem__(self, key):
        self.con.execute("""
            DELETE FROM cache
            WHERE key=?
        """, (key,))
        self.con.commit()

    def keys(self):
        cur = self.con.execute("""
            SELECT key
            FROM cache
        """)
        return [str(r[0]) for r in cur.fetchall()]

    def close(self):
        pass


class DictStore(Store, UserDict.DictMixin):
    def __init__(self):
        Store.__init__(self)

    def close(self):
        pass


class cache_entry(object):
    def __init__(self, value, mtime=None, expires=None):
        self.value = value
        self.mtime = mtime
        self.expires = expires

_SESSION_START = datetime.now()


class cached_wrapper(object):
    """
    TODO: needs documentation
    """
    def __init__(self, function, instance, class_, cache_store,
                 last_modified=None):
        self.function = function
        self.instance = instance
        self.class_ = class_
        self.cache_store = cache_store
        self.last_modified = last_modified

    def has_key(self, key):
        with closing(self.cache_store()) as store:
            return key in store

    def key_from_args(self, args, kwargs=None):
        key = self.function.__name__ + repr(args)
        return key

    def invalidate_key(self, key):
        with closing(self.cache_store()) as store:
            del store[key]

    def last_modified_from_args(self, args, kwargs=None):
        key = self.key_from_args(args, kwargs)
        if self.instance is not None:
            self.instance.last_modified(args)

    def invalidate_args(self, args):
        return self.invalidate_key(self.key_from_args(args))

    def invalidate_all(self):
        prefix = self.key_from_args(()).rstrip(",)")
        with self.cache_store() as store:
            for key in store.keys():
                if key.startswith(prefix):
                    del store[key]

    def memoize(self, args, kwargs, value, timestamp=None):
        key = self.key_from_args(args, kwargs)
        if timestamp is None:
            timestamp = datetime.now()

        with closing(self.cache_store()) as store:
            store[key] = cache_entry(value, mtime=timestamp)

    def __call__(self, *args):
        key = self.key_from_args(args)
        with closing(self.cache_store()) as store:
            valid = True
            if key not in store:
                valid = False
            else:
                entry = store[key]
                rval = entry.value

                if not self.is_entry_valid(entry, args):
                    valid = False
            if not valid:
                rval = self.function(self.instance, *args)
                store[key] = cache_entry(rval, datetime.now(), None)

        return rval

    def min_timestamp(self, args):
        key = self.key_from_args(args)
        return datetime.fromtimestamp(0)

    def is_entry_valid(self, entry, args):
        # Need to check datetime first (it subclasses date)
        if isinstance(entry.mtime, datetime):
            mtime = entry.mtime
        elif isinstance(entry.mtime, date):
            mtime = datetime(entry.mtime.year, entry.mtime.month,
                             entry.mtime.day, 1, 1, 1)
        else:
            return False

        if self.min_timestamp(args) > mtime:
            return False

        last_modified = self.last_modified_from_args(args)

        if isinstance(last_modified, date):
            last_modified = datetime(last_modified.year, last_modified.month,
                                     last_modified.day, 1, 1, 1)
        elif isinstance(last_modified, basestring):
            # Could have different format
            mtime = mtime.strftime("%Y %m %d %H %M %S")

        elif last_modified is None:
            if conf.params["cache.invalidate"] == "always":
                return False
            elif conf.params["cache.invalidate"] == "session":
                last_modified = _SESSION_START
            elif conf.params["cache.invalidate"] == "daily":
                last_modified = datetime.now().replace(hour=0, minute=0,
                                                       second=0, microsecond=0)
            elif conf.params["cache.invalidate"] == "weekly":
                last_modified = datetime.now() - timedelta(7)
            else:  # ???
                pass
        return last_modified <= mtime


class cached_method(object):
    def __init__(self, function):
        self.function = function

    def __get__(self, instance, owner):
        if instance is not None:
            return cached_wrapper(self.function, instance, owner,
                                  self.get_cache_store(instance, owner))
        return self

    def get_cache_store(self, instance, owner):
        if hasattr(instance, "cache_store"):
            return instance.cache_store
        elif not hasattr("_cached_method_cache"):
            instance._cached_method_cache = DictStore()
        return instance._cached_method_cache


class bget_cached_method(cached_method):
    def __get__(self, instance, owner):
        if instance is not None:
            return cached_wrapper(self.function, instance, owner,
                                  self.get_cache_store(instance, owner),
                                  self.get_last_modified(instance, owner))
        return self

    def get_last_modified(self, instance, owner):
        if hasattr(instance, "last_modified"):
            return instance.last_modified


def touch_dir(path):
    path = os.path.expanduser(path)
    if not os.path.exists(path):
        os.makedirs(path)


def clear_cache():
    """Clear all locally cached KEGG data.
    """
    import glob
    path = conf.params["cache.path"]
    if os.path.realpath(path) != os.path.realpath(conf.kegg_dir):
        raise Exception("Non default cache path. Please remove the contents "
                        "of %r manually." % path)

    for cache_filename in glob.glob(os.path.join(path, "*.sqlite3")):
        os.remove(cache_filename)

    for ko_filename in glob.glob(os.path.join(path, "*.keg")):
        os.remove(ko_filename)

    for kgml_filename in glob.glob(os.path.join(path, "*.xml")):
        os.remove(kgml_filename)

    for png_filename in glob.glob(os.path.join(path, "*.png")):
        os.remove(png_filename)