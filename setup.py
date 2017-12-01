#!/usr/bin/env python

import os
from setuptools import setup, find_packages


NAME = 'Orange-Bioinformatics'
DOCUMENTATION_NAME = 'Orange Bioinformatics'

VERSION = '2.6.22'

DESCRIPTION = 'Orange Bioinformatics add-on for Orange data mining software package.'
LONG_DESCRIPTION = open(os.path.join(os.path.dirname(__file__), 'README.rst')).read()
AUTHOR = 'Bioinformatics Laboratory, FRI UL'
AUTHOR_EMAIL = 'contact@orange.biolab.si'
URL = 'http://orange.biolab.si/download'
LICENSE = 'GPLv3'

KEYWORDS = (
    'data mining',
    'machine learning',
    'artificial intelligence',
    'bioinformatics',
    'gene ontology',
    'KEGG',
    'expression profiles',
    'microarray',
    'genomics',
    'orange',
    'orange add-on',
    'orange3 add-on',
)

CLASSIFIERS = (
    'Development Status :: 4 - Beta',
    'Environment :: X11 Applications :: Qt',
    'Environment :: Console',
    'Environment :: Plugins',
    'Programming Language :: Python',
    'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering :: Artificial Intelligence',
    'Topic :: Scientific/Engineering :: Visualization',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Intended Audience :: Education',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
)

PACKAGES = find_packages()
PACKAGE_DATA = {}

requirements = ['requirements.txt']

INSTALL_REQUIRES = sorted(set(
    line.partition('#')[0].strip()
    for file in (os.path.join(os.path.dirname(__file__), file)
                 for file in requirements)
    for line in open(file)
) - {''})


EXTRAS_REQUIRE = {
    ':python_version<="3.4"': ["typing"],
}

ENTRY_POINTS = {
    'orange.addons': (
        'bio = orangecontrib.bio'
    ),
    'orange.widgets': (
        'Bioinformatics = orangecontrib.bio.widgets'
    ),
    'orange.canvas.help': (
        'html-index = orangecontrib.bio.widgets:WIDGET_HELP_PATH'
    )
}

NAMESPACE_PACAKGES = ["orangecontrib", "orangecontrib.bio"]

TEST_SUITE = "orangecontrib.bio.tests.suite"

if __name__ == '__main__':
    setup(
        name=NAME,
        version=VERSION,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        url=URL,
        license=LICENSE,
        keywords=KEYWORDS,
        classifiers=CLASSIFIERS,
        packages=PACKAGES,
        package_data=PACKAGE_DATA,
        install_requires=INSTALL_REQUIRES,
        extras_require=EXTRAS_REQUIRE,
        entry_points=ENTRY_POINTS,
        namespace_packages=NAMESPACE_PACAKGES,
        test_suite=TEST_SUITE,
        include_package_data=True,
        zip_safe=False,
    )
