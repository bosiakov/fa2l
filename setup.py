from setuptools import setup
from os import path
from codecs import open


here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="fa2l",
    version="0.2.1",
    author="Eugene Bosiakov",
    author_email="eugenebosyakov@gmail.com",

    url='https://github.com/bosiakov/fa2l',

    description="Force Atlas 2 graph layout",
    long_description=long_description,

    packages=['fa2l'],

    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Mathematics',

        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    install_requires=[
        "networkx>=2.4,<3",
        "numpy"
    ],

    extras_require={
        'dev': ['matplotlib'],
    },
)
