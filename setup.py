from setuptools import setup, find_packages, Extension

setup(
    name = 'spkmeansmodule',
    version = '0.1.0',
    author = 'ME',
    description=  'kmeans',
    install_requires=['invoke'],
    packages= find_packages(),
    license= 'GPL-2',
    classifiers= ['Development Status :: 3- Alpha',
                  'License :: OSI Approved :: GNU General Public License v2 (GPLv2),'
                  'Natural Language: English',
                  'Programming Language :: Python :: Implementation :: CPython',],
    ext_modules= [
        Extension(
            'spkmeansmodule',
            ['spkmeans.c','spkmeansmodule.c'],
        ),
    ]
)