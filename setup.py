import os
from setuptools import find_packages, setup

# determining the directory containing setup.py
setup_path = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(setup_path, 'README.md'), encoding='utf-8') as f:
    readme = f.read()

setup(
    # package information
    name='pygna',
    packages=find_packages(),
    version='3.1.5-dev',
    description='Geneset Network Analysis',
    long_description=readme,
    license='MIT',
    url='https://github.com/stracquadaniolab/pygna',
    keywords='Bioinformatics Network Statistics',

    #  author information
    author='Viola Fanfani, Giovanni Stracquadanio',
    author_email='v.fanfani@sms.ed.ac.uk',

    # installation info and requirements
    install_requires=[
        'pandas',
        'numpy',
        'scipy',
        'matplotlib',
        'pyyaml',
        'tables>=3.4.4',
        'seaborn>=0.9',
        'palettable',
        'networkx==2.3',
        'statsmodels',
        'argh'
    ],
    setup_requires=[],

    # test info and requirements
    test_suite='tests',
    tests_require=[],
    python_requires='>=3',

    # package deployment info
    include_package_data=True,
    zip_safe=False,

    # all tools have cli interface
    entry_points={
        'console_scripts': [
            'pygna=pygna.cli:main',
        ],
    },

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
