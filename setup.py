from setuptools import setup, find_packages

from warnings import warn

sw_url = 'https://sw.docking.org/search.html'
warn(f'''DISCLAIMER: To use the SmallWorld API please make sure you can legally query a server with it,
for example {sw_url} by John Irwin.
This tool is not affiliated with the great folk at NextMove Software nor the Irwin lab,
but I am grateful for their smashing work and so should you so don't forget to cite them!''')

# ----------- python version check
import sys

if sys.version_info.major != 3 or sys.version_info.minor < 6:
    print(sys.version_info)
    raise SystemError('Module written for Python 3.6+.')

# -------------- fill docstring
import os

this_directory = os.path.abspath(os.path.dirname(__file__))
__doc__ = 'Smallworld API'
if os.path.exists(os.path.join(this_directory, 'README.md')):
    # there is no manifest.in file, so it could be missing
    with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        __doc__ = f.read()

description = f'An (unofficial) Python3 module to query the SmallWorld chemical space search server ({sw_url})'

setup(
    name='smallworld-api',
    version='1.1.3',
    python_requires='>=3.7',
    packages=find_packages(),
    install_requires=['pandas', 'requests', 'ipython'],  # rdkit is optional.
    url='https://github.com/matteoferla/Python_SmallWorld_API',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo@well.ox.ac.uk',
    classifiers=[  # https://pypi.org/classifiers/
        'Development Status :: 5 - Production/Stable',  # Development Status :: 4 - Beta
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description=description,
    long_description=__doc__,
    long_description_content_type='text/markdown'
)
