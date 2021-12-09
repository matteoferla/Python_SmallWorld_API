from setuptools import setup, find_packages

# ----------- python version check
import sys

if sys.version_info.major != 3 or sys.version_info.minor < 6:
    print(sys.version_info)
    raise SystemError('Module written for Python 3.6+.')

# -------------- fill docstring
import os
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    __doc__ = f.read()

description='An (unofficial) Python3 module to query the SmallWorld chemical space search server (https://sw.docking.org/search.html)'

setup(
    name='SmallWorld_API',
    version='0.1',
    packages=find_packages(),
    install_requires=['pandas'],  # rdkit is optional.
    url='https://github.com/matteoferla/Python_SmallWorld_API',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo@well.ox.ac.uk',
    description=description,
    long_description=__doc__,
    long_description_content_type='text/markdown'
)