# Python_SmallWorld_API

> This is Unofficial. So please do not abuse it or use it when you cannot legally use the site!

An (unofficial) Python3 module to query the SmallWorld chemical space search server (https://sw.docking.org/search.html)
.

[![Documentation Status](https://readthedocs.org/projects/python-smallworld-api/badge/?version=latest)](https://python-smallworld-api.readthedocs.io/en/latest/?badge=latest)
[![https img shields io pypi v smallworld api logo python](https://img.shields.io/pypi/v/smallworld--api?logo=python)](https://pypi.org/project/smallworld--api)   [![https img shields io pypi pyversions smallworld api logo python](https://img.shields.io/pypi/pyversions/smallworld--api?logo=python)](https://pypi.org/project/smallworld--api)   [![https img shields io pypi wheel smallworld api logo python](https://img.shields.io/pypi/wheel/smallworld--api?logo=python)](https://pypi.org/project/smallworld--api)   [![https img shields io pypi format smallworld api logo python](https://img.shields.io/pypi/format/smallworld--api?logo=python)](https://pypi.org/project/smallworld--api)   [![https img shields io pypi status smallworld api logo python](https://img.shields.io/pypi/status/smallworld--api?logo=python)](https://pypi.org/project/smallworld--api)   [![https img shields io pypi dm smallworld api logo python](https://img.shields.io/pypi/dm/smallworld--api?logo=python)](https://pypi.org/project/smallworld--api)   [![https img shields io codeclimate maintainability matteoferla Python_SmallWorld_API logo codeclimate](https://img.shields.io/codeclimate/maintainability/matteoferla/Python_SmallWorld_API?logo=codeclimate)](https://codeclimate.com/github/matteoferla/Python_SmallWorld_API)   [![https img shields io codeclimate issues matteoferla Python_SmallWorld_API logo codeclimate](https://img.shields.io/codeclimate/issues/matteoferla/Python_SmallWorld_API?logo=codeclimate)](https://codeclimate.com/github/matteoferla/Python_SmallWorld_API)   [![https img shields io codeclimate tech debt matteoferla Python_SmallWorld_API logo codeclimate](https://img.shields.io/codeclimate/tech-debt/matteoferla/Python_SmallWorld_API?logo=codeclimate)](https://codeclimate.com/github/matteoferla/Python_SmallWorld_API)   [![https img shields io github forks matteoferla Python_SmallWorld_API label Fork style social logo github](https://img.shields.io/github/forks/matteoferla/Python_SmallWorld_API?label=Fork&style=social&logo=github)](https://github.com/matteoferla/Python_SmallWorld_API)   [![https img shields io github stars matteoferla Python_SmallWorld_API style social logo github](https://img.shields.io/github/stars/matteoferla/Python_SmallWorld_API?style=social&logo=github)](https://github.com/matteoferla/Python_SmallWorld_API)   [![https img shields io github watchers matteoferla Python_SmallWorld_API label Watch style social logo github](https://img.shields.io/github/watchers/matteoferla/Python_SmallWorld_API?label=Watch&style=social&logo=github)](https://github.com/matteoferla/Python_SmallWorld_API)   [![https img shields io github last commit matteoferla Python_SmallWorld_API logo github](https://img.shields.io/github/last-commit/matteoferla/Python_SmallWorld_API?logo=github)](https://github.com/matteoferla/Python_SmallWorld_API)   [![https img shields io github license matteoferla Python_SmallWorld_API logo github](https://img.shields.io/github/license/matteoferla/Python_SmallWorld_API?logo=github)](https://github.com/matteoferla/Python_SmallWorld_API/raw/master/LICENCE)   [![https img shields io github release date matteoferla Python_SmallWorld_API logo github](https://img.shields.io/github/release-date/matteoferla/Python_SmallWorld_API?logo=github)](https://github.com/matteoferla/Python_SmallWorld_API)   [![https img shields io github commit activity m matteoferla Python_SmallWorld_API logo github](https://img.shields.io/github/commit-activity/m/matteoferla/Python_SmallWorld_API?logo=github)](https://github.com/matteoferla/Python_SmallWorld_API)   [![https img shields io github issues matteoferla Python_SmallWorld_API logo github](https://img.shields.io/github/issues/matteoferla/Python_SmallWorld_API?logo=github)](https://github.com/matteoferla/Python_SmallWorld_API)   [![https img shields io github issues closed matteoferla Python_SmallWorld_API logo github](https://img.shields.io/github/issues-closed/matteoferla/Python_SmallWorld_API?logo=github)](https://github.com/matteoferla/Python_SmallWorld_API)

## Overview

The SmallWorld server, [sw.docking.org](https://sw.docking.org/search.html), allows one to search for similar compounds
to a give [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
in one of many databases —a very complex feat.

The API points of the site are described in
[wiki.docking.org/index.php/How_to_use_SmallWorld_API](https://wiki.docking.org/index.php/How_to_use_SmallWorld_API).

This Python3 module allows one to search it.

## Install

    pip install smallworld-api

## Usage

```python
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd  # for typehinting below

from smallworld_api import SmallWorld

aspirin = 'O=C(C)Oc1ccccc1C(=O)O'
sw = SmallWorld()
results : pd.DataFrame = sw.search(aspirin, dist=5, db=sw.REAL_dataset)

from IPython.display import display
display(results)
```

The first two import lines are optional. The code works without rdkit, but if pandas gets imported before PandasTools
and Chem imported not in _main_ then display issues happen.

So it's up to you to remember to run:

```python
PandasTools.AddMoleculeColumnToFrame(results, 'smiles', 'molecule', includeFingerprints=True)
```

The argument `db` for `.search` is a string and is the name of the database. These do seem to change, so to update the
class and choose one call:

```python
SmallWorld.retrieve_databases()  #: pd.DataFrame (.db_choices gets updated too)
```

The dynamic properties `.REAL_dataset` and `.ZINC_dataset` simply return the best value from the presets, which may have
become out of date (unless updated).

## Query terms

The first argument passed to `.search` can be:

* a `str` (SMILES)
* a `Chem.Mol` (rdkit is an optional requirement though)
* a list-like (sequence) or a dict-like (mapping) of the above, where the index or key becomes the name in the output
  table.

If one is sure that the correct dataset is being used and any raised `NoMatchError` is due to the SMILES, then once can
add for the last case the argument `tolerate_tolerate_NoMatchError=True`, which makes them ignored bar for a warning.

## Debug

The instantiation is set up so for debugging, namely it has two attributes of interest:

* `sw.last_reply`, a `requests.Response` instance
* `sw.hit_list_id` an integer representing the search (AKA. `hlid` in the server responses)

The errors raised are generally either `requests.HTTPError`
or `smallworld.NoMatchError`. The former is raised by a `requests.Response.raise_for_status` call and means there is a
status code that isn't 200, the latter is raised by one of the various checks in `sw.get_results()`.

For the former errors, i.e. those by a serverside-declared HTML-formatting error (eg. status code 404), if one is in a
Jupyter notebook one can do `sw.show_reply_as_html()`. Generally if you get status code 500, it is best to try again
tomorrow as the server is having a hard time and is probably not okay on the web.

For the latter, the result in `.last_reply` should be a JSON string, therefore should give something like this:

```python3
reply_data: dict = sw.last_reply.json()
```

A common issue is the change in database names, therefore do do and pick a different one
(ATM, the index of the dataframe is the name to use, but in 2021 it was the `name`)

```python3
from IPython.display import display

from smallworld_api import SmallWorld
db_table : pd.DataFrame = SmallWorld.retrieve_databases()
display( db_table )
```

There will be a "ground control to major Tom" warning in the first query. This weird reply means that the stream has
finished, but not closed or something. Ignore it.

Also, as a shorthand, `mol: Chem.Mol = SmallWorld.check_smiles(aspirin)`
can be called to check if the molecules is fine.

## Choices

The database choices can seen with the preset list `SmallWorld.db_choices`. But also this can be recached via the
classmethod `SmallWorld.retrieve_databases()`.

Two databases, `REAL_Space_21Q3_2B(public)` and `REAL_DB_20Q2`, are Enamine REAL databases
(aka. Enamine will make the compound on request). Previously, the
repository, [enamine-real-search](https://github.com/xchem/enamine-real-search) was good for this, but unfortunately
Enamine changed their endpoints. So I wrote this to take its place!
Despite the smaller number of entries, `REAL_DB_20Q2` gives the most hits and is less likely to "Major Tom out".

Likewise, the attribute `SmallWorld.sf_choices` (type list) and 
the classmethod `SmallWorld.retrieve_scorefun_options()` do the same.
The values are less and are: `['Atom Alignment', 'SMARTS Alignment', 'ECFP4', 'Daylight']`, but these
are activated by default and will be visible as columns in the resulting dataframe from a search call.

Here is the full list of databases:

```python
import pandas as pd

choices: pd.DataFrame = SmallWorld.retrieve_databases()

display(choices)
```
Which will return (as of writing on the 9th Dec 2021):

|                                               | name                       |   numEntries |   numMapped |   numUnmapped |   numSkipped | status    |
|:----------------------------------------------|:---------------------------|-------------:|------------:|--------------:|-------------:|:----------|
| REAL_Space_21Q3_All_2B_public.smi.anon        | REAL_Space_21Q3_2B(public) |   1950356098 |  1935062471 |      15293627 |            0 | Available |
| ZINC-All-2020Q2-1.46B.anon                    | ZINC-All-20Q2-1.46B        |   1468554638 |  1467030947 |       1523691 |          231 | Available |
| ZINC-For-Sale-2020Q2-1.46B.anon               | ZINC-For-Sale-20Q2-1.46B   |   1464949146 |  1463519428 |       1429718 |           22 | Available |
| ZINC20-ForSale-21Q3.smi.anon                  | ZINC20-ForSale-21Q3-1.4B   |   1479284919 |  1440784765 |      38500154 |           29 | Available |
| Enamine_REAL_Public_July_2020_Q1-2_1.36B.anon | REAL_DB_20Q2               |   1361198468 |  1350462346 |      10736122 |            0 | Available |
| Wait-OK-2020Q2-1.2B.anon                      | Wait-OK-20Q2-1.2B          |   1174063221 |  1172785190 |       1278031 |            1 | Available |
| WuXi-20Q4.smi.anon                            | WuXi-20Q4-600M             |   2353582875 |   600762581 |    1752820294 |          284 | Available |
| MculeUltimate-20Q2.smi.anon                   | MculeUltimate_20Q2_126M    |    126471523 |   126471523 |             0 |            0 | Available |
| WuXi-2020Q2-120M.anon                         | WuXi-20Q2-120M             |    339132361 |   120400570 |     218731791 |            0 | Available |
| mcule_ultimate_200407_c8bxI4.anon             | Mcule_ultimate_20Q2-126M   |    126471523 |    45589462 |      80882061 |            0 | Available |
| BB-All-2020Q2-26.7M.anon                      | BB-All-20Q2-26.7M          |     26787985 |    26707241 |         80744 |           16 | Available |
| In-Stock-2020Q2-13.8M.anon                    | In-Stock-20Q2-13.8M        |     13842485 |    13829086 |         13399 |            1 | Available |
| ZINC20-InStock-21Q3.smi.anon                  | ZINC20-InStock-21Q3-11M    |     11122445 |    11103910 |         18535 |            5 | Available |
| BBall.smi.anon                                | BB-All-21Q4-3.3M           |      3319960 |     3319705 |           255 |            6 | Available |
| BBnow.smi.anon                                | BB-Now-21Q4-2M             |      2076639 |     2076464 |           175 |            6 | Available |
| BB-Now-2020Q2-1.6M.anon                       | BB-Now-20Q2-1.6M           |      1649789 |     1649386 |           403 |            4 | Available |
| BB_50.smi.anon                                | BB-50-21Q4-1.5M            |      1483551 |     1483434 |           117 |            2 | Available |
| BB_10.smi.anon                                | BB-10-21Q4-1.2M            |      1243321 |     1243241 |            80 |            0 | Available |
| BB_40.smi.anon                                | BB-40-21Q4-590K            |       589959 |      589911 |            48 |            4 | Available |
| interesting.smi.anon                          | ZINC-Interesting-20Q2-320K |       320845 |      320773 |            72 |            1 | Available |
| ZINC-Interesting-2020Q2-300K.anon             | ZINC-Interesting-20Q2-300K |       307854 |      300765 |          7089 |            1 | Available |
| TCNMP-2020Q2-31912.anon                       | TCNMP-20Q2-31912           |        37438 |       31912 |          5526 |            0 | Available |
| BB_30.smi.anon                                | BB-30-21Q4-3K              |         3129 |        3119 |            10 |            0 | Available |
| WorldDrugs-2020Q2-3004.anon                   | WorldDrugs-20Q2-3004       |         3004 |        3003 |             1 |            0 | Available |
| HMDB-2020Q2-584.anon                          | HMDB-20Q2-584              |          585 |         584 |             1 |            0 | Available |

## Names

There is a Python module called [smallworld](https://github.com/benmaier/smallworld),
which implements the small world algorithm.
This is not an API to the [sw.docking.org](https://sw.docking.org/search.html) site.

The blog of the [sw.docking.org](https://sw.docking.org/search.html) site mentions a pysmallworld.
There is no mention of this in Google so I am guessing it is for a future feature?
I however need to use this now as 
I need it as a publicly usable example workflow of [Fragmenstein](https://github.com/matteoferla/Fragmenstein).

Also, there is a great and wacky boardgame called [Small World](https://boardgamegeek.com/boardgame/40692/small-world),
with a curious/agonising dynamic which forces you to not be a collector.






