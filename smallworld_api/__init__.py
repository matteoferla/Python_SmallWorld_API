from __future__ import annotations
"""
The primary class in the smallworld_api module is ``SmallWorld``.
All other classes, in other files are inherited by this.
"""

__all__ = ['SmallWorld', 'NoMatchError']

from warnings import warn
import pandas as pd
from .defaults import Defaults  # class attributes
from .base import Base  # inherits Defaults
from .extras import Extras  # extra methods not required by search
from typing import *
from .nomatcherror import NoMatchError
from .search import Searcher
import time, sys

if TYPE_CHECKING or 'sphinx' in sys.modules:
    from rdkit import Chem


class SmallWorld(Searcher):

    """
    A python3 API based upon https://wiki.docking.org/index.php/How_to_use_SmallWorld_API

    This (``SmallWorld``) is the main class of the module ``smallworld_api`` and is split into separate files
    by functionality. The classes it bases are as follows:
    Defaults -> Common -> Base -> Extras -> Searcher -> SmallWorld

    For example, class attributes are in Defaults.
    """

    def search_smiles(self,
                      smiles: str,
                      db: str,
                      **other_parameters) -> pd.DataFrame:
        """
        Given a smiles and a database return the table of results!

        The optional arguments are:

        * dist = 10 (atom difference distance threshhold)
        * several in `.default_submission`...

        The number of results given are controlled by:

        * length = 10 (number of results)
        * draw = 10 (pointless atm)
        * start = 0

        Which are passed onto `.get_results`.

        Calls ``submit_query`` and then ``get_results`` (code in ``search.py``)
        Returns a pandas dataframe of results.
        The dataframe is not rdkit modified yet.
        """
        start = int(other_parameters['start']) if 'start' in other_parameters else 0
        dist = int(other_parameters['dist']) if 'start' in other_parameters else 0
        length = int(other_parameters['length']) if 'length' in other_parameters else 10
        draw = int(other_parameters['draw']) if 'draw' in other_parameters else 10
        valids = {k: other_parameters[k] for k in self.valid_submit_keys if k in other_parameters}
        if db not in self.db_choices:
            warn(f'{db} is not a valid choice ({self.db_choices}).' +
                 'Check updated with `.retrieve_scorefun_options()`')
        params = {'smi':  smiles,
                  'db':   db,
                  **self.default_submission,
                  'dist': int(dist),
                  **valids}
        self.query_summary: Dict[str, Any] = self.submit_query(params)
        self.hit_list_id: int = self.query_summary['hlid']
        try:
            results = self.get_results(start, length, draw)
        except BaseException as error:
            warn(f'{error.__class__.__name__}: {error} was raised. '+
                 'Retrying in 2 seconds. There may be connection issues')
            time.sleep(2)
            results = self.get_results(start, length, draw)
        return results

    def search_mol(self,
                   mol: Chem.Mol,
                   db: str,
                   **other_parameters) -> pd.DataFrame:
        smiles = self.mol2smiles(mol)
        return self.search_smiles(smiles=smiles, db=db, **other_parameters)

    def search_many(self,
                    query: Union[Sequence[Any], Mapping[str, Any]],
                    db: str,
                    **other_parameters) -> pd.DataFrame:
        """
        search for many SMILES or Chem.Mol.
        """
        results: List[pd.DataFrame] = []
        if isinstance(query, Sequence):  # list or tuple etc.
            iterator = enumerate(query)
        elif isinstance(query, Mapping):  # dict etc.
            iterator = query.items()
        elif isinstance(query, pd.Series):  # a pd.Series is not a sequence or a mapping apparently
            iterator = query.to_dict().items()
        else:
            raise TypeError(f'Unrecognised type: {type(query)} for `.search_many_smiles`')
        tick = 0
        for name, item in iterator:
            # ## what is it?
            if isinstance(item, str):
                # its a smiles
                smiles = item
            elif self.is_this_mol(item):
                # rdkit
                smiles = self.mol2smiles(item)
            elif not item:
                # nothing
                warn(f'Falsy value {item} in the {type(query)} query')
                continue
            else:
                # mystery
                raise TypeError(f'Unrecognised type {type(item)}')
            # ## prevent excessive calls
            tock = time.time()
            if tick > tock - self.speed_threshold:
                time.sleep(self.speed_threshold - (tock - tick))
            tick = time.time()
            # ## what to capture
            if 'tolerate_all_exceptions' in other_parameters:
                tolerated_exceptions = (Exception,)
            elif 'tolerated_exceptions' in other_parameters:
                tolerated_exceptions = other_parameters['tolerated_exceptions']
            elif 'tolerate_NoMatchError' in other_parameters and other_parameters['tolerate_NoMatchError']:
                tolerated_exceptions = (NoMatchError,)
            else:
                tolerated_exceptions = ()
            # ## run!
            try:
                self.reset()
                result: pd.DataFrame = self.search_smiles(smiles=smiles, db=db, **other_parameters)
                result['query_index'] = name
                result['query_smiles'] = smiles
                results.append(result)
            except tolerated_exceptions as error:
                warn(f'{error.__class__.__name__}: {error} for {name}')
        return pd.concat(results, axis='index', ignore_index=True)

    def search(self, query: Any, db: str, **other_parameters) -> pd.DataFrame:
        """
        The query can be 
        
        * a single SMILES,
        * a rdkit.Chem.Mol
        * a list of SMILES or rdkit.Chem.Mol
        * a dictionary of SMILES or rdkit.Chem.Mol
        
        These all lead back to ``.search_smiles``, which functions as follows:
        """
        if isinstance(query, str):
            return self.search_smiles(smiles=query, db=db, **other_parameters)
        elif self.is_this_mol(query):  # rdkit is optional.
            self.search_mol(mol=query, db=db, **other_parameters)
        elif isinstance(query, Sized) and len(query) == 0:
            raise ValueError('Empty query')
        elif isinstance(query, Mapping) or isinstance(query, Sequence):
            return self.search_many(query=query, db=db, **other_parameters)
        else:
            raise TypeError(f'Unknown type {type(query)} for query')


SmallWorld.search.__doc__ += SmallWorld.search_smiles.__doc__
