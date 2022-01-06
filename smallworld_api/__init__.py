__all__ = ['SmallWorld', 'NoMatchError']

import operator, re, json
from warnings import warn
import pandas as pd
import requests
from .defaults import Defaults  # class attributes
from .base import Base  # inherits Defaults
from .extras import Extras  # extra methods not required by search


class NoMatchError(Exception):
    """
    No match for the molecule was found.
    """

    def __str__(self):
        return 'The API returned no matches. ' + \
               '. '.join(map(str, self.args)) + \
               ' Try a different database (cf. `SmallWorld.retrieve_databases()`'


class SmallWorld(Extras):
    """
    A python3 API based upon https://wiki.docking.org/index.php/How_to_use_SmallWorld_API
    """

    # class attributes are in Defaults.

    def search(self,
               smiles: str,
               db: str = 'REAL_DB_20Q2',
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

        Calls ``submit_query`` and then ``get_results``.
        Returns a pandas dataframe of results.
        The dataframe is not rdkit modified yet.
        """
        start = int(other_parameters['start']) if 'start' in other_parameters else 0
        dist = int(other_parameters['dist']) if 'start' in other_parameters else 0
        length = int(other_parameters['length']) if 'length' in other_parameters else 10
        draw = int(other_parameters['draw']) if 'draw' in other_parameters else 10
        valids = {k: other_parameters[k] for k in self.valid_submit_keys if k in other_parameters}
        if db not in self.db_choices:
            warn(f'{db} is not a valid choice ({self.db_choices}).' + \
                 'Check updated with `.retrieve_scorefun_options()`')
        params = {'smi':  smiles,
                  'db':   db,
                  **self.default_submission,
                  'dist': int(dist),
                  **valids}
        self.query_summary = self.submit_query(params)
        self.hit_list_id = self.query_summary['hlid']
        return self.get_results(start, length, draw)

    def submit_query(self, params):
        """
        The first step.
        """
        reply: requests.Response = self._retrieve(url='/search/submit', params=params)
        line_iter = map(bytes.decode, reply.iter_lines())
        line_iter = map(str.strip, line_iter)
        line_iter = filter(lambda line: re.search(r'data:', line), line_iter)
        # using iter_lines + in stream mode does not solve the major tom hanging issue...
        reply_data = [json.loads(re.sub(r'^data:\s?', '', line)) for line in line_iter]
        if 'hlid' not in reply_data[-1]:
            raise ValueError(reply_data[-1])
        if reply_data[-1]['status'] != 'END':
            # "Ground Control to Major Tom" means there is no signal.
            warn(f"No completed return code returned: {reply_data} (generally harmless)")
        return reply_data[-1]

    def get_results(self, start: int = 0, length: int = 10, draw: int = 10) -> pd.DataFrame:
        params = dict(hlid=self.hit_list_id,
                      start=start,
                      length=length,
                      draw=draw
                      )
        params = {**params, **self.valid_export_columns}
        reply = self._retrieve(url='/search/view',
                               params=params)
        if not reply.json()["recordsTotal"]:
            raise NoMatchError(f'There are {reply.json()["recordsTotal"]} hits in the reply')
        reply_data = reply.json()['data']
        if not reply_data:
            raise NoMatchError('There is no `data` in the reply!')
        # expand the first column
        df1 = pd.DataFrame(map(operator.itemgetter(0), reply_data))
        if len(df1) == 0:
            raise NoMatchError('Reply generated an empty table')
        df2 = pd.DataFrame(reply_data).drop(columns=[0])
        columns = [v for p, v in params.items() if re.match(r'columns\[\d+]\[name]', p)]
        df2.columns = columns[1:]
        df = pd.concat([df1, df2], axis=1)
        df['smiles'] = df.hitSmiles.str.split(expand=True)[0]
        # PandasTools.AddMoleculeColumnToFrame(df,'smiles','molecule',includeFingerprints=True)
        return df
