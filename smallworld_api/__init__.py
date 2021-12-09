__all__ = ['SmallWorld']

import operator, re, json
from warnings import warn
import pandas as pd
import requests
from .defaults import Defaults  # class attributes
from .base import Base  # inherits Defaults
from .extras import Extras  # extra methods not required by search


class SmallWorld(Extras):
    """
    A python3 API based upon https://wiki.docking.org/index.php/How_to_use_SmallWorld_API
    """

    # class attributes are in Defaults.

    def search(self,
               smiles: str,
               db: str = 'REAL_DB_20Q2',
               dist: int = 10,
               start: int = 0,
               length: int = 10,
               draw: int = 10,
               **other_parameters) -> pd.DataFrame:
        """
        Given a smiles and a database return the table of results!

        Calls ``submit_query`` and then ``get_results``.
        Returns a pandas dataframe of results.
        The dataframe is not rdkit modified yet.
        """
        valids = {k: other_parameters[k] for k in self.valid_submit_keys if k in other_parameters}
        assert db in self.db_choices, f'{db} is not a valid choice ({self.db_choices}).' + \
                                      'Check updated with `.retrieve_scorefun_options()'
        params = {'smi':  smiles,
                  'db':   db,
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
            warn(f"No completed return code returned: {reply_data}")  # "Ground Control to Major Tom"? Groan.
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
        columns = [v for p, v in params.items() if re.match(r'columns\[\d+]\[name]', p)]
        df1 = pd.DataFrame(map(operator.itemgetter(0), reply.json()['data']))
        df2 = pd.DataFrame(reply.json()['data']).drop(columns=[0])
        df2.columns = columns[1:]
        df = pd.concat([df1, df2], axis=1)
        df['smiles'] = df.hitSmiles.str.split(expand=True)[0]
        # PandasTools.AddMoleculeColumnToFrame(df,'smiles','molecule',includeFingerprints=True)
        return df
