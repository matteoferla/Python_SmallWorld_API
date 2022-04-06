import json
import operator
import re
from typing import *
from warnings import warn
import pandas as pd
import requests
from .extras import Extras  # extra methods not required by search
from .nomatcherror import NoMatchError


class Searcher(Extras):  # Defaults -> Common -> Base -> Extras -> Searcher -> SmallWorld

    def submit_query(self, params) -> Dict[str, Any]:
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
            warn(f"No completed return code returned (generally harmless). See `.query_summary` for actual response")
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
        df['name'] = df.hitSmiles.str.split(expand=True)[1]
        df['smiles'] = df.hitSmiles.str.split(expand=True)[0]
        # PandasTools.AddMoleculeColumnToFrame(df,'smiles','molecule',includeFingerprints=True)
        return df
