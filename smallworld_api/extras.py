import requests
import operator, re, json
from warnings import warn
import pandas as pd
import requests
from typing import *
from IPython.display import display, HTML
from .base import Base  # inherits Defaults
from collections import namedtuple

class Extras(Base):

    def show_reply_as_html(self, reply: Optional[requests.Response] = None):
        """
        The API calls may fail for some reason.
        Generally code 500 due to the server timing out and having issues.
        This prints the reply.
        """
        if reply is None:
            reply = self.last_reply
        display(HTML(reply.text))

    @classmethod
    def retrieve_scorefun_options(cls) -> pd.DataFrame:
        reply: requests.Response = requests.get('https://sw.docking.org/search/config')
        reply.raise_for_status()
        scores = pd.DataFrame(reply.json()['ScoreFuncs'])
        cls.sf_choices = scores.name.to_list()
        return scores

    @classmethod
    def retrieve_databases(cls) -> pd.DataFrame:
        reply = requests.get('https://sw.docking.org/search/maps')
        reply.raise_for_status()
        dbs = (pd.DataFrame.from_dict(reply.json(), orient='index')
               [['name', 'numEntries', 'numMapped', 'numUnmapped', 'numSkipped', 'status']]
               .sort_values('numMapped', ascending=False))
        cls.db_choices = dbs.name.to_list()
        return dbs

    @staticmethod
    def check_smiles(smiles: str): # -> rdkit.Chem.Mol
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None
        return mol
