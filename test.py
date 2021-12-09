import unittest
from smallworld_api import SmallWorld
import pandas as pd

class TestAPI(unittest.TestCase):
    def test_small_db(self):
        from IPython.display import display
        # input smiles
        aspirin = 'O=C(C)Oc1ccccc1C(=O)O'
        SmallWorld.check_smiles(aspirin)   # assert raised if gibberish
        # run!
        sws = SmallWorld()
        results: pd.DataFrame = sws.search(aspirin, dist=5, db='WorldDrugs-20Q2-3004', length=10)
        self.assertEqual(len(results), 10 )

    def test_big_db(self):
        from IPython.display import display
        # input smiles
        melatonin = 'COc1ccc2[nH]cc(CCNC(C)=O)c2c1'
        SmallWorld.check_smiles(melatonin)  # assert raised if gibberish
        # run!
        sws = SmallWorld()
        results: pd.DataFrame = sws.search(melatonin, dist=5, db='REAL_DB_20Q2', length=10)
        self.assertEqual(len(results), 10)

if __name__ == '__main__':
    unittest.main()
