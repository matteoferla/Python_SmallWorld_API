import re
from .defaults import Defaults


class Common(Defaults):  # Defaults -> Common -> Base -> Extras -> Searcher -> SmallWorld

    def _get_year_from_name(self, name) -> float:
        """
        Given a string containing 20Q3 returns 20.5

        :param name:
        :return:
        """
        rex = re.search(r'(?P<year>\d\d)Q(?P<quarter>\d)', name)
        if rex:
            return int(rex.group('year')) + int(rex.group('quarter')) / 4 - 0.25
        rex = re.search(r'20(?P<year>\d\d)', name)  # > 2099 --> sorry captain Kirk your dataset is too new
        if rex:
            return int(rex.group('year'))
        rex = re.search(r'(?P<year>\d\d)', name)
        if rex:  # dodgy?
            return int(rex.group('year'))
        return 0.

    @property
    def REAL_dataset(self) -> str:
        # As per PEP 8 REAL_dataset is better than real_dataset. So ignore PyCharm
        # likewise... This is not an enum.
        return self._latest_dataset(name='real', error_name='Enamine Real')

    @property
    def ZINC_dataset(self) -> str:
        # As per PEP 8 REAL_dataset is better than real_dataset. So ignore PyCharm
        # likewise... This is not an enum.
        return self._latest_dataset(name='zinc', error_name='Zinc')

    def _latest_dataset(self, name: str, error_name: str) -> str:
        options = [db for db in self.db_choices if name in db.lower()]
        if len(options) == 0:
            raise ValueError(f'There is no {error_name} in the options')
        return sorted(options, key=self._get_year_from_name, reverse=True)[0]
