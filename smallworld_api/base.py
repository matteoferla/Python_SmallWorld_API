import requests
from .common_db import Common
from typing import *


class Base(Common):  # Defaults -> Common -> Base -> Extras -> SmallWorld

    def __init__(self):
        """
        The two attributes added during initialisation are for debugging.
        """
        self.last_reply: requests.Response()  # debugging
        self.hit_list_id = -1
        self.session = requests.Session()
        self.query_summary: Dict[str, Any] = {}

    def reset(self):
        self.__init__()

    def _retrieve(self, url: str, params: Dict[str, Any]) -> requests.Response:
        """
        This is convoluted for debugging.
        If something goes wrong... ``.last_reply`` can be inspected.
        """
        self.last_reply: requests.Response = self.session.get(url=self.base_url + url,
                                                              params=params,
                                                              stream=self.stream_response,
                                                              timeout=600,  # 10 minutes
                                                              )
        self.last_reply.raise_for_status()
        return self.last_reply
