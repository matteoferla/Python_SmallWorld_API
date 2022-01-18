class NoMatchError(Exception):
    """
    No match for the molecule was found.
    """

    def __str__(self):
        return 'The API returned no matches. ' + \
               '. '.join(map(str, self.args)) + \
               ' Try a different database (cf. `SmallWorld.retrieve_databases()`'
