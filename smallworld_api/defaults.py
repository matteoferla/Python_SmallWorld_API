class Defaults:  # Defaults -> Common -> Base -> Extras -> SmallWorld

    base_url = 'https://sw.docking.org'
    stream_response = False

    # database choice list is updated with ``df = SmallWorld.retrieve_databases()``
    db_choices = ['WuXi-20Q4.smi.anon', 'REAL_Space_21Q3_All_2B_public.smi.anon',
                  'all-zinc.smi.anon', 'wait-ok.smi.anon', 'MculeUltimate-20Q2.smi.anon',
                  'instock.smi.anon', 'BBall.smi.anon', 'BBnow.smi.anon',
                  'interesting.smi.anon']

    # most of these are enabled by default
    sf_choices = ['Atom Alignment', 'SMARTS Alignment', 'ECFP4', 'Daylight']

    # these are the keys that are valid according to the documentation for a submit request:
    valid_submit_keys = ['smi', 'db', 'dist', 'tdn', 'tup', 'rdn', 'rup', 'ldn', 'lup', 'scores']
    default_submission = {'dist': 8,
                          'tdn':  6,
                          'rdn':  6,
                          'rup':  2,
                          'ldn':  2,
                          'lup':  2,
                          'maj':  6,
                          'min':  6,
                          'sub':  6}

    # this is the required information (kind of)
    valid_export_columns = {'columns[0][data]':           '0',
                            'columns[0][name]':           'alignment',
                            'columns[0][searchable]':     'true',
                            'columns[0][orderable]':      'false',
                            'columns[0][search][value]':  '',
                            'columns[0][search][regex]':  'false',
                            'columns[1][data]':           '1',
                            'columns[1][name]':           'dist',
                            'columns[1][searchable]':     'true',
                            'columns[1][orderable]':      'true',
                            'columns[1][search][value]':  '0-12',
                            'columns[1][search][regex]':  'false',
                            'columns[2][data]':           '2',
                            'columns[2][name]':           'ecfp4',
                            'columns[2][searchable]':     'true',
                            'columns[2][orderable]':      'true',
                            'columns[2][search][value]':  '',
                            'columns[2][search][regex]':  'false',
                            'columns[3][data]':           '3',
                            'columns[3][name]':           'daylight',
                            'columns[3][searchable]':     'true',
                            'columns[3][orderable]':      'true',
                            'columns[3][search][value]':  '',
                            'columns[3][search][regex]':  'false',
                            'columns[4][data]':           '4',
                            'columns[4][name]':           'topodist',
                            'columns[4][searchable]':     'true',
                            'columns[4][orderable]':      'true',
                            'columns[4][search][value]':  '0-8',
                            'columns[4][search][regex]':  'false',
                            'columns[5][data]':           '5',
                            'columns[5][name]':           'mces',
                            'columns[5][searchable]':     'true',
                            'columns[5][orderable]':      'true',
                            'columns[5][search][value]':  '',
                            'columns[5][search][regex]':  'false',
                            'columns[6][data]':           '6',
                            'columns[6][name]':           'tdn',
                            'columns[6][searchable]':     'true',
                            'columns[6][orderable]':      'true',
                            'columns[6][search][value]':  '0-6',
                            'columns[6][search][regex]':  'false',
                            'columns[7][data]':           '7',
                            'columns[7][name]':           'tup',
                            'columns[7][searchable]':     'true',
                            'columns[7][orderable]':      'true',
                            'columns[7][search][value]':  '0-6',
                            'columns[7][search][regex]':  'false',
                            'columns[8][data]':           '8',
                            'columns[8][name]':           'rdn',
                            'columns[8][searchable]':     'true',
                            'columns[8][orderable]':      'true',
                            'columns[8][search][value]':  '0-6',
                            'columns[8][search][regex]':  'false',
                            'columns[9][data]':           '9',
                            'columns[9][name]':           'rup',
                            'columns[9][searchable]':     'true',
                            'columns[9][orderable]':      'true',
                            'columns[9][search][value]':  '0-2',
                            'columns[9][search][regex]':  'false',
                            'columns[10][data]':          '10',
                            'columns[10][name]':          'ldn',
                            'columns[10][searchable]':    'true',
                            'columns[10][orderable]':     'true',
                            'columns[10][search][value]': '0-2',
                            'columns[10][search][regex]': 'false',
                            'columns[11][data]':          '11',
                            'columns[11][name]':          'lup',
                            'columns[11][searchable]':    'true',
                            'columns[11][orderable]':     'true',
                            'columns[11][search][value]': '0-2',
                            'columns[11][search][regex]': 'false',
                            'columns[12][data]':          '12',
                            'columns[12][name]':          'mut',
                            'columns[12][searchable]':    'true',
                            'columns[12][orderable]':     'true',
                            'columns[12][search][value]': '',
                            'columns[12][search][regex]': 'false',
                            'columns[13][data]':          '13',
                            'columns[13][name]':          'maj',
                            'columns[13][searchable]':    'true',
                            'columns[13][orderable]':     'true',
                            'columns[13][search][value]': '0-6',
                            'columns[13][search][regex]': 'false',
                            'columns[14][data]':          '14',
                            'columns[14][name]':          'min',
                            'columns[14][searchable]':    'true',
                            'columns[14][orderable]':     'true',
                            'columns[14][search][value]': '0-6',
                            'columns[14][search][regex]': 'false',
                            'columns[15][data]':          '15',
                            'columns[15][name]':          'hyb',
                            'columns[15][searchable]':    'true',
                            'columns[15][orderable]':     'true',
                            'columns[15][search][value]': '0-6',
                            'columns[15][search][regex]': 'false',
                            'columns[16][data]':          '16',
                            'columns[16][name]':          'sub',
                            'columns[16][searchable]':    'true',
                            'columns[16][orderable]':     'true',
                            'columns[16][search][value]': '0-6',
                            'columns[16][search][regex]': 'false',
                            'order[0][column]':           '0',
                            'order[0][dir]':              'asc',
                            'search[value]':              '',
                            'search[regex]':              'false'}
