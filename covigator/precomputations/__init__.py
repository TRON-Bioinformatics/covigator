GENES_DICT = {
            'orf10': 'ORF10',
            's': 'S',
            'orf3a': 'ORF3a',
            'e': 'E',
            'm': 'M',
            'orf6': 'ORF6',
            'orf7a': 'ORF7a',
            'orf8': 'ORF8',
            'orf7b': 'ORF7b',
            'nc': 'N',  # tricky!
            # polyprotein...
            'nsp1': 'ORF1ab',
            'nsp10': 'ORF1ab',
            'nsp15': 'ORF1ab',
            'nsp2': 'ORF1ab',
            'nsp3': 'ORF1ab',
            'nsp4': 'ORF1ab',
            'nsp6': 'ORF1ab',
            'nsp7': 'ORF1ab',
            'nsp8': 'ORF1ab',
            'nsp9': 'ORF1ab'}

# transform protein coordinates in the polyprotein ORF1ab
# these values have been extracted from the CoVigator domains table derived from Pfam domains
# domains = pd.read_csv('../data/domain.csv')
POLYPROTEIN_OFFSET = {
    'nsp1': 9,
    'nsp10': 4262,
    'nsp15': 5929,
    'nsp2': 183,
    'nsp3': 880,
    'nsp4': 2788,
    'nsp6': 3598,
    'nsp7': 3860,
    'nsp8': 3943,
    'nsp9': 4141}