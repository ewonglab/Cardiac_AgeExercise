# ==== PUBLICATION WORKFLOW SUMMARY ====
# Purpose: Shared TF-family taxonomy (ordered, first-match-wins rules) imported by the snATAC co-binding scripts.
# Steps: 1) Define gene-name prefix / membership rules per family; 2) expose assign_family() for consistent labelling across analyses.
# Provenance: verbatim copy of TOBIAS_Cobinding/tf_families.py; only this header was added.
# =======================================
#!/usr/bin/env python3
"""
Shared TF family taxonomy for the TOBIAS cobinding project.

Single source of truth so family labels are identical across
jaccard_clustering.py, density_controlled_cobinding.py and
enriched_motifs_analysis.py.

Rules are ordered, FIRST MATCH WINS — specific gene-name families come before
broad prefix families, and the generic ZF/Zbtb catch-all sits near the end so it
doesn't swallow more specific zinc-finger factors (KLF/Sp, GLI/ZIC, EGR, etc.).

Mouse gene symbols are Title-case (Hoxa1, Foxc1, Klf4); matching uses exact
case-sensitive prefixes / membership as they appear in the data.
"""

# AP-1 genes are excluded upstream in every analysis; listed here only so the
# taxonomy never mislabels them if they ever slip through.
_AP1 = {'Fos', 'Fosb', 'Fosl1', 'Fosl2', 'FosB', 'c-Fos', 'FOSL1', 'FOSL2',
        'Jun', 'Junb', 'Jund', 'JunB', 'JunD', 'c-Jun', 'JUNB', 'JUND'}

# Homeodomain (non-HOX) prefixes — kept separate from HOX clusters
_HOMEOBOX = ('Dlx', 'Lhx', 'Vsx', 'Nkx', 'Alx', 'Vax', 'Emx', 'En1', 'En2',
             'Otp', 'Otx', 'Pou', 'Lmx', 'Prrx', 'Pitx', 'Msx', 'Six', 'Gsx',
             'Gbx', 'Hmx', 'Isl', 'Meox', 'Evx', 'Nobox', 'Mnx', 'Mixl',
             'Sebox', 'Esx', 'Dmbx', 'Uncx', 'Barx', 'Barhl', 'Cdx', 'Shox',
             'Arx', 'Rax', 'Phox', 'Dprx', 'Tlx', 'Lbx', 'Hlx', 'Noto', 'Pdx',
             'Cux', 'Onecut', 'Hnf1', 'Gsc', 'Hesx', 'Pknox', 'Meis', 'Pbx',
             'Tgif', 'Prep', 'Crx', 'Vent', 'Dbx',
             'Irx', 'Bsx', 'Homez', 'Isx', 'Obox', 'Nanog', 'Hdx',  # Iroquois + strays
             'Hoxa')  # Hoxa caught earlier

_BHLH = ('Myc', 'Max', 'Mxd', 'Mnt', 'Mga', 'Bhlh', 'Neurod', 'Neurog', 'Ascl',
         'Atoh', 'Hand', 'Tcf3', 'Tcf4', 'Tcf12', 'Myod', 'Myf', 'Myog', 'Twist',
         'Id1', 'Id2', 'Id3', 'Id4', 'Usf', 'Srebf', 'Mlx', 'Tfe3', 'Tfeb',
         'Mitf', 'Hif', 'Arnt', 'Ahr', 'Clock', 'Bmal', 'Npas', 'Sim', 'Epas',
         'Olig', 'Ptf1a', 'Tal', 'Lyl', 'Scl', 'Nhlh', 'Ferd',
         'Tfap4', 'Msc', 'Msgn', 'Tcfl5')

_NR = ('Nr0', 'Nr1', 'Nr2', 'Nr3', 'Nr4', 'Nr5', 'Nr6', 'Rxr', 'Ppar', 'Esrr',
       'Esr', 'Thr', 'Rar', 'Vdr', 'Pgr', 'Hnf4', 'Rora', 'Rorb', 'Rorc')

_BZIP = ('Cebp', 'Atf', 'Creb', 'Crebl', 'Crem', 'Batf', 'Nfil3', 'Dbp', 'Tef',
         'Hlf', 'Bach', 'Maf', 'Mafb', 'Maff', 'Mafg', 'Mafk', 'Nfe2', 'Xbp1',
         'Ddit3')

_ETS = ('Ets', 'Elf', 'Elk', 'Etv', 'Erg', 'Fli', 'Gabp', 'Spdef', 'Ehf',
        'Spib', 'Spic', 'Spi1', 'Erf', 'Fev', 'Et2')

_ZF = ('Zfp', 'Zbtb', 'Zscan', 'Zkscan', 'Zic', 'Gli', 'Glis', 'Sall', 'Snai',
       'Scrt', 'Insm', 'Egr', 'Wt1', 'Ctcf', 'Yy1', 'Yy2', 'Bcl11', 'Ikzf',
       'Prdm', 'Plag', 'Rreb', 'Mtf1', 'Osr', 'Bcl6',
       'Zeb', 'Ovol', 'Gfi', 'Vezf', 'Zik', 'Thap',   # EMT-ZF + other zinc fingers
       'Klf')  # Klf normally caught earlier

# Ordered, first-match-wins. (family_name, predicate)
TF_FAMILY_RULES = [
    ('HOX',              lambda n: n.startswith('Hox')),
    ('FOX',              lambda n: n.startswith('Fox')),
    ('Homeobox',         lambda n: n.startswith(_HOMEOBOX)),
    ('bZIP/CEBP-ATF',    lambda n: n.startswith(_BZIP)),
    ('bHLH',             lambda n: n.startswith(_BHLH)),
    ('Nuclear-receptor', lambda n: n.startswith(_NR) or n in ('Ar', 'Nr')),
    ('ETS',              lambda n: n.startswith(_ETS)),
    ('KLF/Sp',           lambda n: n.startswith(('Klf', 'Sp')) and not n.startswith('Spi')),
    ('GATA',             lambda n: n.startswith('Gata') or n.startswith('Trps') or n.startswith('Gatad')),
    ('SOX',              lambda n: n.startswith('Sox') or n == 'Sry'),
    ('TBX',              lambda n: n.startswith(('Tbx', 'Tbr', 'Eomes')) or n == 'T'),
    ('PAX',              lambda n: n.startswith('Pax')),
    ('NFI',              lambda n: n.startswith(('Nfia', 'Nfib', 'Nfic', 'Nfix'))),
    ('ARID',             lambda n: n.startswith('Arid')),
    ('DMRT',             lambda n: n.startswith('Dmrt')),
    ('HMG',              lambda n: n.startswith('Hmg') or n in ('Bbx', 'Cic', 'Tox')),
    ('HSF',              lambda n: n.startswith('Hsf')),
    ('EBF',              lambda n: n.startswith('Ebf')),
    ('SRF',              lambda n: n == 'Srf'),
    ('CXXC',             lambda n: n.startswith('Cxxc') or n == 'Kdm2b'),
    ('GCM',              lambda n: n.startswith('Gcm')),
    ('RUNX',             lambda n: n.startswith('Runx') or n.startswith('Cbf')),
    ('TEAD',             lambda n: n.startswith('Tead')),
    ('MEF2',             lambda n: n.startswith('Mef2')),
    ('NFAT',             lambda n: n.startswith('Nfat')),
    ('NF-kB',            lambda n: n in ('Rela', 'Relb', 'Rel', 'Nfkb1', 'Nfkb2')),
    ('STAT',             lambda n: n.startswith('Stat')),
    ('IRF',              lambda n: n.startswith('Irf')),
    ('NF-Y',             lambda n: n.startswith('Nfy')),
    ('Notch/HES',        lambda n: n.startswith(('Hes', 'Hey')) or n == 'Rbpj'),
    ('SMAD',             lambda n: n.startswith('Smad')),
    ('TCF/LEF',          lambda n: n in ('Tcf7', 'Tcf7l1', 'Tcf7l2', 'Lef1')),
    ('E2F',              lambda n: n.startswith('E2f') or n.startswith('Tfdp')),
    ('ZF/Zbtb',          lambda n: n.startswith(_ZF)),
    ('Other',            lambda n: True),
]


def assign_family(tf):
    """Return the family name for a TF symbol (first-match-wins)."""
    if tf in _AP1:
        return 'AP-1'
    for fam, fn in TF_FAMILY_RULES:
        try:
            if fn(tf):
                return fam
        except Exception:
            pass
    return 'Other'


# Canonical plotting order (most-relevant pioneer/dynamic families first)
FAMILY_ORDER = [
    'HOX', 'FOX', 'Homeobox', 'PAX', 'bZIP/CEBP-ATF', 'bHLH', 'Nuclear-receptor',
    'ETS', 'KLF/Sp', 'GATA', 'SOX', 'TBX', 'RUNX', 'TEAD', 'MEF2', 'SRF', 'NFAT',
    'NF-kB', 'STAT', 'IRF', 'NF-Y', 'NFI', 'Notch/HES', 'SMAD', 'TCF/LEF', 'E2F',
    'HMG', 'ARID', 'DMRT', 'HSF', 'EBF', 'CXXC', 'GCM', 'ZF/Zbtb', 'Other',
]

# Stable colour per family for consistent figures across scripts
FAMILY_COLORS = {
    'HOX':              '#FF9800',
    'FOX':              '#8E24AA',
    'Homeobox':         '#FFC107',
    'bZIP/CEBP-ATF':    '#795548',
    'bHLH':             '#3F51B5',
    'Nuclear-receptor': '#009688',
    'ETS':              '#4CAF50',
    'KLF/Sp':           '#E53935',
    'GATA':             '#00ACC1',
    'SOX':              '#7CB342',
    'TBX':              '#F06292',
    'RUNX':             '#5D4037',
    'TEAD':             '#26A69A',
    'MEF2':             '#AB47BC',
    'NFAT':             '#EC407A',
    'NF-kB':            '#6D4C41',
    'STAT':             '#42A5F5',
    'IRF':              '#1E88E5',
    'NF-Y':             '#8D6E63',
    'Notch/HES':        '#00BCD4',
    'SMAD':             '#9CCC65',
    'TCF/LEF':          '#FFB300',
    'E2F':              '#FF5722',
    'ZF/Zbtb':          '#9C27B0',
    'PAX':              '#FFB74D',
    'NFI':              '#4DB6AC',
    'ARID':             '#BA68C8',
    'DMRT':             '#A1887F',
    'HMG':              '#90A4AE',
    'HSF':              '#FF7043',
    'EBF':              '#9575CD',
    'SRF':              '#26C6DA',
    'CXXC':             '#78909C',
    'GCM':              '#C0CA33',
    'Other':            '#BDBDBD',
    'AP-1':             '#212121',
}
