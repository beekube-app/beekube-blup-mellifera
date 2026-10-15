"""Mini-base Index Mellifera pour les tests des volets 2 et 3 (sous-ensemble du schéma de mellifera/db.py)."""
import sqlite3

from blup2.schemas import BlupV2Input, OptionsIn
from tests.builders import criterion

# Colonnes lues par blup2/external.py et blup2/external_records.py ; le contrat avec la vraie base du crawler est vérifié par
# tests/mellifera/test_mellifera_external_contract.py.
SCHEMA = """
CREATE TABLE queen (
    id INTEGER PRIMARY KEY AUTOINCREMENT, code TEXT NOT NULL UNIQUE, name TEXT, birth_year INTEGER,
    mother_code TEXT, drone_mother_code TEXT, n_drone_queens INTEGER, mating_code TEXT,
    single_drone INTEGER NOT NULL DEFAULT 0, race_description TEXT, breeder_id INTEGER);
CREATE INDEX queen_name ON queen (name, birth_year);
CREATE TABLE mating_code (code TEXT PRIMARY KEY, mating_type TEXT NOT NULL, mating_type_override TEXT);
CREATE VIEW queen_mating AS
    SELECT q.*, COALESCE(m.mating_type_override, m.mating_type, 'unknown') AS mating_type
    FROM queen q LEFT JOIN mating_code m ON m.code = q.mating_code;
INSERT INTO mating_code VALUES ('ins', 'insemination', NULL), ('nat', 'natural', NULL), ('xx', 'unknown', NULL);
"""

NOTES_SCHEMA = """
CREATE TABLE trait (id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, label TEXT NOT NULL, kind TEXT NOT NULL,
    scale_min REAL, scale_max REAL, direction INTEGER NOT NULL);
CREATE TABLE criterion (id INTEGER PRIMARY KEY AUTOINCREMENT, mode TEXT NOT NULL, code TEXT NOT NULL,
    trait_id INTEGER REFERENCES trait (id), UNIQUE (mode, code));
CREATE TABLE evaluation (
    id INTEGER PRIMARY KEY AUTOINCREMENT, queen_name_raw TEXT NOT NULL, queen_birth_year_raw INTEGER,
    queen_code TEXT, queen_key TEXT NOT NULL, page_breeder_id INTEGER, tester TEXT NOT NULL,
    tester_breeder_id INTEGER, eval_date TEXT NOT NULL, eval_year INTEGER, mode TEXT NOT NULL,
    criterion_code TEXT NOT NULL, value REAL, value_raw TEXT NOT NULL, source TEXT NOT NULL, url TEXT NOT NULL,
    UNIQUE (queen_key, source, mode, tester, eval_date, criterion_code));
INSERT INTO trait VALUES (1, 'swarming', 'Essaimage', 'note', 0, 6, 1), (2, 'mildness', 'Douceur', 'note', 0, 6, 1),
    (19, 'hyg24h', 'Couvain nettoyé en 24 h (%)', 'percentage', NULL, NULL, 1),
    (21, 'honeyAmt', 'Quantité de miel (kg)', 'harvest', NULL, NULL, 1);
INSERT INTO criterion (mode, code, trait_id) VALUES ('evjp', 'mildness', 2), ('evsp', 'mildness', 2),
    ('radar', 'mildness', 2), ('evjp', 'swarming', 1), ('evjp', 'hyg24h', 19), ('evjp', 'honeyAmt', 21),
    ('evjp', 'comment', NULL);
"""
BREEDERS = {"CBP": 145, "NG": 200}         # éleveur Index Mellifera d'après le trigramme du code

# Lignée CBP : 0276 ← mère 0100 (× drones 0009 NG) ← mère 0050 ← mère 0010 ; 0277 sœur de 0276.
LINEAGE = [
    # code, année, mère, « Drones », nom, code d'accouplement
    ("FR-CBP-B-0276-2022", 2022, "FR-CBP-B-0100-2020", "FR-NG-B-0009-2018", "B276(CBP)", "ins"),
    ("FR-CBP-B-0277-2022", 2022, "FR-CBP-B-0100-2020", "FR-NG-B-0009-2018", "B277(CBP)", "ins"),
    ("FR-CBP-B-0100-2020", 2020, "FR-CBP-B-0050-2018", None, "B100(CBP)", "nat"),
    ("FR-CBP-B-0050-2018", 2018, "FR-CBP-B-0010-2016", None, None, None),
    ("FR-CBP-B-0010-2016", 2016, None, None, None, None),
    ("FR-NG-B-0009-2018", 2018, None, None, "B9(NG)", "nat"),
    ("FR-ZZ-B-0001-2020", 2020, None, None, "TWIN", None),          # nom porté par deux reines
    ("FR-ZZ-B-0002-2020", 2020, None, None, "TWIN", None),
]


# Contemporaines : reines notées par les mêmes testeurs, les mêmes années, que des ancêtres de la lignée CBP.
PEERS = [
    ("FR-T-B-0001-2022", 2022, None, None, None, None),
    ("FR-T-B-0002-2022", 2022, None, None, None, None),
    ("FR-T-B-0003-2020", 2020, "FR-CBP-B-0050-2018", None, None, None),        # fille d'un ancêtre
    ("FR-CBP-B-0300-2023", 2023, None, None, None, None),                      # radar du même éleveur (145)
]

# Notes : code, source, mode, critère du site, valeur, date, testeur. Ancêtres 0276, 0100 et 0009 puis
# contemporaines ; 0009 (testeur 301 en 2019) n'a pas de contemporaine.
NOTES = [
    ("FR-CBP-B-0276-2022", "season_table", "evjp", "mildness", 5, "2023-07-06", 300),
    ("FR-CBP-B-0276-2022", "season_table", "evjp", "swarming", 3, "2023-07-06", 300),
    ("FR-CBP-B-0276-2022", "radar", "radar", "mildness", 4, "2023-01-01", None),     # couverte par la saison 2023
    ("FR-CBP-B-0276-2022", "radar", "radar", "mildness", 3, "2024-01-01", None),     # rucher : éleveur CBP (145)
    ("FR-CBP-B-0100-2020", "season_table", "evjp", "mildness", 6, "2021-07-01", 300),
    ("FR-CBP-B-0100-2020", "season_table", "evsp", "mildness", 4, "2021-07-01", 300),  # même clé : moyennée
    ("FR-NG-B-0009-2018", "season_table", "evjp", "mildness", 2, "2019-06-01", 301),
    ("FR-T-B-0001-2022", "season_table", "evjp", "mildness", 3, "2023-07-06", 300),
    ("FR-T-B-0001-2022", "season_table", "evjp", "swarming", 4, "2023-07-06", 300),
    ("FR-T-B-0002-2022", "season_table", "evjp", "mildness", 4, "2023-07-07", 300),
    ("FR-T-B-0003-2020", "season_table", "evjp", "mildness", 2, "2021-07-01", 300),
    ("FR-T-B-0003-2020", "season_table", "evjp", "mildness", 1, "2022-07-01", 300),  # autre année : ignorée
    ("FR-CBP-B-0300-2023", "radar", "radar", "mildness", 5, "2024-01-01", None),
    (None, "season_table", "evjp", "mildness", 4, "2023-07-08", 300),               # reine non résolue
]


def make_db(path, rows=LINEAGE, notes=None):
    """notes : None = base du volet 2 (sans tables de notes) ; sinon lignes au format de NOTES (PEERS ajoutées)."""
    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA)
    conn.executemany("INSERT INTO queen (code, birth_year, mother_code, drone_mother_code, name, mating_code) "
                     "VALUES (?, ?, ?, ?, ?, ?)", list(rows) + (PEERS if notes is not None else []))
    for trigram, breeder in BREEDERS.items():
        conn.execute("UPDATE queen SET breeder_id = ? WHERE code LIKE ?", (breeder, f"FR-{trigram}-%"))
    if notes is not None:
        conn.executescript(NOTES_SCHEMA)
        conn.executemany(
            "INSERT INTO evaluation (queen_name_raw, queen_code, queen_key, tester, tester_breeder_id, eval_date, "
            "mode, criterion_code, value, value_raw, source, url) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, '/x/')",
            [(code or "X", code, code or f"name:X|{k}", "" if tester is None else f"T{tester}", tester, date, mode,
              crit, value, "" if value is None else str(value), source)
             for k, (code, source, mode, crit, value, date, tester) in enumerate(notes)])
    conn.commit()
    conn.close()
    return str(path)


def inp(queens, depth=7, records=(), criteria=None, **kw):
    return BlupV2Input(exploitation=1, criteria=criteria or [criterion(1)], queens=queens, records=list(records),
                       options=OptionsIn(external_pedigree_depth=depth), **kw)


def warning(exp_or_out, code):
    return [w.message for w in exp_or_out.warnings if w.code == code]


CYCLES = [
    ("FR-C-B-0001-2020", 2020, "FR-C-B-0002-2019", "FR-C-B-0001-2020", None, "ins"),   # sa propre RPM
    ("FR-C-B-0002-2019", 2019, "FR-C-B-0003-2018", None, None, None),
    ("FR-C-B-0003-2018", 2018, "FR-C-B-0002-2019", None, None, None),                  # cycle 2 ↔ 3
    ("FR-C-B-0004-2017", 2017, "FR-C-B-0004-2017", None, None, None),                  # sa propre mère
]
