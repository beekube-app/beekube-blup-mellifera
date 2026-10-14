"""Mini-base Index Mellifera pour les tests du volet 2 (sous-ensemble du schéma de mellifera/db.py)."""
import sqlite3

from blup2.schemas import BlupV2Input, OptionsIn
from tests.builders import criterion

# Colonnes lues par blup2/external.py ; le contrat avec la vraie base du crawler est vérifié par
# tests/mellifera/test_mellifera_external_contract.py.
SCHEMA = """
CREATE TABLE queen (
    id INTEGER PRIMARY KEY AUTOINCREMENT, code TEXT NOT NULL UNIQUE, name TEXT, birth_year INTEGER,
    mother_code TEXT, drone_mother_code TEXT, n_drone_queens INTEGER, mating_code TEXT,
    single_drone INTEGER NOT NULL DEFAULT 0, race_description TEXT);
CREATE INDEX queen_name ON queen (name, birth_year);
CREATE TABLE mating_code (code TEXT PRIMARY KEY, mating_type TEXT NOT NULL, mating_type_override TEXT);
CREATE VIEW queen_mating AS
    SELECT q.*, COALESCE(m.mating_type_override, m.mating_type, 'unknown') AS mating_type
    FROM queen q LEFT JOIN mating_code m ON m.code = q.mating_code;
INSERT INTO mating_code VALUES ('ins', 'insemination', NULL), ('nat', 'natural', NULL), ('xx', 'unknown', NULL);
"""

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


def make_db(path, rows=LINEAGE):
    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA)
    conn.executemany("INSERT INTO queen (code, birth_year, mother_code, drone_mother_code, name, mating_code) "
                     "VALUES (?, ?, ?, ?, ?, ?)", rows)
    conn.commit()
    conn.close()
    return str(path)


def inp(queens, depth=7, records=()):
    return BlupV2Input(exploitation=1, criteria=[criterion(1)], queens=queens, records=list(records),
                       options=OptionsIn(external_pedigree_depth=depth))


def warning(exp_or_out, code):
    return [w.message for w in exp_or_out.warnings if w.code == code]


CYCLES = [
    ("FR-C-B-0001-2020", 2020, "FR-C-B-0002-2019", "FR-C-B-0001-2020", None, "ins"),   # sa propre RPM
    ("FR-C-B-0002-2019", 2019, "FR-C-B-0003-2018", None, None, None),
    ("FR-C-B-0003-2018", 2018, "FR-C-B-0002-2019", None, None, None),                  # cycle 2 ↔ 3
    ("FR-C-B-0004-2017", 2017, "FR-C-B-0004-2017", None, None, None),                  # sa propre mère
]
