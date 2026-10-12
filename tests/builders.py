"""Constructeurs d'objets d'entrée pour les tests."""
from __future__ import annotations

from blup2.schemas import CriterionIn, MatingIn, QueenIn, RecordIn


def mating(type="station", drone_queens=(), mother=None, n_rpm=None, n_drones=None, station=None,
           year=None, ext=None) -> MatingIn:
    return MatingIn(type=type, drone_queens=list(drone_queens), drone_queens_mother=mother,
                    n_drone_queens=n_rpm, n_drones=n_drones, station=station, year=year, external_ref=ext)


def queen(id, born=2020, mother=None, mating=None, **kw) -> QueenIn:
    return QueenIn(id=id, born=born, mother=mother, mating=mating, **kw)


def record(queen, value, criterion=1, date="2024-05-01", apiary=1, evaluator=1, **kw) -> RecordIn:
    return RecordIn(queen=queen, criterion=criterion, value=value, date=date, apiary=apiary,
                    evaluator=evaluator, **kw)


def criterion(id=1, **kw) -> CriterionIn:
    return CriterionIn(id=id, name=f"c{id}", **kw)


def beebreed_queens(S=8, D=12) -> list[QueenIn]:
    """Pedigree de docs/analyses/verify_paper.py::beebreed_pedigree au format v2.

    Sires de base S0a / S0b = stations 1 / 2 sans mère connue ; S1 = groupe de RPM filles de Q0_1
    (station 3) ; S2 = groupe de RPM filles de Q1_4 (station 4). Identifiants : Q0_k = k,
    Q1_k = 10 + k, Q2_k = 20 + k, Q3_k = 30 + k.
    """
    s0a = dict(type="station", station=1, year=2000, n_rpm=S, n_drones=D)
    s0b = dict(type="station", station=2, year=2000, n_rpm=S, n_drones=D)
    s1 = dict(type="station", station=3, year=2001, mother=1, n_rpm=S, n_drones=D)
    s2 = dict(type="station", station=4, year=2002, mother=14, n_rpm=S, n_drones=D)
    qs = [queen(1, 2000, mating=mating(**s0a)), queen(2, 2000, mating=mating(**s0a)),
          queen(3, 2000, mating=mating(**s0b)), queen(4, 2000, mating=mating(**s0b))]
    for k in range(1, 4):
        qs.append(queen(10 + k, 2001, mother=2, mating=mating(**s1)))
    for k in range(4, 7):
        qs.append(queen(10 + k, 2001, mother=3, mating=mating(**s1)))
    qs.append(queen(17, 2001, mother=4, mating=mating(**s1)))
    for k in range(1, 4):
        qs.append(queen(20 + k, 2002, mother=11, mating=mating(**s2)))
    for k in range(4, 7):
        qs.append(queen(20 + k, 2002, mother=15, mating=mating(**s2)))
    for k in range(1, 3):
        qs.append(queen(30 + k, 2003, mother=21, mating=mating(type="natural", n_drones=D)))
    return qs
