"""Pedigree abeille : nœuds reine (Q), sire (S) et groupe d'ouvrières (W).

Règles (spec §1) :
- mère de Q(x) = Q(mother(x)) ;
- père de Q(x) = nœud sire de l'accouplement de sa mère ;
- W(x) : parents Q(x) et le nœud sire de l'accouplement de x ;
- fécondation naturelle / inconnue : pas de nœud sire.
Auto-références (spec §5.1) : une reine qui est sa propre mère, figure dans ses propres
drone_queens ou est sa propre drone_queens_mother perd ce lien avant la détection de cycle
(avertissement self_reference_removed) ; les cycles de longueur ≥ 2 restent une erreur.
Nœud sire d'un accouplement contrôlé (station, insémination) :
- une seule RPM listée (drone_queens = [P]) : le sire est le nœud Q(P) lui-même (S = 1) ;
- sinon nœud S « groupe de RPM sœurs » : mère = drone_queens_mother (ou mère commune des RPM
  listées), S = len(drone_queens) sinon n_drone_queens sinon 1 ; les accouplements de même clé
  (external_ref, sinon (station, année, mère des RPM, RPM listées, S)) partagent le même nœud.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:  # pragma: no cover
    from .schemas import MatingIn, QueenIn, WarningOut


MAX_LISTED_IDS = 20   # identifiants cités dans un message ; au-delà, « … et N autres »


def listed_ids(ids: list[int]) -> str:
    """« 3, 5, 8 » ; au-delà de MAX_LISTED_IDS identifiants : « 1, 2, …, 20 … et N autres »."""
    head = ", ".join(str(i) for i in ids[:MAX_LISTED_IDS])
    rest = len(ids) - MAX_LISTED_IDS
    return f"{head} … et {rest} autres" if rest > 0 else head


def remove_self_references(queens: list) -> tuple[list, list[int]]:
    """Retire les liens d'une reine vers elle-même : mother → None, id retiré de drone_queens (liste
    vide = aucune RPM connue), drone_queens_mother → None. Retourne (reines, ids concernés triés)."""
    out, ids = [], []
    for q in queens:
        update = {}
        if q.mother == q.id:
            update["mother"] = None
        m = q.mating
        if m is not None and (q.id in m.drone_queens or m.drone_queens_mother == q.id):
            update["mating"] = m.model_copy(update={
                "drone_queens": [p for p in m.drone_queens if p != q.id],
                "drone_queens_mother": None if m.drone_queens_mother == q.id else m.drone_queens_mother,
            })
        if update:
            q = q.model_copy(update=update)
            ids.append(q.id)
        out.append(q)
    return out, sorted(ids)


@dataclass(frozen=True)
class _BaseQueen:
    """Reine absente de l'entrée, ajoutée comme reine de base (sans parents ni données)."""

    id: int
    mother: None = None
    mating: None = None


def _sire_spec(mating: Optional["MatingIn"], queens: dict[int, Any]) -> Optional[tuple]:
    """('Q', id_rpm) | ('S', clé, mère_des_rpm, S) | None (sire inconnu)."""
    if mating is None or mating.type in ("natural", "unknown"):
        return None
    dq = list(dict.fromkeys(mating.drone_queens))
    if len(dq) == 1:
        return ("Q", dq[0])
    dam = mating.drone_queens_mother
    if dam is None and len(dq) >= 2:
        mothers = {getattr(queens.get(p), "mother", None) for p in dq}
        if len(mothers) == 1 and None not in mothers:
            dam = mothers.pop()
    if dam is None and not dq and mating.external_ref is None and mating.station is None:
        return None
    n_rpm = len(dq) if dq else (mating.n_drone_queens or 1)
    if mating.external_ref is not None:
        key = ("ext", mating.external_ref)
    else:
        key = ("grp", mating.station, mating.year, dam, tuple(sorted(dq)), n_rpm)
    return ("S", key, dam, n_rpm)


def _sire_node_key(spec: Optional[tuple]) -> Optional[tuple]:
    if spec is None:
        return None
    return ("Q", spec[1]) if spec[0] == "Q" else ("S", spec[1])


def _referenced_ids(q) -> list[int]:
    ids = []
    if q.mother is not None:
        ids.append(q.mother)
    if q.mating is not None and q.mating.type not in ("natural", "unknown"):
        ids.extend(q.mating.drone_queens)
        if q.mating.drone_queens_mother is not None:
            ids.append(q.mating.drone_queens_mother)
    return ids


def _graph(queens: list) -> tuple[dict[int, Any], list[int], dict, dict, dict]:
    """Construit le graphe des nœuds. Retourne (reines, ids ajoutés, specs, dépendances, infos S)."""
    by_id: dict[int, Any] = {q.id: q for q in queens}
    added: list[int] = []
    for q in queens:
        for pid in _referenced_ids(q):
            if pid not in by_id:
                by_id[pid] = _BaseQueen(pid)
                added.append(pid)
    specs = {qid: _sire_spec(q.mating, by_id) for qid, q in by_id.items()}
    deps: dict[tuple, list[tuple]] = {}
    s_info: dict[tuple, tuple] = {}

    def block_deps(mother_id: int) -> list[tuple]:
        out = [("Q", mother_id)]
        spec = specs.get(mother_id)
        if spec is not None:
            out += [("W", mother_id), _sire_node_key(spec)]
        return out

    for qid, q in by_id.items():
        spec = specs[qid]
        if spec is not None and spec[0] == "S":
            skey = ("S", spec[1])
            if skey not in deps:
                s_info[skey] = (spec[2], spec[3])
                deps[skey] = block_deps(spec[2]) if spec[2] is not None else []
        deps[("Q", qid)] = block_deps(q.mother) if q.mother is not None else []
        wdeps = [("Q", qid)]
        if spec is not None:
            wdeps.append(_sire_node_key(spec))
        deps[("W", qid)] = wdeps
    return by_id, added, specs, deps, s_info


def _toposort(deps: dict[tuple, list[tuple]]) -> tuple[list[tuple], list[tuple]]:
    """Tri de Kahn ; retourne (ordre, nœuds pris dans un cycle)."""
    children: dict[tuple, list[tuple]] = {k: [] for k in deps}
    indeg = {k: 0 for k in deps}
    for k, ds in deps.items():
        for d in ds:
            children[d].append(k)
            indeg[k] += 1
    ready = deque(k for k in deps if indeg[k] == 0)
    order: list[tuple] = []
    while ready:
        k = ready.popleft()
        order.append(k)
        for c in children[k]:
            indeg[c] -= 1
            if indeg[c] == 0:
                ready.append(c)
    if len(order) == len(deps):
        return order, []
    rest = {k for k in deps if indeg[k] > 0}
    changed = True
    while changed:  # retire les nœuds qui ne font que dépendre d'un cycle
        changed = False
        for k in list(rest):
            if not any(c in rest for c in children[k]):
                rest.discard(k)
                changed = True
    return order, sorted(rest, key=repr)


def _queen_of_key(key: tuple, s_info: dict) -> Optional[int]:
    if key[0] in ("Q", "W"):
        return key[1]
    return s_info.get(key, (None,))[0]


def find_cycle(queens: list) -> list[int]:
    """Identifiants des reines prises dans un cycle du pedigree (liste vide si aucun).

    Les auto-références sont retirées d'abord : seuls les cycles de longueur ≥ 2 sont signalés."""
    _, _, _, deps, s_info = _graph(remove_self_references(queens)[0])
    _, cyc = _toposort(deps)
    return sorted({q for q in (_queen_of_key(k, s_info) for k in cyc) if q is not None})


@dataclass
class BeePedigree:
    kind: list[str] = field(default_factory=list)
    dam: list[Optional[int]] = field(default_factory=list)
    sire: list[Optional[int]] = field(default_factory=list)
    S: list[int] = field(default_factory=list)
    D: list[int] = field(default_factory=list)
    label: list[tuple] = field(default_factory=list)
    queen_ids: list[int] = field(default_factory=list)
    added_queen_ids: list[int] = field(default_factory=list)
    _q: dict[int, int] = field(default_factory=dict)
    _w: dict[int, int] = field(default_factory=dict)
    _sire: dict[int, Optional[int]] = field(default_factory=dict)
    cache: dict = field(default_factory=dict, repr=False)

    @property
    def n_nodes(self) -> int:
        return len(self.kind)

    def q_node(self, queen_id: int) -> int:
        return self._q[queen_id]

    def w_node(self, queen_id: int) -> Optional[int]:
        return self._w.get(queen_id)

    def sire_of_mating(self, queen_id: int) -> Optional[int]:
        return self._sire.get(queen_id)

    def queen_of_node(self, node: int) -> Optional[int]:
        lab = self.label[node]
        return lab[1] if lab[0] in ("Q", "W") else None


def build_pedigree(queens: list["QueenIn"], default_n_drones: int) -> tuple[BeePedigree, list["WarningOut"]]:
    from .schemas import WarningOut

    queens, self_refs = remove_self_references(queens)
    by_id, added, specs, deps, s_info = _graph(queens)
    order, cyc = _toposort(deps)
    if cyc:
        qs = sorted({q for q in (_queen_of_key(k, s_info) for k in cyc) if q is not None})
        raise ValueError(f"pedigree cyclique : reines concernées {listed_ids(qs)}")

    def n_drones(qid: Optional[int]) -> int:
        if qid is None:
            return default_n_drones
        m = by_id[qid].mating
        return m.n_drones if (m is not None and m.n_drones) else default_n_drones

    ped = BeePedigree()
    index: dict[tuple, int] = {}
    for key in order:
        i = len(ped.kind)
        index[key] = i
        if key[0] == "Q":
            q = by_id[key[1]]
            dam = index[("Q", q.mother)] if q.mother is not None else None
            sk = _sire_node_key(specs.get(q.mother)) if q.mother is not None else None
            ped.kind.append("Q")
            ped.dam.append(dam)
            ped.sire.append(index[sk] if sk is not None else None)
            ped.S.append(1)
            ped.D.append(n_drones(q.mother))
            ped._q[key[1]] = i
        elif key[0] == "W":
            sk = _sire_node_key(specs[key[1]])
            ped.kind.append("W")
            ped.dam.append(index[("Q", key[1])])
            ped.sire.append(index[sk] if sk is not None else None)
            ped.S.append(1)
            ped.D.append(n_drones(key[1]))
            ped._w[key[1]] = i
        else:
            mother, n_rpm = s_info[key]
            ped.kind.append("S")
            ped.dam.append(index[("Q", mother)] if mother is not None else None)
            sk = _sire_node_key(specs.get(mother)) if mother is not None else None
            ped.sire.append(index[sk] if sk is not None else None)
            ped.S.append(int(n_rpm))
            ped.D.append(n_drones(mother))
        ped.label.append(key)
    for qid in by_id:
        sk = _sire_node_key(specs[qid])
        ped._sire[qid] = index[sk] if sk is not None else None
    ped.queen_ids = [q.id for q in queens]
    ped.added_queen_ids = list(added)
    warnings: list[WarningOut] = []
    if self_refs:
        warnings.append(
            WarningOut(
                code="self_reference_removed",
                message=f"{len(self_refs)} reine(s) citée(s) comme sa propre mère, sa propre RPM ou la mère "
                f"de ses RPM : lien retiré du pedigree ({listed_ids(self_refs)})",
            )
        )
    if added:
        warnings.append(
            WarningOut(
                code="unknown_parent_added",
                message=f"{len(added)} reine(s) citée(s) comme parent mais absente(s) de queens, "
                f"ajoutée(s) comme reines de base : ({listed_ids(sorted(added))})",
            )
        )
    return ped, warnings
