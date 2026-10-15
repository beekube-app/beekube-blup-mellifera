"""Moteur v2 : pedigree -> parenté -> modèle par critère -> index -> accouplements."""
from __future__ import annotations

import logging
from typing import Optional

from . import ENGINE_VERSION
from .connectivity import connectivity
from .external import expand_external
from .indexes import compute_indexes, queen_outputs
from .kinship import inbreeding, inverse_relationship
from .matings import evaluate_pairs, suggest_matings
from .model_standard import CriterionFit, FitParams, fit_standard
from .pedigree import build_pedigree
from .records import prepare_criterion
from .schemas import BlupV2Input, BlupV2Output, MatingV2Input, MatingV2Output, WarningOut
from .trend import genetic_trend

logger = logging.getLogger(__name__)

MIN_DISCONNECTED_SHARE = 0.10   # disconnected_groups si au moins 10 % des reines notées sont hors du bloc principal


def _fitter(profile: str, warnings: list[WarningOut]):
    if profile == "breeder":
        try:
            from .model_breeder import fit_breeder  # plan « profil breeder »
            return fit_breeder
        except ModuleNotFoundError as exc:
            if exc.name != "blup2.model_breeder":
                raise
            warnings.append(WarningOut(code="breeder_fallback", message="profil breeder non disponible dans "
                                       "cette version du moteur : profil standard utilisé"))
    return fit_standard


def _evaluate(inp: BlupV2Input):
    ext = expand_external(inp)      # ancêtres Index Mellifera : dans le pedigree seulement (spec IM §8.3)
    extra_queens = ext.ancestors + ext.contemporaries             # et leurs notes (spec volet 3 §5)
    extra_records = ext.records
    try:
        ped, warnings = build_pedigree(ext.queens, inp.options.default_n_drones, extra_queens=extra_queens)
        warnings.extend(ext.warnings)
    except ValueError:
        if not extra_queens:
            raise
        logger.exception("pedigree invalide avec les ancêtres Index Mellifera : calcul sans eux")
        ped, warnings = build_pedigree(inp.queens, inp.options.default_n_drones)
        extra_queens, extra_records = [], []
        warnings.extend(w for w in ext.warnings if w.code not in ("external_ref_resolved", "external_records_added"))
        warnings.append(WarningOut(code="external_pedigree_unavailable", message=(
            "ascendance Index Mellifera incompatible avec le pedigree de l'exploitation : calcul sans les ancêtres")))
    Ainv = inverse_relationship(ped)
    F = inbreeding(ped)
    f_workers = {q: float(F[ped.w_node(q)]) for q in (*ped.queen_ids, *ped.added_queen_ids)}
    fit_fn = _fitter(inp.profile, warnings)
    external_nodes = 0
    if ext.contemporaries and extra_queens:    # nœuds des seules contemporaines (spec volet 3 §5)
        external_nodes = ped.n_nodes - build_pedigree(ext.queens, inp.options.default_n_drones,
                                                      extra_queens=ext.ancestors)[0].n_nodes
    params = FitParams(estimate_variances=inp.options.estimate_variances, external_nodes=external_nodes)
    fits: dict[int, Optional[CriterionFit]] = {}
    for c in inp.criteria:
        try:
            cdata = prepare_criterion(inp, c.id, f_workers=f_workers, extra_queens=extra_queens,
                                      extra_records=extra_records)
            fit = fit_fn(cdata, ped, Ainv, params)
            fits[c.id] = fit
            warnings.extend(fit.warnings)
        except Exception as exc:  # un critère en échec n'empêche pas les autres
            logger.exception("échec du calcul du critère %s", c.id)
            fits[c.id] = None
            warnings.append(WarningOut(code="criterion_failed", criterion=c.id,
                                       message=f"échec du calcul : {type(exc).__name__}: {exc}"))
    table = compute_indexes(inp, fits)
    return ped, F, fits, table, warnings


def run_blup(inp: BlupV2Input) -> BlupV2Output:
    ped, F, fits, table, warnings = _evaluate(inp)
    conn = connectivity(inp, ped)
    if conn.n_components > 1 and conn.minor_share >= MIN_DISCONNECTED_SHARE:
        warnings.append(WarningOut(
            code="disconnected_groups",
            message=f"{conn.n_components} blocs de ruchers non reliés (ni reine commune ni apparentées "
                    f"a ≥ 0,125), {round(100 * conn.minor_share)} % des reines notées hors du bloc principal : "
                    f"leurs index ne sont pas comparables"))
    return BlupV2Output(
        engine_version=ENGINE_VERSION,
        profile=inp.profile,
        parameters={str(cid): f.parameters for cid, f in fits.items() if f is not None},
        base=table.base,
        queens=queen_outputs(inp, ped, fits, table, F),
        matings=suggest_matings(inp, ped, table),
        trend=genetic_trend(inp, table),
        connectivity=conn,
        warnings=warnings,
    )


def run_matings(inp: MatingV2Input) -> MatingV2Output:
    ped, F, fits, table, warnings = _evaluate(inp)
    return MatingV2Output(engine_version=ENGINE_VERSION, profile=inp.profile,
                          matings=evaluate_pairs(inp, ped, table, inp.pairs), warnings=warnings)
