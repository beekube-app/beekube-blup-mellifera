"""Moteur v2 : pedigree -> parenté -> modèle par critère -> index -> accouplements."""
from __future__ import annotations

from typing import Optional

from . import ENGINE_VERSION
from .connectivity import connectivity
from .indexes import compute_indexes, queen_outputs
from .kinship import inbreeding, inverse_relationship
from .matings import evaluate_pairs, suggest_matings
from .model_standard import CriterionFit, FitParams, fit_standard
from .pedigree import build_pedigree
from .records import prepare_criterion
from .schemas import BlupV2Input, BlupV2Output, MatingV2Input, MatingV2Output, WarningOut
from .trend import genetic_trend

MIN_DISCONNECTED_SHARE = 0.10   # disconnected_groups si au moins 10 % des reines notées sont hors du bloc principal


def _fitter(profile: str, warnings: list[WarningOut]):
    if profile == "breeder":
        try:
            from .model_breeder import fit_breeder  # plan « profil breeder »
            return fit_breeder
        except ImportError:
            warnings.append(WarningOut(code="breeder_fallback", message="profil breeder non disponible dans "
                                       "cette version du moteur : profil standard utilisé"))
    return fit_standard


def _evaluate(inp: BlupV2Input):
    ped, warnings = build_pedigree(inp.queens, inp.options.default_n_drones)
    Ainv = inverse_relationship(ped)
    F = inbreeding(ped)
    f_workers = {q: float(F[ped.w_node(q)]) for q in ped.queen_ids}
    fit_fn = _fitter(inp.profile, warnings)
    params = FitParams(estimate_variances=inp.options.estimate_variances)
    fits: dict[int, Optional[CriterionFit]] = {}
    for c in inp.criteria:
        try:
            cdata = prepare_criterion(inp, c.id, f_workers=f_workers)
            fit = fit_fn(cdata, ped, Ainv, params)
            fits[c.id] = fit
            warnings.extend(fit.warnings)
        except Exception as exc:  # un critère en échec n'empêche pas les autres
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
