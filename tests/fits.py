"""CriterionFit construits à la main pour tester index, accouplements et tendance sans modèle."""
from __future__ import annotations

from blup2.model_standard import CriterionFit
from blup2.schemas import CriterionParameters


def make_fit(cid, ebv, ebv_queen=None, status=None, sigma=1.0, sigma_queen=None, reliability=None):
    params = CriterionParameters(source="input", h2=0.25, repeatability=0.4, sigma2_a=sigma ** 2, sigma2_pe=0.1,
                                 sigma2_e=1.0, n_records=0, n_queens=0, n_contemporary_groups=1,
                                 reliability_method="exact")
    return CriterionFit(
        criterion_id=cid, parameters=params, ebv=dict(ebv), ebv_queen=dict(ebv_queen or ebv),
        reliability=dict(reliability or {q: 0.5 for q in ebv}),
        status=dict(status or {q: "evaluated" for q in ebv}),
        sigma_index=sigma, sigma_queen=sigma_queen if sigma_queen is not None else sigma,
    )
