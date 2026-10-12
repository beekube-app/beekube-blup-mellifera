# -*- coding: utf-8 -*-
"""
Simulation : précision des index selon la structure de parenté utilisée.
Vraies valeurs génétiques ~ N(0, A_BB sigma2_A) (reines, sires, groupes d'ouvrières).
Phénotype d'une colonie (2-5 notes) : y = mu + Abar_w(colonie) + pe(colonie) + rucher + e.
Méthodes (variances connues, équations du modèle mixte en numpy) :
  a     BLUP abeille, A de Brascamp & Bijma, sires et accouplements connus
  a_ph  BLUP abeille, sires/accouplements INCONNUS -> sire fantôme de base par reine (faisable avec les données Beekube)
  b     modèle de répétabilité : effet reine iid (= blup.r actuel, has_drone_info FALSE)
  b2    effet reine iid + effet 'drone_parent' iid (= blup.r avec drone_parent renseigné)
  c     modèle animal diploïde classique : drone_parent traité comme un père diploïde
  d0    modèle animal diploïde 'mère seule' (sires inconnus)
"""
import sys, os, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bee_kinship import BeePedigree, getA_diploid


# ----------------------------------------------------------------------------
def gen_pedigree(rng, G=4, n0=30, ns=6, n_dams=15, kmin=3, kmax=5):
    """Retourne la liste ordonnée d'individus (dict) : kind Q/S/W, dam, sire, gen, mate."""
    inds = []
    queens_by_gen = {}
    sires_by_gen = {}
    mate = {}

    def add(id, kind, dam=None, sire=None, gen=0):
        inds.append(dict(id=id, kind=kind, dam=dam, sire=sire, gen=gen))
        return id

    # génération 0
    sires_by_gen[0] = [add(f"S0_{k}", "S", gen=0) for k in range(ns)]
    queens_by_gen[0] = [add(f"Q0_{k}", "Q", gen=0) for k in range(n0)]
    for q in queens_by_gen[0]:
        mate[q] = rng.choice(sires_by_gen[0])
        add("W_" + q, "W", q, mate[q], 0)
    for g in range(1, G + 1):
        prev = queens_by_gen[g - 1]
        # sires de la génération g : filles d'une reine de g-1 (et de son accouplement)
        sdams = rng.choice(prev, size=ns, replace=False)
        sires_by_gen[g] = [add(f"S{g}_{k}", "S", d, mate[d], g) for k, d in enumerate(sdams)]
        dams = rng.choice(prev, size=n_dams, replace=False)
        queens_by_gen[g] = []
        for d in dams:
            for j in range(rng.integers(kmin, kmax + 1)):
                q = add(f"Q{g}_{d}_{j}", "Q", d, mate[d], g)
                queens_by_gen[g].append(q)
        for q in queens_by_gen[g]:
            mate[q] = rng.choice(sires_by_gen[g])
            add("W_" + q, "W", q, mate[q], g)
    return inds, mate


def build_bb(inds, S, D, drop_sires=False, phantom=False):
    """A_BB. drop_sires : les sires sont inconnus. phantom : sire fantôme de base par reine-mère/accouplée."""
    P = BeePedigree(S_default=S, D_default=D)
    if not drop_sires:
        for x in inds:
            P.add(x["id"], x["kind"], x["dam"], x["sire"])
        return P
    # sires inconnus : on ne garde que reines (et groupes d'ouvrières si fantôme)
    for x in inds:
        if x["kind"] == "Q":
            if phantom:
                if x["dam"] is not None and ("PH_" + x["dam"]) not in P.inds:
                    P.add("PH_" + x["dam"], "S")
                P.add(x["id"], "Q", x["dam"], ("PH_" + x["dam"]) if (phantom and x["dam"] is not None) else None)
        elif x["kind"] == "W" and phantom:
            q = x["dam"]
            if ("PH_" + q) not in P.inds:
                P.add("PH_" + q, "S")
            P.add(x["id"], "W", q, "PH_" + q)
    return P


def mme_solve(y, X, blocks, s2e):
    """blocks : liste de (Z, Ginv) avec Ginv = inverse de la matrice de covariance (déjà /sigma2).
    Retourne (beta, [u_k])."""
    Zs = [Z for Z, _ in blocks]
    W = np.hstack([X] + Zs)
    LHS = W.T @ W
    RHS = W.T @ y
    off = X.shape[1]
    for Z, Ginv in blocks:
        k = Z.shape[1]
        LHS[off:off + k, off:off + k] += Ginv * s2e
        off += k
    sol = np.linalg.solve(LHS, RHS)
    out, off = [], X.shape[1]
    beta = sol[:off]
    for Z, _ in blocks:
        k = Z.shape[1]
        out.append(sol[off:off + k]); off += k
    return beta, out


def incidence(levels_of_obs, levels):
    pos = {l: i for i, l in enumerate(levels)}
    Z = np.zeros((len(levels_of_obs), len(levels)))
    for r, l in enumerate(levels_of_obs):
        Z[r, pos[l]] = 1.0
    return Z


def stats(ebv, tbv):
    if np.std(ebv) < 1e-12:
        return np.nan, np.nan
    c = np.corrcoef(ebv, tbv)[0, 1]
    b = np.cov(tbv, ebv)[0, 1] / np.var(ebv, ddof=1)
    return c, b


def one_rep(rng, S=8, D=12, s2a=1.0, s2pe=0.5, s2ap=0.5, s2e=1.5, p_notes=0.7, n_ap=4,
            analysis=None):
    """analysis : dict de variantes supplémentaires {name: (S_an, D_an, drop, phantom)}"""
    inds, mate = gen_pedigree(rng)
    P = build_bb(inds, S, D)
    A = P.build()
    n = len(inds)
    L = np.linalg.cholesky(A + 1e-12 * np.eye(n))
    a_true = L @ rng.standard_normal(n) * np.sqrt(s2a)
    idx = {x["id"]: i for i, x in enumerate(inds)}
    queens = [x["id"] for x in inds if x["kind"] == "Q"]
    # notes
    rec_q = [q for q in queens if rng.random() < p_notes]
    apiary_of = {q: rng.integers(n_ap) for q in queens}
    ap_eff = rng.standard_normal(n_ap) * np.sqrt(s2ap)
    pe = {q: rng.standard_normal() * np.sqrt(s2pe) for q in rec_q}
    obs_q, y = [], []
    for q in rec_q:
        for _ in range(rng.integers(2, 6)):
            obs_q.append(q)
            y.append(5.0 + a_true[idx["W_" + q]] + pe[q] + ap_eff[apiary_of[q]] + rng.standard_normal() * np.sqrt(s2e))
    y = np.array(y)
    X = np.ones((len(y), 1))
    Zap = incidence([apiary_of[q] for q in obs_q], list(range(n_ap)))
    Zpe = incidence(obs_q, rec_q)
    var_w = s2a * np.mean([A[idx["W_" + q], idx["W_" + q]] for q in queens])   # variance génétique d'une colonie
    res = {}
    tw = np.array([a_true[idx["W_" + q]] for q in queens])   # vraie valeur de la colonie / des filles (éq. 2)
    tq = np.array([a_true[idx[q]] for q in queens])          # vraie valeur propre de la reine
    has = np.array([q in rec_q for q in queens])

    def record(name, ebv_col, ebv_q=None):
        if ebv_q is None:
            ebv_q = ebv_col
        r = {}
        r["c_col_notes"], r["b_col_notes"] = stats(ebv_col[has], tw[has])
        r["c_q_notes"], _ = stats(ebv_q[has], tq[has])
        r["c_col_sans"], _ = stats(ebv_col[~has], tw[~has])
        r["c_q_sans"], _ = stats(ebv_q[~has], tq[~has])
        # réponse : moyenne vraie (colonie) des 10 % meilleures reines notées selon l'index
        k = max(1, int(0.1 * has.sum()))
        top = np.argsort(-ebv_col[has])[:k]
        r["resp_top10"] = tw[has][top].mean()
        res[name] = r

    # --- a : BLUP abeille complet ---
    def run_bb(Pm, name, scale=1.0):
        Am_inv = Pm.ainv()
        ids = [i.id for i in Pm.order]
        Zg = incidence(["W_" + q for q in obs_q], ids)
        _, (g, p, r_) = mme_solve(y, X, [(Zg, Am_inv / (s2a * scale)), (Zpe, np.eye(len(rec_q)) / s2pe),
                                        (Zap, np.eye(n_ap) / s2ap)], s2e)
        ebv_col = np.array([g[Pm.inds["W_" + q].idx] for q in queens])
        ebv_q = np.array([g[Pm.inds[q].idx] for q in queens])
        record(name, ebv_col, ebv_q)

    run_bb(P, "a")
    # --- a_ph : sires inconnus, fantôme par reine ---
    run_bb(build_bb(inds, S, D, drop_sires=True, phantom=True), "a_ph")
    # --- b : reine iid (pe confondu) ---
    Zq = incidence(obs_q, rec_q)
    _, (u, r_) = mme_solve(y, X, [(Zq, np.eye(len(rec_q)) / (var_w + s2pe)), (Zap, np.eye(n_ap) / s2ap)], s2e)
    posr = {q: i for i, q in enumerate(rec_q)}
    ebv = np.array([u[posr[q]] if q in posr else 0.0 for q in queens])
    record("b", ebv)
    # --- b2 : reine iid + groupe drone_parent (= sire de la reine) iid ---
    sire_of = {x["id"]: x["sire"] for x in inds}
    groups = sorted({sire_of[q] for q in rec_q if sire_of[q] is not None})
    if groups:
        Zs = incidence([sire_of[q] if sire_of[q] is not None else groups[0] for q in obs_q], groups)
        for r_i, q in enumerate(obs_q):   # reines sans sire : pas d'effet
            if sire_of[q] is None:
                Zs[r_i, :] = 0
        a_ss = s2a * np.mean([A[idx[s], idx[s]] for s in groups])
        v_s = 0.25 * a_ss
        _, (u, sg, r_) = mme_solve(y, X, [(Zq, np.eye(len(rec_q)) / (var_w + s2pe - v_s)),
                                         (Zs, np.eye(len(groups)) / v_s), (Zap, np.eye(n_ap) / s2ap)], s2e)
        pg = {g_: i for i, g_ in enumerate(groups)}
        ebv = np.array([(u[posr[q]] if q in posr else 0.0) + (sg[pg[sire_of[q]]] if sire_of[q] in pg else 0.0)
                        for q in queens])
        record("b2", ebv)
    # --- c : diploïde, drone_parent = père diploïde ---
    ids_c = [x["id"] for x in inds if x["kind"] != "W"]
    dam_c = {x["id"]: x["dam"] for x in inds if x["kind"] != "W"}
    sire_c = {x["id"]: x["sire"] for x in inds if x["kind"] != "W"}
    Ac = getA_diploid(ids_c, dam_c, sire_c)
    Zc = incidence(obs_q, ids_c)
    _, (g, p, r_) = mme_solve(y, X, [(Zc, np.linalg.inv(Ac) / var_w), (Zpe, np.eye(len(rec_q)) / s2pe),
                                    (Zap, np.eye(n_ap) / s2ap)], s2e)
    pc = {k: i for i, k in enumerate(ids_c)}
    record("c", np.array([g[pc[q]] for q in queens]))
    # --- d0 : diploïde mère seule ---
    Ad = getA_diploid(queens, {q: dam_c[q] for q in queens}, {})
    Zd = incidence(obs_q, queens)
    _, (g, p, r_) = mme_solve(y, X, [(Zd, np.linalg.inv(Ad) / var_w), (Zpe, np.eye(len(rec_q)) / s2pe),
                                    (Zap, np.eye(n_ap) / s2ap)], s2e)
    record("d0", g)
    # --- variantes d'analyse (partie 4) ---
    if analysis:
        for name, (S_an, D_an, drop, ph) in analysis.items():
            run_bb(build_bb(inds, S_an, D_an, drop_sires=drop, phantom=ph), name)
    return res


def summarize(all_res, names, keys, title):
    print("\n" + title)
    hdr = f"{'méthode':8s}" + "".join(f"{k:>18s}" for k in keys)
    print(hdr)
    for nm in names:
        line = f"{nm:8s}"
        for k in keys:
            v = np.array([r[nm][k] for r in all_res if nm in r], float)
            v = v[~np.isnan(v)]
            if len(v) == 0:
                line += f"{'n.c.':>18s}"
            else:
                line += f"{v.mean():9.3f} ±{v.std(ddof=1)/np.sqrt(len(v)):6.3f}  "
        print(line)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--s2e", type=float, default=1.5)
    ap.add_argument("--part4", action="store_true")
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    analysis = None
    if args.part4:
        analysis = {
            "aD4": (8, 4, False, False), "aD8": (8, 8, False, False), "aD16": (8, 16, False, False),
            "aD24": (8, 24, False, False), "aD48": (8, 48, False, False),
            "aS1": (1, 12, False, False), "aS4": (4, 12, False, False), "aS16": (16, 12, False, False),
            "aS100": (100, 12, False, False),
            "phS1": (1, 12, True, True), "phS4": (4, 12, True, True), "phD6": (8, 6, True, True),
            "phD24": (8, 24, True, True),
        }
    all_res = []
    for r in range(args.reps):
        all_res.append(one_rep(rng, s2e=args.s2e, analysis=analysis))
    keys = ["c_col_notes", "b_col_notes", "c_q_notes", "c_col_sans", "c_q_sans", "resp_top10"]
    names = ["a", "a_ph", "b", "b2", "c", "d0"]
    print(f"Simulation : {args.reps} répétitions, S=8, D=12, s2A=1, s2pe=0.5, s2rucher=0.5, s2e={args.s2e}, "
          f"70 % des reines notées (2-5 notes)")
    print("Colonnes : c_col_notes = corr(index, vraie valeur colonie/filles) reines notées ; b_col_notes = pente vraie~index ;")
    print("           c_q_notes = corr avec la valeur propre de la reine ; *_sans = reines SANS notes ; resp_top10 = vraie valeur")
    print("           moyenne (colonie) des 10 % meilleures reines notées selon l'index (sigma_A = 1).")
    summarize(all_res, names, keys, "Partie 3 : méthodes")
    if args.part4:
        summarize(all_res, ["a", "aD4", "aD8", "aD16", "aD24", "aD48"], keys,
                  "Partie 4a : A_BB avec D supposé (vrai D = 12)")
        summarize(all_res, ["a", "aS1", "aS4", "aS16", "aS100"], keys,
                  "Partie 4b : A_BB avec S supposé (vrai S = 8)")
        summarize(all_res, ["a", "a_ph", "phS1", "phS4", "phD6", "phD24", "d0", "b"], keys,
                  "Partie 4c : sires inconnus : fantôme (S,D supposés) vs mère seule (d0) vs iid (b)")


if __name__ == "__main__":
    main()
