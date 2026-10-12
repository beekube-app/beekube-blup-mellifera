# -*- coding: utf-8 -*-
"""
Vérification de bee_kinship.py contre les formules fermées de Brascamp & Bijma 2014
(l'article ne contient PAS d'exemple numérique de matrice A : on vérifie les équations).
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bee_kinship import BeePedigree, getA_diploid

S, D = 8, 12
np.set_printoptions(precision=5, suppress=True, linewidth=150)


def beebreed_pedigree(S=8, D=12, poisson=True):
    """Pedigree de type Beebreed (Figure 1, p.2) sur 3 générations avec consanguinité."""
    P = BeePedigree(S_default=S, D_default=D, poisson=poisson)
    # génération 0 : sires de base et reines de base
    for s in ["S0a", "S0b"]:
        P.add(s, "S")
    for q in ["Q0_1", "Q0_2", "Q0_3", "Q0_4"]:
        P.add(q, "Q")
    # groupes d'ouvrières des reines de base (accouplées aux sires de base)
    P.add("W_Q0_1", "W", "Q0_1", "S0a")
    P.add("W_Q0_2", "W", "Q0_2", "S0a")
    P.add("W_Q0_3", "W", "Q0_3", "S0b")
    P.add("W_Q0_4", "W", "Q0_4", "S0b")
    # génération 1 : un sire issu de Q0_1 x S0a ; reines pleines soeurs Q0_2 x S0a et Q0_3 x S0b
    P.add("S1", "S", "Q0_1", "S0a")
    for k in range(1, 4):
        P.add(f"Q1_{k}", "Q", "Q0_2", "S0a")
    for k in range(4, 7):
        P.add(f"Q1_{k}", "Q", "Q0_3", "S0b")
    P.add("Q1_7", "Q", "Q0_4", "S0b")
    # les reines de génération 1 sont accouplées au sire S1 (apparenté à Q1_1..3 : même sire S0a)
    for k in range(1, 8):
        P.add(f"W_Q1_{k}", "W", f"Q1_{k}", "S1")
    # génération 2 : sire S2 issu de Q1_4 x S1 ; reines issues de Q1_1 x S1 (consanguines) et Q1_5 x S1
    P.add("S2", "S", "Q1_4", "S1")
    for k in range(1, 4):
        P.add(f"Q2_{k}", "Q", "Q1_1", "S1")
    for k in range(4, 7):
        P.add(f"Q2_{k}", "Q", "Q1_5", "S1")
    for k in range(1, 7):
        P.add(f"W_Q2_{k}", "W", f"Q2_{k}", "S2")
    # génération 3 : reines issues de Q2_1 x S2 (fortement consanguines)
    for k in range(1, 3):
        P.add(f"Q3_{k}", "Q", "Q2_1", "S2")
    return P


def check(label, got, exp, tol=1e-10):
    ok = abs(got - exp) < tol
    print(f"  [{'OK' if ok else 'ECART'}] {label:70s} calc={got:.6f} attendu={exp:.6f}")
    return ok


def main():
    P = beebreed_pedigree(S, D)
    A = P.build()
    Ainv = P.ainv()
    n = len(P.order)
    print(f"Pedigree Beebreed-like : {n} individus "
          f"({sum(i.kind=='Q' for i in P.order)} reines, {sum(i.kind=='S' for i in P.order)} sires, "
          f"{sum(i.kind=='W' for i in P.order)} groupes d'ouvrières)")
    allok = True

    print("\n1. Identités matricielles")
    A2 = P.T @ P.D @ P.T.T
    e1 = np.abs(A - A2).max()
    print(f"  max|A - T D T'|                 = {e1:.2e}")
    e2 = np.abs(Ainv - np.linalg.inv(A)).max()
    print(f"  max|(I-M)'D^-1(I-M) - inv(A)|   = {e2:.2e}")
    e3 = np.abs(A @ Ainv - np.eye(n)).max()
    print(f"  max|A A^-1 - I|                 = {e3:.2e}")
    ev = np.linalg.eigvalsh(A).min()
    print(f"  plus petite valeur propre de A  = {ev:.4f}  (définie positive : {ev > 0})")
    allok &= e1 < 1e-10 and e2 < 1e-8 and e3 < 1e-8

    a = lambda x, y: A[P.inds[x].idx, P.inds[y].idx]
    print("\n2. Formules fermées (parents de base, F=0, ass=0)")
    allok &= check("sire de base : a_ss = 1/S  (Annexe 1 cas 2, éq. 12)", a("S0a", "S0a"), 1 / S)
    allok &= check("reine de base : a_ii = 1", a("Q0_1", "Q0_1"), 1.0)
    allok &= check("mère-fille : 1/2", a("Q0_2", "Q1_1"), 0.5)
    allok &= check("reine - son sire (groupe RPM) : 1/2 a_ss = 1/(2S)", a("Q1_1", "S0a"), 0.5 / S)
    aFS = 0.25 + 1 / (2 * D) + 1 / (4 * S)
    allok &= check("pleines soeurs : 1/4 + 1/(2D) + 1/(4S)  (éq. 24a, p.7)", a("Q1_1", "Q1_2"), aFS)
    allok &= check("groupe d'ouvrières de Q0_2 : a_ww = a_FS (var(Abar_w))", a("W_Q0_2", "W_Q0_2"), aFS)
    allok &= check("reine - groupe d'ouvrières de sa mère : = a_FS", a("Q1_1", "W_Q0_2"), aFS)
    allok &= check("demi-soeurs paternelles (Q1_1, Q1_7 ? non: S0a vs S0b) -> 0", a("Q1_1", "Q1_7"), 0.0)
    allok &= check("demi-soeurs paternelles Q1_1 (Q0_2xS0a) / S1-RPM... via Q0_1xS0a : 1/(4S)",
                   a("Q1_1", "S1"), 0.5 * a("Q1_1", "Q0_1") + 0.5 * a("Q1_1", "S0a"))
    allok &= check("demi-soeurs paternelles (même sire, mères différentes) : 1/(4S) [p.8 apHS, ass=0]",
                   a("Q1_4", "Q1_7"), (1 + (S - 1) * 0) / (4 * S))
    allok &= check("D_queen éq.13 base : 1/4 + 1/4 + (S-1)/(4S)", P.D[P.inds["Q1_1"].idx, P.inds["Q1_1"].idx],
                   0.25 + 0.25 + 0.25 * (S - 1) / S)
    covFS = (2 - 0) / (4 * D)
    allok &= check("cov(dFS) Poisson éq.21a = (2-ass)/(4D)", P.block_cov[("Q0_2", "S0a")], covFS)
    allok &= check("D_worker = cov(dFS) (éq.19-20)", P.D[P.inds["W_Q0_2"].idx, P.inds["W_Q0_2"].idx], covFS)
    allok &= check("D_sire éq.17 = var(d)/S + (S-1)/S cov(dFS)", P.D[P.inds["S1"].idx, P.inds["S1"].idx],
                   (0.5 + 0.25 * (S - 1) / S) / S + (S - 1) / S * covFS)
    allok &= check("hors-diag. D pleines soeurs = cov(dFS)", P.D[P.inds["Q1_1"].idx, P.inds["Q1_2"].idx], covFS)
    allok &= check("hors-diag. D reine / groupe d'ouvrières de même mère = cov(dFS)",
                   P.D[P.inds["Q1_1"].idx, P.inds["W_Q0_2"].idx], covFS)
    allok &= check("hors-diag. D entre blocs = 0", P.D[P.inds["Q1_1"].idx, P.inds["Q1_4"].idx], 0.0)

    print("\n3. Récursion de a_ss (éq. 24a) et consanguinité (éq. 25) sur les générations suivantes")
    # S1 : RPM filles de Q0_1 x S0a : Fd=0, Fs=0, ass_{i-1}=0, asd=0
    Fs1, ass1, _ = P.sire_attr["S1"]
    ass1_exp = (1 + 0) / 4 + 1 / (2 * D) + (1 + 0) / (4 * S) + (D * S - D - S) * 0 / (4 * D * S) + 0 / 2
    allok &= check("ass(S1) éq.24a", ass1, ass1_exp)
    allok &= check("a_ss(S1) = (1+Fs+(S-1)ass)/S  (éq.12)", a("S1", "S1"), (1 + Fs1 + (S - 1) * ass1) / S)
    # S2 : RPM filles de Q1_4 x S1. Fd = F(Q1_4) = 0 ; Fs = Fs1 ; ass_{i-1} = ass1 ; asd = a(Q1_4, S1)
    Fs2, ass2, _ = P.sire_attr["S2"]
    Fd = a("Q1_4", "Q1_4") - 1
    asd = a("Q1_4", "S1")
    ass2_exp = (1 + Fd) / 4 + 1 / (2 * D) + (1 + Fs1) / (4 * S) + (D * S - D - S) * ass1 / (4 * D * S) + asd / 2
    allok &= check("ass(S2) éq.24a (avec asd != 0)", ass2, ass2_exp)
    allok &= check("Fs(S2) = asd/2 (éq.25)", Fs2, asd / 2)
    allok &= check("a_ss(S2) = (1+Fs+(S-1)ass)/S", a("S2", "S2"), (1 + Fs2 + (S - 1) * ass2) / S)
    # Q2_1..3 : filles de Q1_1 x S1 ; Q1_1 et S1 apparentés (même sire S0a) -> consanguines
    F_Q2 = a("Q2_1", "Q2_1") - 1
    allok &= check("F(Q2_1) = 1/2 a(Q1_1,S1) (éq.25)", F_Q2, 0.5 * a("Q1_1", "S1"))
    # pleines soeurs Q2_1,Q2_2 = formule 24a appliquée à (Q1_1, S1)
    exp = (1 + 0) / 4 + 1 / (2 * D) + (1 + Fs1) / (4 * S) + (D * S - D - S) * ass1 / (4 * D * S) + a("Q1_1", "S1") / 2
    allok &= check("a(Q2_1,Q2_2) = éq.24a appliquée au couple (Q1_1,S1)", a("Q2_1", "Q2_2"), exp)
    # éq. 13 avec consanguinité pour Q3_1
    Fd3 = a("Q2_1", "Q2_1") - 1
    exp13 = 0.25 * (1 - Fd3) + 0.25 * (1 - Fs2) + 0.25 * (S - 1) / S * (1 + Fs2 - ass2)
    allok &= check("D(Q3_1) éq.13 avec Fd, Fs, ass non nuls", P.D[P.inds["Q3_1"].idx, P.inds["Q3_1"].idx], exp13)
    allok &= check("F(Q3_1) = 1/2 a(Q2_1,S2)", a("Q3_1", "Q3_1") - 1, 0.5 * a("Q2_1", "S2"))
    print(f"  valeurs : ass(S1)={ass1:.5f}  ass(S2)={ass2:.5f}  Fs(S2)={Fs2:.5f}  F(Q2_1)={F_Q2:.5f}  F(Q3_1)={a('Q3_1','Q3_1')-1:.5f}")

    print("\n4. Variante BER des probabilités (éq. 23b -> 21b, 24b)")
    Pb = beebreed_pedigree(S, D, poisson=False)
    Ab = Pb.build()
    ab = lambda x, y: Ab[Pb.inds[x].idx, Pb.inds[y].idx]
    exp24b = 0.25 + 1 / (2 * D) + (D - 1) * (1 + 0) / (4 * D * S) + 0
    allok &= check("pleines soeurs (base) éq.24b", ab("Q1_1", "Q1_2"), exp24b)
    allok &= check("cov(dFS) éq.21b", Pb.block_cov[("Q0_2", "S0a")], (1 - 0) / (4 * D) + (1 + 0 - 0) * (S - 1) / (4 * D * S))
    # valeur asymptotique BER (p.8) : a_FS equilibre de 24b avec Fs=Fd=asd=0
    aeq = (0.25 + 1 / (2 * D) + (D - 1) / (4 * D * S)) / (1 - (D - 1) * (S - 1) / (4 * D * S))
    apHS = (1 + (S - 1) * aeq) / (4 * S)
    q = np.sqrt(((aeq - 0.25) + apHS) / 2)
    print(f"  asymptote BER a_FS = {aeq:.4f}, a_pFS = {aeq-0.25:.4f}, a_pHS = {apHS:.4f}, q = {q:.4f}  (éq. 27-28, p.8)")

    print("\n5. Sire inconnu : cas 3 de l'Annexe 1 == matrice diploïde 'mère seule' (getA)")
    P3 = BeePedigree(S, D)
    for q in ["A", "B"]:
        P3.add(q, "Q")
    P3.add("C", "Q", "A", None)
    P3.add("E", "Q", "A", None)
    P3.add("G", "Q", "C", None)
    A3 = P3.build()
    Ad = getA_diploid(["A", "B", "C", "E", "G"], {"C": "A", "E": "A", "G": "C"}, {})
    e = np.abs(A3 - Ad).max()
    print(f"  max|A_BB(sire inconnu) - getA_diploïde| = {e:.2e}")
    allok &= e < 1e-12

    print("\n6. Gene-dropping Monte-Carlo (1 mère x 1 sire de S RPM, D mâles ~ Poisson)")
    rng = np.random.default_rng(1)
    nrep = 20000
    fs_vals, hs_vals = [], []
    for r in range(nrep):
        # nombre de mâles par RPM ~ Poisson(D/S), contribution de chaque mâle ~ Poisson(grand) -> ~ égale
        ndr = rng.poisson(D / S, size=S)
        if ndr.sum() == 0:
            continue
        drones_rpm = np.repeat(np.arange(S), ndr)           # RPM d'origine de chaque mâle
        w = rng.poisson(200, size=len(drones_rpm)).astype(float)
        w /= w.sum()
        # deux filles de la même mère : allèles maternels (0/1), mâle tiré selon w
        pat_x, pat_y = rng.choice(len(drones_rpm), size=2, p=w)
        same_drone = pat_x == pat_y
        same_rpm = drones_rpm[pat_x] == drones_rpm[pat_y]
        # allèle du mâle = un des deux allèles de sa RPM (haploïde : 1 allèle)
        al_x = rng.integers(2); al_y = al_x if same_drone else rng.integers(2)
        pat_ibd = same_drone or (same_rpm and al_x == al_y)
        mat_ibd = rng.integers(2) == rng.integers(2)
        # coancestrie = 1/4 [P(mat=mat) + P(pat=pat)] (allèles mat/pat distincts entre fondateurs)
        fs_vals.append(0.5 * (int(mat_ibd) + int(pat_ibd)))   # a = 2 * coancestrie
        # demi-soeurs paternelles : mères différentes, même sire, mâles différents (sauf hasard : mères
        # différentes => mâles différents forcément)
        ndr2 = rng.poisson(D / S, size=S)
        if ndr2.sum() == 0:
            continue
        rpm_x = rng.choice(np.repeat(np.arange(S), ndr))
        rpm_y = rng.choice(np.repeat(np.arange(S), ndr2))
        hs_vals.append(0.5 * int((rpm_x == rpm_y) and (rng.integers(2) == rng.integers(2))))
    fs, hs = np.mean(fs_vals), np.mean(hs_vals)
    se_fs, se_hs = np.std(fs_vals) / np.sqrt(len(fs_vals)), np.std(hs_vals) / np.sqrt(len(hs_vals))
    print(f"  pleines soeurs : simulé {fs:.4f} ± {se_fs:.4f}   théorie éq.24a {aFS:.4f}   (BER 24b {exp24b:.4f})")
    print(f"  demi-soeurs pat.: simulé {hs:.4f} ± {se_hs:.4f}   théorie 1/(4S) = {1/(4*S):.4f}")

    print("\nRESULTAT GLOBAL :", "toutes les vérifications passent" if allok else "DES ECARTS EXISTENT")
    return allok


if __name__ == "__main__":
    main()
