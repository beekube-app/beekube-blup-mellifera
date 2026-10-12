# -*- coding: utf-8 -*-
"""
Matrice de parenté abeille (Brascamp & Bijma 2014, GSE 46:53) et comparateurs.

Trois types d'individus : 'Q' reine, 'S' « sire » = groupe de S reines productrices
de mâles (RPM) sœurs, 'W' groupe d'ouvrières d'une colonie.
Chaque individu a dam (reine) et sire (groupe) éventuellement inconnus (None).

A = (I - M)^-1 D (I - M)^-T   (éq. 7-8), D bloc-diagonale par couple (dam, sire)
A^-1 = (I - M)' D^-1 (I - M)   (éq. 8)
"""
import numpy as np


class Ind:
    __slots__ = ("id", "kind", "dam", "sire", "S", "D", "idx")

    def __init__(self, id, kind, dam=None, sire=None, S=None, D=None):
        self.id, self.kind, self.dam, self.sire, self.S, self.D = id, kind, dam, sire, S, D
        self.idx = None


class BeePedigree:
    """
    Paramètres :
      S_default : nombre de RPM par groupe de mâles (sire)
      D_default : nombre de mâles accouplés à une reine
      poisson   : True -> p1 = 1/D, p2 = 1/D + 1/S (éq. 22a, 23a)  => cov(dFS) = (2 - ass)/(4D) (éq. 21a)
                  sinon on utilise p2_BER = (S + D - 1)/(DS) (éq. 23b) -> éq. 21b
    """

    def __init__(self, S_default=8, D_default=12, poisson=True):
        self.inds = {}
        self.order = []
        self.S_default = S_default
        self.D_default = D_default
        self.poisson = poisson

    def add(self, id, kind, dam=None, sire=None, S=None, D=None):
        if id in self.inds:
            raise ValueError(f"id dupliqué {id}")
        if kind == "S" and S is None:
            S = self.S_default
        if D is None:
            D = self.D_default
        ind = Ind(id, kind, dam, sire, S, D)
        ind.idx = len(self.order)
        self.inds[id] = ind
        self.order.append(ind)
        return ind

    # ------------------------------------------------------------------
    def _check_order(self):
        seen = set()
        for ind in self.order:
            for p in (ind.dam, ind.sire):
                if p is not None and p not in seen:
                    raise ValueError(f"parent {p} de {ind.id} non défini avant lui (trier le pedigree)")
            seen.add(ind.id)

    def _p12(self, S, D):
        p1 = 1.0 / D
        if self.poisson:
            p2 = min(1.0, 1.0 / D + 1.0 / S)   # S=1 : p2 = 1 exactement (éq. 23, sigma2_cS = 0)
        else:
            p2 = (S + D - 1.0) / (D * S)
        return p1, p2

    def build(self):
        """Construit A, M, D (bloc-diag.), Ainv, F, et les attributs (Fs, ass) des sires."""
        self._check_order()
        n = len(self.order)
        A = np.zeros((n, n))
        T = np.zeros((n, n))      # T = (I - M)^-1, triangulaire inférieure
        M = np.zeros((n, n))
        Dm = np.zeros((n, n))     # matrice D (bloc-diagonale)
        # attributs par sire : (Fs, ass) ; par reine : F (= A_ii - 1)
        sire_attr = {}
        # covariance des termes de ségrégation du bloc (dam, sire) -> covFS
        block_cov = {}
        block_members = {}

        for ind in self.order:
            i = ind.idx
            d = self.inds[ind.dam] if ind.dam is not None else None
            s = self.inds[ind.sire] if ind.sire is not None else None
            Fd = A[d.idx, d.idx] - 1.0 if d is not None else 0.0
            if s is not None:
                Fs, ass, S = sire_attr[s.id]
            else:
                Fs, ass, S = 0.0, 0.0, None
            D = ind.D

            # --- variance de ségrégation d'une reine (éq. 13 et Annexe 1) ---
            if d is not None and s is not None:
                var_q = 0.25 * (1 - Fd) + 0.25 * (1 - Fs) + 0.25 * (S - 1) / S * (1 + Fs - ass)
                p1, p2 = self._p12(S, D)
                covFS = p1 * (1 - Fs) / 4.0 + (p2 - 1.0 / S) * (1 + Fs - ass) / 4.0   # éq. 21
            elif d is not None:   # cas 3 : mère connue, sire inconnu
                var_q = 0.5 + 0.25 * (1 - Fd)
                covFS = 0.0
            elif s is not None:   # cas 4 : sire connu, mère inconnue
                var_q = 0.5 + 0.25 * (1 - Fs) + 0.25 * (S - 1) / S * (1 + Fs - ass)
                covFS = (1 - Fs) / (4.0 * D)   # cas 5 de l'Annexe 1
            else:                 # cas 1 : base
                var_q = 1.0
                covFS = 0.0

            if ind.kind == "Q":
                Dii = var_q
            elif ind.kind == "S":
                Si = ind.S
                if d is None and s is None:
                    Dii = 1.0 / Si                                  # cas 2
                elif d is not None and s is None:
                    Dii = var_q / Si                                # cas 6
                else:
                    Dii = var_q / Si + (Si - 1.0) / Si * covFS      # éq. 17 (et cas 5)
            elif ind.kind == "W":
                # n ouvrières -> infini : var(δ̄) = cov(δFS) (éq. 19)
                if d is None and s is None:
                    raise ValueError("groupe d'ouvrières sans parents")
                Dii = covFS
            else:
                raise ValueError(ind.kind)

            # --- lignes de M et T ---
            if d is not None:
                M[i, d.idx] = 0.5
                T[i, :] += 0.5 * T[d.idx, :]
            if s is not None:
                M[i, s.idx] = 0.5
                T[i, :] += 0.5 * T[s.idx, :]
            T[i, i] = 1.0

            # --- bloc de D : covariance avec les autres descendants du même couple (dam, sire) ---
            key = (ind.dam, ind.sire) if ind.dam is not None else None
            Dm[i, i] = Dii
            if key is not None and covFS != 0.0:
                for k in block_members.get(key, []):
                    Dm[i, k] = Dm[k, i] = covFS
            if key is not None:
                block_members.setdefault(key, []).append(i)
                block_cov[key] = covFS

            # --- ligne i de A : a_ij = 1/2 a_dj + 1/2 a_sj + sum_k D_ik T_jk ---
            row = np.zeros(n)
            if d is not None:
                row += 0.5 * A[d.idx, :]
            if s is not None:
                row += 0.5 * A[s.idx, :]
            row += T @ Dm[i, :]          # sum_k T_jk D_ik  (T_jk = 0 pour j non encore défini)
            A[i, :] = row
            A[:, i] = row
            # a_ii = (T D T')_ii
            A[i, i] = T[i, :] @ Dm @ T[i, :]

            # attributs si c'est un sire : Fs = F des RPM, ass = parenté entre RPM
            if ind.kind == "S":
                if d is not None and s is not None:
                    F_rpm = 0.5 * A[d.idx, s.idx]
                    a_rpm = 0.25 * A[d.idx, d.idx] + 0.25 * A[s.idx, s.idx] + 0.5 * A[d.idx, s.idx] + covFS
                elif d is not None:
                    F_rpm, a_rpm = 0.0, 0.25 * A[d.idx, d.idx]
                elif s is not None:
                    F_rpm, a_rpm = 0.0, 0.25 * A[s.idx, s.idx] + covFS
                else:
                    F_rpm, a_rpm = 0.0, 0.0
                sire_attr[ind.id] = (F_rpm, a_rpm, ind.S)

        self.A, self.M, self.D, self.T = A, M, Dm, T
        self.sire_attr = sire_attr
        self.block_cov = block_cov
        self.block_members = block_members
        return A

    def ainv(self):
        if not hasattr(self, "A"):
            self.build()
        """A^-1 = (I-M)' D^-1 (I-M), D^-1 bloc par bloc (éq. 8)."""
        n = len(self.order)
        Dinv = np.zeros((n, n))
        done = np.zeros(n, bool)
        for key, members in self.block_members.items():
            idx = np.array(members)
            Dinv[np.ix_(idx, idx)] = np.linalg.inv(self.D[np.ix_(idx, idx)])
            done[idx] = True
        for i in np.where(~done)[0]:
            Dinv[i, i] = 1.0 / self.D[i, i]
        IM = np.eye(n) - self.M
        self.Ainv = IM.T @ Dinv @ IM
        return self.Ainv

    def index(self, ids):
        return np.array([self.inds[i].idx for i in ids])

    def sub(self, ids):
        ix = self.index(ids)
        return self.A[np.ix_(ix, ix)]


# ----------------------------------------------------------------------
# Comparateurs : getA (pedigreemm, diploïde classique) et calculate_bee_kinship (blup.r)
# ----------------------------------------------------------------------
def getA_diploid(ids, dam, sire):
    """Méthode tabulaire de Henderson (pedigreemm::getA). ids triés parents avant enfants.
    dam/sire : dict id -> parent ou None."""
    n = len(ids)
    pos = {k: i for i, k in enumerate(ids)}
    A = np.zeros((n, n))
    for i, k in enumerate(ids):
        d, s = dam.get(k), sire.get(k)
        di = pos[d] if d is not None else None
        si = pos[s] if s is not None else None
        for j in range(i):
            v = 0.0
            if di is not None:
                v += 0.5 * A[di, j]
            if si is not None:
                v += 0.5 * A[si, j]
            A[i, j] = A[j, i] = v
        A[i, i] = 1.0 + (0.5 * A[di, si] if (di is not None and si is not None) else 0.0)
    return A


def blup_r_kinship(ids, dam, sire, n_drones=12, has_drone_info=True):
    """Réplique ligne à ligne de calculate_bee_kinship (blup.r l.220-267).
    - has_drone_info FALSE : getA
    - sinon : A * (1 + 1/(4 n)), puis pour tout couple d'ids ayant la même mère (non NA),
      quel que soit le père : 0.25 + 0.5/n."""
    A = getA_diploid(ids, dam, sire)
    if not has_drone_info:
        return A
    Aadj = A * (1 + 1.0 / (4 * n_drones))
    pos = {k: i for i, k in enumerate(ids)}
    queen_ids = [k for k in ids if dam.get(k) is not None]
    for a in range(len(queen_ids) - 1):
        for b in range(a + 1, len(queen_ids)):
            if dam[queen_ids[a]] == dam[queen_ids[b]]:
                i, j = pos[queen_ids[a]], pos[queen_ids[b]]
                Aadj[i, j] = Aadj[j, i] = 0.25 + 0.5 / n_drones
    return Aadj
