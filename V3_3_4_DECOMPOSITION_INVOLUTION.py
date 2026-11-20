#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V3.3.4 — DECOMPOSITION EXACTE DE L'INVOLUTION S SUR 59 DEGRES DE LIBERTE

But
---
Construire exactement les espaces propres H+ = ker(S-I) et H- = ker(S+I),
puis verifier que les matrices mecaniques invariantes se decomposent en deux
blocs independants dans cette base.

Aucune modification des fichiers V3.3.3.
Toutes les operations algebriques sont effectuees avec SymPy/Rationals.

ENTREES
-------
1) MATRICE_S_59_V3_3_3.csv
   Colonnes attendues : i,j,S_ij_exact
   (S_ij_flottant peut etre present mais n'est jamais utilise.)

2) Matrices M, K0, Kd, Kedge sous forme CSV :
   - format long : colonnes i,j,value_exact (ou exact/value/coefficient)
   - ou format dense : 59 lignes x 59 colonnes de coefficients rationnels.

Les matrices doivent utiliser le meme ordre de degres de liberte que S.
Les indices du format long peuvent commencer a 0 ou a 1.

EXEMPLE
-------
python V3_3_4_DECOMPOSITION_INVOLUTION.py \
  --S MATRICE_S_59_V3_3_3.csv \
  --M MATRICE_M_59.csv \
  --K0 MATRICE_K0_59.csv \
  --Kd MATRICE_Kd_59.csv \
  --Kedge MATRICE_Kedge_59.csv \
  --out V3_3_4_SORTIES

SORTIES
-------
- RAPPORT_V3_3_4.txt
- DIMENSIONS_H_PLUS_H_MOINS.csv
- BASE_H_PLUS_EXACTE.csv
- BASE_H_MOINS_EXACTE.csv
- POURSUITE_PLUS_<matrice>.csv
- POURSUITE_MOINS_<matrice>.csv
- BLOC_CROISE_<matrice>.csv
- MANIFESTE_SHA256_V3_3_4.txt

STATUT
------
Ce script est un protocole de calcul. Il ne certifie aucun resultat tant qu'il
n'a pas ete execute sur les fichiers sources authentiques et que ses sorties
n'ont pas ete examinees.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import sympy as sp

N = 59
S_EXACT_COL = "S_ij_exact"


def rational(value: str) -> sp.Rational:
    """Parse strictement un coefficient rationnel, sans conversion float."""
    value = str(value).strip()
    if value == "":
        raise ValueError("Coefficient rationnel vide.")
    if value.lower() in {"nan", "inf", "+inf", "-inf", "none"}:
        raise ValueError(f"Coefficient non fini/invalide : {value}")
    try:
        q = sp.Rational(value)
    except Exception as exc:
        raise ValueError(f"Coefficient non rationnel {value!r}") from exc
    if q.is_Rational is not True:
        raise ValueError(f"Coefficient non rationnel : {value!r}")
    return q


def read_csv_rows(path: Path) -> List[List[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.reader(f))


def read_s_matrix(path: Path) -> sp.Matrix:
    rows = read_csv_rows(path)
    if not rows:
        raise ValueError(f"Fichier vide : {path}")

    header = [x.strip() for x in rows[0]]
    required = {"i", "j", S_EXACT_COL}
    if not required.issubset(set(header)):
        raise ValueError(
            f"{path}: colonnes requises {sorted(required)}, trouvees {header}"
        )

    ix, jx, vx = header.index("i"), header.index("j"), header.index(S_EXACT_COL)
    entries: Dict[Tuple[int, int], sp.Rational] = {}
    raw_indices = []

    for line_no, row in enumerate(rows[1:], start=2):
        if not row or all(not x.strip() for x in row):
            continue
        try:
            i_raw, j_raw = int(row[ix]), int(row[jx])
            v = rational(row[vx])
        except Exception as exc:
            raise ValueError(f"{path}, ligne {line_no}: {exc}") from exc
        raw_indices.extend([i_raw, j_raw])
        key = (i_raw, j_raw)
        if key in entries:
            raise ValueError(f"{path}: entree dupliquee {key}")
        entries[key] = v

    if not raw_indices:
        raise ValueError(f"{path}: aucune entree.")

    base = 0 if min(raw_indices) == 0 else 1
    M = sp.zeros(N, N)
    expected = {(i + base, j + base) for i in range(N) for j in range(N)}
    if set(entries) != expected:
        missing = len(expected - set(entries))
        extra = len(set(entries) - expected)
        raise ValueError(
            f"{path}: matrice incomplete/incoherente; "
            f"{missing} entrees manquantes, {extra} inattendues."
        )
    for i in range(N):
        for j in range(N):
            M[i, j] = entries[(i + base, j + base)]
    return M


def read_matrix(path: Path) -> sp.Matrix:
    """Lit une matrice rationnelle au format long ou dense."""
    rows = read_csv_rows(path)
    if not rows:
        raise ValueError(f"Fichier vide : {path}")

    header = [x.strip() for x in rows[0]]
    long_value_candidates = ("value_exact", "exact", "coefficient", "value")
    value_col = next((c for c in long_value_candidates if c in header), None)

    if {"i", "j"}.issubset(set(header)) and value_col is not None:
        ix, jx, vx = header.index("i"), header.index("j"), header.index(value_col)
        entries: Dict[Tuple[int, int], sp.Rational] = {}
        raw_indices = []
        for line_no, row in enumerate(rows[1:], start=2):
            if not row or all(not x.strip() for x in row):
                continue
            try:
                i_raw, j_raw = int(row[ix]), int(row[jx])
                val = rational(row[vx])
            except Exception as exc:
                raise ValueError(f"{path}, ligne {line_no}: {exc}") from exc
            raw_indices.extend([i_raw, j_raw])
            key = (i_raw, j_raw)
            if key in entries:
                raise ValueError(f"{path}: entree dupliquee {key}")
            entries[key] = val

        base = 0 if min(raw_indices) == 0 else 1
        expected = {(i + base, j + base) for i in range(N) for j in range(N)}
        if set(entries) != expected:
            raise ValueError(
                f"{path}: matrice longue incomplete; "
                f"{len(expected-set(entries))} entrees manquantes."
            )
        A = sp.zeros(N, N)
        for i in range(N):
            for j in range(N):
                A[i, j] = entries[(i + base, j + base)]
        return A

    # Format dense : accepte un CSV sans en-tete, ou un en-tete de 59 colonnes.
    dense_rows = rows
    if len(dense_rows) == N + 1:
        # Une ligne d'en-tete est toleree uniquement si elle n'est pas numerique.
        try:
            [rational(x) for x in dense_rows[0]]
        except Exception:
            dense_rows = dense_rows[1:]

    if len(dense_rows) != N:
        raise ValueError(
            f"{path}: format non reconnu. Attendu 59x59 dense ou format long."
        )

    parsed = []
    for line_no, row in enumerate(dense_rows, start=1):
        if len(row) != N:
            raise ValueError(
                f"{path}, ligne {line_no}: {len(row)} colonnes au lieu de {N}."
            )
        parsed.append([rational(x) for x in row])
    return sp.Matrix(parsed)


def write_matrix(path: Path, A: sp.Matrix) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow([f"c{j}" for j in range(A.cols)])
        for i in range(A.rows):
            w.writerow([str(A[i, j]) for j in range(A.cols)])


def write_basis(path: Path, B: sp.Matrix) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["dof"] + [f"basis_{j+1}" for j in range(B.cols)])
        for i in range(B.rows):
            w.writerow([i] + [str(B[i, j]) for j in range(B.cols)])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def exact_zero(A: sp.Matrix) -> bool:
    return all(sp.cancel(x) == 0 for x in A)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--S", required=True, type=Path)
    parser.add_argument("--M", required=True, type=Path)
    parser.add_argument("--K0", required=True, type=Path)
    parser.add_argument("--Kd", required=True, type=Path)
    parser.add_argument("--Kedge", required=True, type=Path)
    parser.add_argument("--out", default=Path("V3_3_4_SORTIES"), type=Path)
    args = parser.parse_args()

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    report: List[str] = []

    def log(msg: str) -> None:
        print(msg)
        report.append(msg)

    log("V3.3.4 — DECOMPOSITION EXACTE DE L'INVOLUTION S")
    log(f"Dimension : {N}")
    log("Arithmetique : rationnels exacts SymPy")
    log("")

    inputs = {
        "S": args.S,
        "M": args.M,
        "K0": args.K0,
        "Kd": args.Kd,
        "Kedge": args.Kedge,
    }
    for name, path in inputs.items():
        require(path.is_file(), f"Fichier source absent : {path}")
        log(f"SOURCE {name}: {path} | SHA256={sha256(path)}")

    S = read_s_matrix(args.S)
    mats = {
        "M": read_matrix(args.M),
        "K0": read_matrix(args.K0),
        "Kd": read_matrix(args.Kd),
        "Kedge": read_matrix(args.Kedge),
    }

    I = sp.eye(N)
    log("")
    log("=== CONTROLES DE L'INVOLUTION ===")
    err_s2 = S * S - I
    require(exact_zero(err_s2), "ECHEC : S^2 != I.")
    log("[DEMONTRÉ EXACT] S^2 = I.")

    trS = sp.trace(S)
    require(trS.q == 1, "Trace de S non entiere.")
    dim_plus_formula = sp.Rational(N + trS, 2)
    dim_minus_formula = sp.Rational(N - trS, 2)
    require(dim_plus_formula.q == 1 and dim_minus_formula.q == 1,
            "Dimensions calculees non entieres.")

    Pp = (I + S) / 2
    Pm = (I - S) / 2
    require(exact_zero(Pp * Pp - Pp), "ECHEC : P+ non projecteur.")
    require(exact_zero(Pm * Pm - Pm), "ECHEC : P- non projecteur.")
    require(exact_zero(Pp * Pm), "ECHEC : P+P- != 0.")
    require(exact_zero(Pp + Pm - I), "ECHEC : P++P- != I.")
    log("[DEMONTRÉ EXACT] P+ et P- sont des projecteurs complementaires.")
    log(f"Trace(S) = {trS}")
    log(f"Dimension attendue H+ = {dim_plus_formula}")
    log(f"Dimension attendue H- = {dim_minus_formula}")

    # Les colonnes de nullspace forment des bases rationnelles exactes.
    Bp = (S - I).nullspace()
    Bm = (S + I).nullspace()
    Uplus = sp.Matrix.hstack(*Bp) if Bp else sp.zeros(N, 0)
    Uminus = sp.Matrix.hstack(*Bm) if Bm else sp.zeros(N, 0)
    dp, dm = Uplus.cols, Uminus.cols

    require(dp == int(dim_plus_formula), "Dimension de ker(S-I) incoherente.")
    require(dm == int(dim_minus_formula), "Dimension de ker(S+I) incoherente.")
    require(dp + dm == N, "Les deux sous-espaces ne couvrent pas R^59.")
    require(exact_zero(S * Uplus - Uplus), "ECHEC : S U+ != U+.")
    require(exact_zero(S * Uminus + Uminus), "ECHEC : S U- != -U-.")
    T = Uplus.row_join(Uminus)
    require(T.rank() == N, "La concaténation des bases n'est pas inversible.")

    log("[DEMONTRÉ EXACT] H+ = ker(S-I), H- = ker(S+I).")
    log(f"Dimension calculee H+ = {dp}")
    log(f"Dimension calculee H- = {dm}")
    log("[DEMONTRÉ EXACT] Rang([U+ U-]) = 59.")
    write_basis(out / "BASE_H_PLUS_EXACTE.csv", Uplus)
    write_basis(out / "BASE_H_MOINS_EXACTE.csv", Uminus)

    with (out / "DIMENSIONS_H_PLUS_H_MOINS.csv").open(
        "w", encoding="utf-8", newline=""
    ) as f:
        w = csv.writer(f)
        w.writerow(["dimension_totale", "trace_S", "dim_H_plus", "dim_H_moins"])
        w.writerow([N, str(trS), dp, dm])

    log("")
    log("=== CONTROLES DES MATRICES MECANIQUES ===")
    for name, A in mats.items():
        require(A.shape == (N, N), f"{name}: dimensions incorrectes.")
        require(exact_zero(A - A.T), f"{name}: matrice non symetrique.")
        inv_err = S.T * A * S - A
        require(exact_zero(inv_err), f"{name}: S^T {name} S != {name}.")
        cross = Uplus.T * A * Uminus
        require(exact_zero(cross), f"{name}: bloc croise H+/H- non nul.")
        plus = Uplus.T * A * Uplus
        minus = Uminus.T * A * Uminus
        write_matrix(out / f"POURSUITE_PLUS_{name}.csv", plus)
        write_matrix(out / f"POURSUITE_MOINS_{name}.csv", minus)
        write_matrix(out / f"BLOC_CROISE_{name}.csv", cross)
        log(f"[DEMONTRÉ EXACT] {name}: symetrie, invariance et bloc croise nul.")

    # La base concatenee doit rendre explicite la structure par blocs.
    log("")
    log("=== MATRICES DANS LA BASE ADAPTEE ===")
    for name, A in mats.items():
        adapted = T.T * A * T
        expected = sp.diag(
            Uplus.T * A * Uplus,
            Uminus.T * A * Uminus
        )
        require(exact_zero(adapted - expected),
                f"{name}: la matrice adaptee n'est pas bloc diagonale.")
        write_matrix(out / f"MATRICE_ADAPTEE_{name}.csv", adapted)
        log(f"[DEMONTRÉ EXACT] {name}: T^T {name} T est bloc diagonale.")

    # Rapport synthétique, puis manifeste des fichiers produits.
    log("")
    log("=== CONCLUSION ===")
    log("La decomposition H+ / H- est verifiee exactement pour les fichiers fournis.")
    log("Les matrices mecaniques M, K0, Kd et Kedge sont bloc-diagonales")
    log("dans la base adaptee, avec blocs croises exactement nuls.")
    log("Cela etablit une decouplage algebrique pour ces matrices invariantes.")
    log("Cela ne demontre pas a lui seul la positivite, la stabilite dynamique,")
    log("la non-linearite, le confinement, les solitons ou les breathers.")
    log("Les faces individuelles y/z ne sont pas supposees invariantes : elles")
    log("sont permutees par S et doivent etre traitees comme une famille covariante.")

    report_path = out / "RAPPORT_V3_3_4.txt"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")

    manifest_path = out / "MANIFESTE_SHA256_V3_3_4.txt"
    files = sorted(p for p in out.iterdir() if p.is_file() and p != manifest_path)
    with manifest_path.open("w", encoding="utf-8") as f:
        for p in files:
            f.write(f"{sha256(p)}  {p.name}\n")

    print(f"\nSorties ecrites dans : {out.resolve()}")
    print(f"Manifeste : {manifest_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ECHEC V3.3.4 : {exc}", file=sys.stderr)
        raise
