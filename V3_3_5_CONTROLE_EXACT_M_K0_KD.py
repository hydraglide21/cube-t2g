#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V3.3.5 — Controle exact de l'involution S sur M, K0 et Kd.

Sources attendues (dans le repertoire courant, ou indiquees en arguments) :
  MATRICE_S_59_V3_3_3.csv
  matrices_completes_Eg_T2g_exactes.txt Matrices du modèle MATRICE_S

Les matrices mecaniques globales sont assemblees par blocs diag(Eg,T2g),
avec dimensions 20 + 39 = 59. Tous les calculs sont faits avec SymPy et des
rationnels exacts. Aucun flottant n'est utilise.

Exemple :
  python V3_3_5_CONTROLE_EXACT_M_K0_KD.py
  python V3_3_5_CONTROLE_EXACT_M_K0_KD.py --S chemin/S.csv --archive chemin/matrices.txt

Dependance : sympy (pip install sympy)
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path
from typing import List

import sympy as sp

N = 59


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def exact_zero(matrix: sp.MatrixBase) -> bool:
    return all(sp.cancel(value) == 0 for value in matrix)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read_s_matrix(path: Path) -> sp.Matrix:
    entries = {}
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"i", "j", "S_ij_exact"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                f"{path}: colonnes requises {sorted(required)} absentes. "
                f"Entete trouvee : {reader.fieldnames}"
            )
        for line_number, row in enumerate(reader, start=2):
            i, j = int(row["i"]), int(row["j"])
            if not (0 <= i < N and 0 <= j < N):
                raise ValueError(f"{path}, ligne {line_number}: indice hors limites {(i, j)}")
            key = (i, j)
            if key in entries:
                raise ValueError(f"{path}, ligne {line_number}: entree dupliquee {key}")
            value_text = (row["S_ij_exact"] or "").strip()
            if not value_text:
                raise ValueError(f"{path}, ligne {line_number}: coefficient vide")
            entries[key] = sp.Rational(value_text)

    expected = {(i, j) for i in range(N) for j in range(N)}
    if set(entries) != expected:
        missing = len(expected - set(entries))
        extra = len(set(entries) - expected)
        raise ValueError(f"Matrice S incomplete/incoherente : {missing} manquantes, {extra} inattendues")

    matrix = sp.zeros(N, N)
    for (i, j), value in entries.items():
        matrix[i, j] = value
    return matrix


def read_archived_matrix(archive_text: str, name: str, size: int) -> sp.Matrix:
    """Lit une matrice SymPy imprimee dans l'archive, sans eval/float."""
    header = re.compile(rf"^\s*{re.escape(name)}\s*\({size}x{size}\)\s*=\s*$", re.M)
    match = header.search(archive_text)
    if match is None:
        raise ValueError(f"Matrice {name} ({size}x{size}) introuvable dans l'archive")

    # Extraire le contenu de Matrix([ ... ]) puis lire chaque ligne [a,b,...].
    tail = archive_text[match.end():]
    matrix_start = re.search(r"Matrix\(\[", tail)
    if matrix_start is None:
        raise ValueError(f"Bloc Matrix([ absent apres l'entete {name}")
    tail = tail[matrix_start.end():]
    matrix_end = re.search(r"\]\s*\)", tail, flags=re.S)
    if matrix_end is None:
        raise ValueError(f"Fermeture Matrix introuvable pour {name}")
    body = tail[:matrix_end.start()]
    row_strings = re.findall(r"\[([^\[\]]*)\]", body, flags=re.S)
    rows: List[List[sp.Rational]] = []
    for row_number, row in enumerate(row_strings, start=1):
        tokens = [token.strip() for token in row.split(",")]
        if len(tokens) != size:
            raise ValueError(
                f"{name}: ligne {row_number}, {len(tokens)} coefficients au lieu de {size}"
            )
        try:
            rows.append([sp.Rational(token) for token in tokens])
        except Exception as exc:
            raise ValueError(f"{name}: coefficient rationnel invalide dans la ligne {row_number}") from exc
    if len(rows) != size:
        raise ValueError(f"{name}: {len(rows)} lignes lues au lieu de {size}")
    return sp.Matrix(rows)


def write_matrix(path: Path, matrix: sp.MatrixBase) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow([f"c{j}" for j in range(matrix.cols)])
        for i in range(matrix.rows):
            writer.writerow([str(matrix[i, j]) for j in range(matrix.cols)])


def write_basis(path: Path, basis: sp.MatrixBase) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["dof"] + [f"basis_{j+1}" for j in range(basis.cols)])
        for i in range(basis.rows):
            writer.writerow([i] + [str(basis[i, j]) for j in range(basis.cols)])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--S", type=Path, default=Path("MATRICE_S_59_V3_3_3.csv"))
    parser.add_argument(
        "--archive", type=Path,
        default=Path("matrices_completes_Eg_T2g_exactes.txt Matrices du modèle MATRICE_S"),
    )
    parser.add_argument("--out", type=Path, default=Path("V3_3_5_SORTIES"))
    args = parser.parse_args()

    for source in (args.S, args.archive):
        if not source.is_file():
            raise FileNotFoundError(f"Fichier source absent : {source}")
    args.out.mkdir(parents=True, exist_ok=True)

    S = read_s_matrix(args.S)
    archive_text = args.archive.read_text(encoding="utf-8")

    mechanical = {}
    for name in ("M", "K0", "Kd"):
        eg = read_archived_matrix(archive_text, name + "_E", 20)
        t2g = read_archived_matrix(archive_text, name + "_T", 39)
        mechanical[name] = sp.diag(eg, t2g)

    identity = sp.eye(N)
    require(exact_zero(S * S - identity), "ECHEC : S^2 != I_59")
    trace = sp.trace(S)
    require(trace == 13, f"ECHEC : trace(S) = {trace}, valeur attendue 13")

    projector_plus = (identity + S) / 2
    projector_minus = (identity - S) / 2
    require(exact_zero(projector_plus**2 - projector_plus), "P+ n'est pas un projecteur")
    require(exact_zero(projector_minus**2 - projector_minus), "P- n'est pas un projecteur")
    require(exact_zero(projector_plus * projector_minus), "P+ P- != 0")
    require(exact_zero(projector_plus + projector_minus - identity), "P+ + P- != I")

    basis_plus = sp.Matrix.hstack(*((S - identity).nullspace()))
    basis_minus = sp.Matrix.hstack(*((S + identity).nullspace()))
    require(basis_plus.cols == 36, f"dim(H+) inattendue : {basis_plus.cols}")
    require(basis_minus.cols == 23, f"dim(H-) inattendue : {basis_minus.cols}")
    require(exact_zero(S * basis_plus - basis_plus), "S U+ != U+")
    require(exact_zero(S * basis_minus + basis_minus), "S U- != -U-")
    adapted_basis = basis_plus.row_join(basis_minus)
    require(adapted_basis.rank() == N, "Rang([U+ U-]) != 59")

    report = [
        "V3.3.5 — CONTROLE EXACT DE L'INVOLUTION S",
        "Calculs : rationnels exacts SymPy ; aucun flottant.",
        f"Source S : {args.S} | SHA256={sha256(args.S)}",
        f"Archive matrices : {args.archive} | SHA256={sha256(args.archive)}",
        "",
        "INVOLUTION ET ESPACES PROPRES",
        "[DEMONTRÉ EXACT] S^2 = I_59.",
        f"[DEMONTRÉ EXACT] Tr(S) = {trace}.",
        "[DEMONTRÉ EXACT] P+ et P- sont des projecteurs complementaires.",
        f"[DEMONTRÉ EXACT] dim(H+) = {basis_plus.cols}.",
        f"[DEMONTRÉ EXACT] dim(H-) = {basis_minus.cols}.",
        "[DEMONTRÉ EXACT] rang([U+ U-]) = 59.",
        "",
        "CONTROLES DES MATRICES MECANIQUES",
    ]

    for name, matrix in mechanical.items():
        require(exact_zero(matrix - matrix.T), f"{name} n'est pas symetrique")
        require(exact_zero(S.T * matrix * S - matrix), f"S^T {name} S != {name}")
        cross = basis_plus.T * matrix * basis_minus
        require(exact_zero(cross), f"Bloc croise H+/H- non nul pour {name}")
        plus_block = basis_plus.T * matrix * basis_plus
        minus_block = basis_minus.T * matrix * basis_minus
        adapted = adapted_basis.T * matrix * adapted_basis
        require(exact_zero(adapted - sp.diag(plus_block, minus_block)),
                f"La matrice adaptee {name} n'est pas bloc-diagonale")

        write_matrix(args.out / f"POURSUITE_PLUS_{name}.csv", plus_block)
        write_matrix(args.out / f"POURSUITE_MOINS_{name}.csv", minus_block)
        write_matrix(args.out / f"BLOC_CROISE_{name}.csv", cross)
        write_matrix(args.out / f"MATRICE_ADAPTEE_{name}.csv", adapted)
        report.append(
            f"[DEMONTRÉ EXACT] {name} : symetrie, S^T {name} S = {name}, "
            "bloc croise nul et matrice adaptee bloc-diagonale."
        )

    # Dimensions par secteur, calculees independamment sur les deux blocs de S.
    s_eg, s_t2g = S[:20, :20], S[20:, 20:]
    require(exact_zero(S[:20, 20:]) and exact_zero(S[20:, :20]),
            "S n'est pas bloc-diagonale selon les secteurs Eg/T2g")
    eg_plus = (s_eg - sp.eye(20)).nullspace()
    eg_minus = (s_eg + sp.eye(20)).nullspace()
    t_plus = (s_t2g - sp.eye(39)).nullspace()
    t_minus = (s_t2g + sp.eye(39)).nullspace()
    require((len(eg_plus), len(eg_minus), len(t_plus), len(t_minus)) == (10, 10, 26, 13),
            "Dimensions propres sectorielles inattendues")

    write_basis(args.out / "BASE_H_PLUS_EXACTE.csv", basis_plus)
    write_basis(args.out / "BASE_H_MOINS_EXACTE.csv", basis_minus)
    with (args.out / "DIMENSIONS_H_PLUS_H_MOINS.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["secteur", "dimension", "H_plus", "H_moins"])
        writer.writerow(["Eg", 20, len(eg_plus), len(eg_minus)])
        writer.writerow(["T2g", 39, len(t_plus), len(t_minus)])
        writer.writerow(["Total", 59, basis_plus.cols, basis_minus.cols])

    report += [
        "",
        "DIMENSIONS PAR SECTEUR",
        f"Eg : H+={len(eg_plus)}, H-={len(eg_minus)}.",
        f"T2g : H+={len(t_plus)}, H-={len(t_minus)}.",
        "",
        "PORTEE ET LIMITES",
        "Le decouplage algebrique H+/H- est demontre ici pour M, K0 et Kd.",
        "Kedge n'est pas teste : aucune matrice source Kedge 59x59 n'est fournie.",
        "Ces controles ne demontrent pas la positivite, la stabilite dynamique,",
        "le confinement, les non-linearites, les solitons ou les breathers.",
    ]

    report_path = args.out / "RAPPORT_V3_3_5.txt"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    manifest_path = args.out / "MANIFESTE_SHA256_V3_3_5.txt"
    with manifest_path.open("w", encoding="utf-8") as stream:
        for path in sorted(p for p in args.out.iterdir() if p.is_file() and p != manifest_path):
            stream.write(f"{sha256(path)}  {path.name}\n")

    print("\n".join(report))
    print(f"\nFichiers produits dans : {args.out.resolve()}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ECHEC V3.3.5 : {exc}", file=sys.stderr)
        raise
