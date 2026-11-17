```
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
V17.5 — THEOREME DE NON-CONFINEMENT DYNAMIQUE DU MODELE REDUIT.
Sources (lecture seule) : V16.1 + V17.4, SHA verifies avant toute operation.
M = M13, K = K0_13 + (1/2) Kd_13, A = M^{-1} K (resolution exacte M X = K).
O = empilement des blocs (A^m)_LN, m = 0..12 (39x10 ; le bloc m = 0 est nul, donc
la version m = 1..12, 36x10, a meme rang et meme noyau).
Tout calcul en Q, arithmetique rationnelle exacte. Aucun fichier anterieur modifie.
"""
import hashlib, json, pathlib, pickle, re, time
import numpy as np
import sympy as sp

HERE = pathlib.Path('/mnt/data/V31_repro')
V161, V171, V172, V173, V174 = HERE/'V16_1', HERE/'V17_1', HERE/'V17_2', HERE/'V17_3', HERE/'V17_4'
OUT = HERE/'V17_5'; OUT.mkdir(exist_ok=True)
LOG = []
def note(s):
    LOG.append(s); print(s, flush=True)
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
t0 = time.time()

# ---------- 0. SHA-256 des sources, AVANT toute operation ----------
EXP = {
 'V16_1_MATRICES_REDUITES.pkl': '43f620d47a1535d5a474a67add1e269c079942873d6f7f1ef71108c57bf80d07',
 'V16_1_RAPPORT.json': '76adc689c2764512d4e58e7c47b423ebdfc81329d58a8cb03a6e734f66152154',
 'V17_4_SCRIPT.py': '6025fa5bfe829026819d47abb49cdd6d73a3abddf93907544233eac8567912fc',
 'V17_4_CERTIFICAT_EXACT.txt': '565ae2eaac0f627f7c6a7f3fe316dc71683ec940adf8fdb5a08bc159e2519ec7',
 'V17_4_RAPPORT.json': 'a0eb732b1b8d1ac047b25a237e2fdab95d1cb429613f88649c3ee3e6955494d7',
}
SRC_PATHS = {
 'V16_1_MATRICES_REDUITES.pkl': V161/'V16_1_MATRICES_REDUITES.pkl',
 'V16_1_RAPPORT.json': V161/'V16_1_RAPPORT.json',
 'V17_4_SCRIPT.py': V174/'V17_4_SCRIPT.py',
 'V17_4_CERTIFICAT_EXACT.txt': V174/'V17_4_CERTIFICAT_EXACT.txt',
 'V17_4_RAPPORT.json': V174/'V17_4_RAPPORT.json',
}
src_sha = {k: sha(p) for k, p in SRC_PATHS.items()}
for k, v in src_sha.items():
    ok = (v == EXP[k]); note(f"[SHA] {k} : {v} -> {'OK' if ok else 'ECHEC'}")
    assert ok, f'SHA mismatch {k}'
# Journal V17.4 : aucun manifeste anterieur n'existe pour lui ; on ne peut pas le contre-verifier
# contre un enregistrement prealable, donc on n'invente pas de SHA attendu. On verifie en revanche
# que les SHA inscrits DANS le journal correspondent aux SHA actuels des fichiers V17.4 qu'il cite.
j174 = (V174/'V17_4_JOURNAL.txt').read_text(encoding='utf-8')
sha_j174_now = sha(V174/'V17_4_JOURNAL.txt')
inside = {}
for m in re.finditer(r'V17_4_(SCRIPT\.py|CERTIFICAT_EXACT\.txt|RAPPORT\.json)\s+([0-9a-f]{64})', j174):
    inside['V17_4_'+m.group(1)] = m.group(2)
jcheck = {k: (inside.get(k) == src_sha.get(k)) for k in ('V17_4_SCRIPT.py', 'V17_4_CERTIFICAT_EXACT.txt', 'V17_4_RAPPORT.json')}
note(f"[SHA] Journal V17.4 actuel (sans manifeste anterieur, enregistre a titre de tracabilite) : {sha_j174_now}")
note(f"[SHA] SHA inscrits DANS le journal V17.4 == SHA actuels : {jcheck}")
assert all(jcheck.values())
# Empreintes "avant" de tous les fichiers anterieurs (pour le controle "inchangé" final)
OLD_FILES = {}
for d in (V161, V171, V172, V173, V174):
    for p in sorted(d.glob('*')):
        if p.is_file(): OLD_FILES[str(p.relative_to(HERE))] = sha(p)
note(f"[SHA] Empreintes avant calcul enregistrees pour {len(OLD_FILES)} fichiers anterieurs.")

# ---------- 1. Matrices exactes, operateur A ----------
data = pickle.loads(SRC_PATHS['V16_1_MATRICES_REDUITES.pkl'].read_bytes())
M = sp.Matrix(data['M13']); K0 = sp.Matrix(data['K0_13']); Kd = sp.Matrix(data['Kd_13'])
K = K0 + sp.Rational(1, 2)*Kd
assert M == M.T and K == K.T, 'symetrie'
A = M.solve(K)  # resolution exacte M X = K (M inversible : la resolution reussit)
assert M*A == K
ALN = A[10:13, 0:10]
note(f"[MATRICE] M, K symetriques exactes ; A = M^-1 K par resolution exacte ; M*A == K : True ; A_LN : 3x10.")

# ---------- 2. Contre-verification V17.4 ----------
rep174 = json.loads(SRC_PATHS['V17_4_RAPPORT.json'].read_text())
aln174 = [[sp.Rational(e) for e in row] for row in rep174['A_LN']['coefficients_exacts']]
assert sp.Matrix(aln174) == ALN, 'A_LN != V17.4 (JSON)'
cert174 = SRC_PATHS['V17_4_CERTIFICAT_EXACT.txt'].read_text(encoding='utf-8')
sec = cert174.split('2. BLOC A_LN')[1].split('Rang exact')[0]
rows = [sp.Matrix([[sp.Rational(x.strip()) for x in ln.strip().strip('[]').split(',')]]) for ln in sec.splitlines() if ln.strip().startswith('[')]
assert sp.Matrix.vstack(*rows) == ALN, 'A_LN != V17.4 (certificat txt)'
note('[XCHECK] A_LN identique au certificat V17.4 (JSON + txt) : True')
# Ordres elementaires V17.4
powers = [sp.eye(13)]
for m in range(1, 13): powers.append(powers[-1]*A)
first_order = {}
for i in range(10):
    mf = None
    for m in range(1, 13):
        if powers[m][10:13, 0:10][:, i] != sp.zeros(3, 1): mf = m; break
    first_order[i] = mf
ord174 = {int(k): v for k, v in rep174['ordres']['premier_ordre_non_nul_par_i'].items()}
note(f"[XCHECK] Ordres elementaires == V17.4 : {first_order == ord174} (recalcules : {first_order})")
assert first_order == ord174

# ---------- 3. Matrice d'observabilite O (m = 0..12, 39x10) ----------
O = sp.Matrix.vstack(*[powers[m][10:13, 0:10] for m in range(13)])
assert O[0:3, :] == sp.zeros(3, 10), 'le bloc m=0 doit etre nul (I_LN = 0)'
O36 = O[3:, :]  # version m = 1..12 (36x10) : meme rang, meme noyau (le bloc m=0 est nul)
rank_O = O.rank(); ker_O = O.nullspace()
rank_O36 = O36.rank(); ker_O36 = O36.nullspace()
note(f"[O] 39x10, rang exact = {rank_O}, dim ker O = {len(ker_O)} ; version 36x10 (m=1..12) : rang = {rank_O36}, dim ker = {len(ker_O36)}")
assert rank_O == 10 and len(ker_O) == 0 and rank_O36 == 10 and len(ker_O36) == 0

# ---------- 4. Certificat de rang : mineur 10x10 non nul (lignes pivots du RREF de O^T) ----------
_, pivots = O.T.rref()
pivots = [int(p) for p in pivots]
assert len(pivots) == 10
minor = O[pivots, :]
note(f"[MINEUR] Lignes retenues (parmi les 39) : {pivots}  (bloc m = ligne//3, ligne dans le bloc = ligne%3)")
det_bareiss = minor.det()
det_domain = minor.det(method='domain-ge')
inv_check = minor*minor.inv() == sp.eye(10)
note(f"[MINEUR] det (Bareiss) == det (domain-ge) : {det_bareiss == det_domain} ; minor * minor^-1 == I : {inv_check}")
assert det_bareiss != 0 and det_bareiss == det_domain and inv_check
det_minor = det_bareiss
note(f"[MINEUR] Determinant exact : {str(det_minor)[:120]}... (longueur : {len(str(det_minor))} caracteres)")

# ---------- 5. Cayley-Hamilton : p(A) = 0 (Horner, exact) ----------
lam = sp.symbols('lam'); pch = A.charpoly(lam); coeffs = pch.all_coeffs()
pA = sp.zeros(13, 13)
for c in coeffs: pA = pA*A + c*sp.eye(13)
ch_ok = (pA == sp.zeros(13, 13))
note(f"[CAYLEY-HAMILTON] p(A) = 0 verifie exactement : {ch_ok}")
assert ch_ok

# ---------- 6. Illustrations exactes (non des preuves : le theoreme est ker O = {0} + Cayley-Hamilton) ----------
# e_0 : activation differee, premier ordre m = 2
lead_e0 = powers[2][10:13, 0:10][:, 0]
# v = (1,...,1) : combinaison generique
v_gen = sp.ones(10, 1)
m_gen = next(m for m in range(1, 13) if powers[m][10:13, 0:10]*v_gen != sp.zeros(3, 1))
lead_gen = powers[m_gen][10:13, 0:10]*v_gen
note(f"[ILLUSTRATION] e_0 : premier ordre m = 2, coefficient dominant (A^2)_LN e_0 = {[str(e) for e in lead_e0]}")
note(f"[ILLUSTRATION] v = (1,...,1) : premier ordre m = {m_gen}, coefficient dominant = {[str(e) for e in lead_gen]}")

# ---------- 7. Verification "inchangé" des fichiers anterieurs ----------
after = {}
for d in (V161, V171, V172, V173, V174):
    for p in sorted(d.glob('*')):
        if p.is_file(): after[str(p.relative_to(HERE))] = sha(p)
unchanged = {k: (after.get(k) == v) for k, v in OLD_FILES.items()}
note(f"[SHA-APRES] Tous les fichiers anterieurs inchanges : {all(unchanged.values())} ({len(unchanged)} fichiers)")
assert all(unchanged.values())

# ---------- 8. Livrables ----------
theo = []
theo.append('V17.5 — THEOREME DE NON-CONFINEMENT DYNAMIQUE DU MODELE REDUIT')
theo.append('='*80)
theo.append('')
theo.append('ENONCE')
theo.append('Soient M = M13 et K = K0_13 + (1/2) Kd_13 les matrices reduites exactes V16.1')
theo.append('(13x13, rationnelles), A = M^{-1} K, partition N = {0,...,9}, L = {10,11,12}.')
theo.append('')
theo.append('THEOREME. A r = 1/2, dans le modele lineaire reduit M q\'\'(t) + K q(t) = 0, pour toute')
theo.append('condition initiale q(0) = (v, 0) avec v dans Q^10 non nul, et q\'(0) = 0, la composante')
theo.append('q_L(t) n\'est pas identiquement nulle ; en particulier la trajectoire ne reste pas')
theo.append('confinee dans N. Plus precisement, il existe m dans {1,...,12} tel que')
theo.append('(A^m)_LN v != 0, et')
theo.append('    q_L(t) = (-1)^m t^(2m)/(2m)! * (A^m)_LN v + O(t^(2m+2)).')
theo.append('')
theo.append('DEMONSTRATION')
theo.append('1. M est inversible (resolution exacte M X = K reussie ; verification M A = K).')
theo.append('   L\'equation s\'ecrit q\'\' = -A q.')
theo.append('2. Par recurrence sur m : q^(2m) = (-A)^m q et q^(2m+1) = (-A)^m q\'.')
theo.append('   En t = 0 avec q\'(0) = 0 : q^(2m)(0) = (-1)^m A^m q(0) et q^(2m+1)(0) = 0.')
theo.append('3. Les solutions d\'une EDO lineaire a coefficients constants sont analytiques (entieres).')
theo.append('   Serie de Taylor (convergente pour tout t) :')
theo.append('       q(t) = somme_{m>=0} t^(2m)/(2m)! * (-1)^m A^m q(0).')
theo.append('   Projection sur L :')
theo.append('       q_L(t) = somme_{m>=0} t^(2m)/(2m)! * (-1)^m (A^m)_LN v.')
theo.append('   Le terme m = 0 est nul car (A^0)_LN = I_LN = 0.')
theo.append('4. Cayley-Hamilton : p(A) = 0 ou p est le polynome caracteristique de A')
theo.append('   (verifie exactement par evaluation de Horner, voir certificat).')
theo.append('   Donc pour m >= 13, A^m est combinaison Q-lineaire de I, A, ..., A^12.')
theo.append('   Par consequent : si (A^k)_LN v = 0 pour k = 0,...,12, alors (A^m)_LN v = 0 pour tout m >= 0.')
theo.append('   Reciproquement : si (A^m)_LN v = 0 pour tout m >= 0, alors en particulier pour m = 0..12,')
theo.append('   donc O v = 0, ou O est la matrice d\'observabilite ci-dessous.')
theo.append('5. O est la matrice 39x10 formee des blocs (A^m)_LN, m = 0..12 (le bloc m = 0 est nul ;')
theo.append('   la version 36x10, m = 1..12, a meme rang et meme noyau).')
theo.append('   Certificat de rang (voir V17_5_CERTIFICAT_MINEUR.txt) : le mineur 10x10 forme des lignes')
theo.append(f'   {pivots} a pour determinant')
theo.append(f'   {det_minor}')
theo.append('   qui est non nul. Donc rang(O) >= 10 ; comme O a 10 colonnes, rang(O) = 10, donc ker(O) = {0}.')
theo.append('6. Pour v != 0, v n\'est pas dans ker(O), donc il existe m dans {0,...,12} avec')
theo.append('   (A^m)_LN v != 0. Le bloc m = 0 etant nul, m appartient a {1,...,12}.')
theo.append('   (m = 0 est exclu de la borne finale : le bloc identite a un bloc LN nul.)')
theo.append('7. Soit m* le plus petit tel m. Alors')
theo.append('       q_L(t) = (-1)^m* t^(2m*)/(2m*)! * (A^m*)_LN v + O(t^(2m*+2))')
theo.append('   avec un coefficient dominant non nul. Par analyticite, q_L(t) != 0 pour tout t > 0')
theo.append('   suffisamment petit. Donc q_L n\'est pas identiquement nulle : la trajectoire n\'est pas')
theo.append('   confinee dans N.  C.Q.F.D.')
theo.append('')
theo.append('REMARQUES')
theo.append('- Le resultat vaut pour tout v dans Q^10 non nul. L\'argument de rang etant independant du')
theo.append('  corps de base (les matrices sont rationnelles), il vaut aussi pour v dans R^10 non nul.')
theo.append('- Le theoreme n\'affirme PAS que q_L(t) != 0 a tous les instants ulterieurs : des annulations')
theo.append('  ponctuelles sont possibles. Il affirme la non-confinement : q_L n\'est pas identiquement nulle.')
theo.append('- Le theoreme porte sur les conditions initiales q(0) = (v,0), q\'(0) = 0.')
theo.append('  (Extension a q\'(0) != 0 : hors perimetre ; elle demanderait l\'observabilite complete')
theo.append('  incluant les derivees impaires.)')
theo.append('- Portee : exact sur Q ; dix coordonnees N, trois coordonnees L ; universel pour toute')
theo.append('  combinaison non nulle de N (pas seulement les dix excitations elementaires).')
theo.append('- Hors perimetre : cube physique complet, non-linearites, amortissement, autres valeurs de r,')
theo.append('  interpretation physique de resonance ou de soliton.')
theo.append('')
theo.append('ILLUSTRATIONS EXACTES (non des preuves : le theoreme est ker O = {0} + Cayley-Hamilton)')
theo.append(f'- e_0 (activation differee) : premier ordre m = 2, coefficient dominant (A^2)_LN e_0 = {[str(e) for e in lead_e0]}.')
theo.append(f'  Donc q_L(t) = t^4/24 * (A^2)_LN e_0 + O(t^6) pour la condition initiale e_0.')
theo.append(f'- v = (1,...,1) : premier ordre m = {m_gen}, coefficient dominant = {[str(e) for e in lead_gen]}.')
theo.append(f'  Donc q_L(t) = (-1)^{m_gen} t^{2*m_gen}/(2*{m_gen})! * (A^{m_gen})_LN v + O(t^{2*m_gen+2}).')
theo.append('')
theo.append('CONTROLES INFORMATIQUES')
theo.append('- SHA-256 des sources verifies avant toute operation (voir rapport JSON).')
theo.append('- M A = K verifie exactement ; symetrie exacte de M et K ; partition 10+3 reconstituee.')
theo.append('- A_LN identique au certificat V17.4 (JSON + txt) ; dix ordres elementaires coherents.')
theo.append('- Cayley-Hamilton p(A) = 0 verifie exactement (Horner).')
theo.append('- Determinant du mineur verifie par deux methodes (Bareiss, domain-ge) + minor * minor^-1 = I.')
theo.append('- Concordance numerique V17.3 : validation complementaire, NON une preuve du theoreme.')
theo.append('- Fichiers anterieurs verifies inchanges apres calcul (voir journal).')
(OUT/'V17_5_THEOREME_NON_CONFINEMENT.txt').write_text('\n'.join(theo)+'\n', encoding='utf-8')

certm = []
certm.append('V17.5 — CERTIFICAT DU MINEUR 10x10 DE LA MATRICE D\'OBSERVABILITE O')
certm.append('='*80)
certm.append('')
certm.append('Matrice d\'observabilite O : 39x10, blocs (A^m)_LN, m = 0..12.')
certm.append('Le bloc m = 0 est la matrice nulle (I_LN = 0) ; la version 36x10 (m = 1..12) a meme rang et meme noyau.')
certm.append('')
certm.append('Lignes retenues (parmi les 39 ; ligne = 3*m + j, m = ordre, j = 0..2 composante L) :')
certm.append(f'  {pivots}')
certm.append('Selection : lignes pivots du RREF exact de O^T (une base de lignes de O, ordre croissant).')
certm.append('')
certm.append('Matrice 10x10 correspondante (coefficients exacts, lignes dans l\'ordre retenu) :')
for r in range(10):
    certm.append('  [' + ', '.join(str(minor[r, c]) for c in range(10)) + ']')
certm.append('')
certm.append('Determinant rationnel exact (sans troncature) :')
certm.append(f'  {det_minor}')
certm.append('')
certm.append('Verification independante :')
certm.append(f'  - Methode de Bareiss : {det_bareiss}')
certm.append(f'  - Methode domain-ge  : {det_domain}')
certm.append(f'  - Egalite des deux methodes : {det_bareiss == det_domain}')
certm.append(f'  - minor * minor^-1 == I_10 : {inv_check}')
certm.append('  - Le determinant est non nul, donc les 10 lignes sont independantes : rang(O) >= 10.')
certm.append('  - O possede 10 colonnes, donc rang(O) <= 10. Conclusion : rang(O) = 10, ker(O) = {0}.')
certm.append('')
certm.append('Borne 12 (Cayley-Hamilton) : p(A) = 0 verifie exactement ; pour m >= 13, A^m est combinaison')
certm.append('Q-lineaire de I, A, ..., A^12. Donc O v = 0 implique (A^m)_LN v = 0 pour tout m >= 0,')
certm.append('et ker(O) = {0} implique que pour tout v != 0, il existe m dans {1,...,12} avec (A^m)_LN v != 0.')
(OUT/'V17_5_CERTIFICAT_MINEUR.txt').write_text('\n'.join(certm)+'\n', encoding='utf-8')

report = {
 'titre': 'V17.5 — Theoreme de non-confinement dynamique du modele reduit',
 'methode': 'Arithmetique rationnelle exacte (SymPy). A = M^-1 K par resolution exacte M X = K. O = empilement des blocs (A^m)_LN, m = 0..12 (39x10). Mineur 10x10 par lignes pivots du RREF de O^T. Cayley-Hamilton verifie par Horner. Aucun calcul numerique V17.1-V17.3 refait.',
 'sources_sha256': src_sha, 'sources_attendues': EXP, 'verification_sources': 'OK',
 'journal_V17_4': {'sha256_actuel': sha_j174_now, 'manifeste_anterieur': 'AUCUN (pas de manifeste V17.4) — SHA enregistre a titre de tracabilite, non contre-verifie contre un enregistrement prealable',
   'sha_inscrits_dans_journal_egaux_actuels': jcheck},
 'definitions': {'M': 'M13 (V16.1)', 'K': 'K0_13 + (1/2) Kd_13', 'A': 'M^-1 K (resolution exacte)', 'partition': 'N = indices 0..9, L = indices 10..12'},
 'controles': {'M_fois_A_egale_K': True, 'M_symetrique_exact': True, 'K_symetrique_exact': True, 'partition_10_3': True,
   'A_LN_identique_V17_4_JSON_txt': True, 'ordres_elementaires_coherents_V17_4': True},
 'observabilite': {'dimensions_O': '39x10 (blocs m = 0..12 ; bloc m = 0 nul)', 'dimensions_O_m_1_12': '36x10 (meme rang, meme noyau)',
   'rang_exact_O': rank_O, 'ker_O': '{0}', 'rang_O_m_1_12': rank_O36, 'ker_O_m_1_12': '{0}'},
 'mineur': {'lignes_retenues_parmi_39': pivots,
   'matrice_10x10_exacte': [[str(minor[r, c]) for c in range(10)] for r in range(10)],
   'determinant_exact': str(det_minor),
   'verification': {'bareiss': str(det_bareiss), 'domain_ge': str(det_domain), 'egalite': bool(det_bareiss == det_domain),
     'minor_fois_inverse_egale_I': bool(inv_check)},
   'conclusion_rang': 'det != 0 => rang(O) >= 10 ; O a 10 colonnes => rang(O) = 10, ker(O) = {0}'},
 'cayley_hamilton': {'p_A_egale_0_verifie_exactement': ch_ok, 'consequence': 'Pour m >= 13, A^m est combinaison Q-lineaire de I, A, ..., A^12.'},
 'theoreme': 'A r = 1/2, dans le modele lineaire reduit de dimension 13 defini par les matrices V16.1, aucune condition initiale non nulle portee par N (avec q\'(0) = 0) ne reste confinee dans N pour toute l\'evolution.',
 'demonstration': {'etapes': ['M inversible, q\'\' = -A q', 'q^(2m)(0) = (-1)^m A^m q(0), q^(2m+1)(0) = 0', 'Solutions analytiques (entieres) : serie de Taylor en puissances paires',
   'Cayley-Hamilton : A^m (m>=13) combinaison de I..A^12', 'O v = 0 <=> (A^m)_LN v = 0 pour tout m >= 0',
   'ker(O) = {0} (mineur) => pour tout v != 0, il existe m dans {1,...,12} avec (A^m)_LN v != 0',
   'Premier ordre non nul m* : q_L(t) = (-1)^m* t^(2m*)/(2m*)! (A^m*)_LN v + O(t^(2m*+2)) != 0 pour t > 0 petit'],
   'remarque': 'q_L(t) n\'est pas identiquement nul ; des annulations ponctuelles ulterieures restent possibles.'},
 'illustrations_exactes': {'e_0': {'premier_ordre': 2, 'coefficient_dominant': [str(e) for e in lead_e0]},
   'v_ones': {'premier_ordre': int(m_gen), 'coefficient_dominant': [str(e) for e in lead_gen]}},
 'comparaison_V17_3': {'statut': 'Validation complementaire uniquement, NON une preuve du theoreme',
   'rappel': 'V17.3 : les 10 essais elementaires activent L (E_L^max > seuil), y compris e_0 (activation differee, ordre 2). Coherent avec le theoreme.'},
 'contre_verification_V17_4': {'A_LN_identique': True, 'ordres_coherents': True, 'observabilite_borne_12_coherentes': True, 'resultats_anterieurs_modifies': False},
 'conclusions': {
   '1_rang_exact_O': 'OUI, rang exact de O = 10.',
   '2_mineur_retenu': f'Lignes {pivots} (parmi les 39).',
   '3_determinant_exact': str(det_minor),
   '4_cayley_hamilton_justifie_borne_12': 'OUI, p(A) = 0 verifie exactement ; A^m (m>=13) combinaison de I..A^12.',
   '5_theoreme_universel_demonte': 'OUI (mineur exact + Cayley-Hamilton + Taylor).',
   '6_limites_physiques_ouvertes': ['Cube physique complet', 'Non-linearites', 'Amortissement', 'Autres valeurs de r', 'Interpretation physique de resonance ou de soliton'],
 },
 'separation': {'demonstration_mathematique_exacte': 'Theoreme + preuve (fichier theoreme) ; certificat du mineur (fichier certificat).',
   'controles_informatiques': 'SHA, M*A = K, symetrie, A_LN == V17.4, determinant par deux methodes + inverse, fichiers anterieurs inchanges.',
   'comparaison_numerique_historique': 'V17.3 rappelee comme validation complementaire, non comme preuve.',
   'interpretation_physique': 'Hors perimetre.'},
 'fichiers_anterieurs_inchanges': unchanged,
 'livrables_sha256': {'V17_5_SCRIPT.py': sha(__file__),
   'V17_5_THEOREME_NON_CONFINEMENT.txt': sha(OUT/'V17_5_THEOREME_NON_CONFINEMENT.txt'),
   'V17_5_CERTIFICAT_MINEUR.txt': sha(OUT/'V17_5_CERTIFICAT_MINEUR.txt')},
 'duree_s': round(time.time()-t0, 1),
}
(OUT/'V17_5_RAPPORT.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
(OUT/'V17_5_JOURNAL.txt').write_text('\n'.join(LOG) + '\n' +
  f"[SHA] Livrables : V17_5_SCRIPT.py {sha(__file__)} ; V17_5_THEOREME_NON_CONFINEMENT.txt {sha(OUT/'V17_5_THEOREME_NON_CONFINEMENT.txt')} ; "
  f"V17_5_CERTIFICAT_MINEUR.txt {sha(OUT/'V17_5_CERTIFICAT_MINEUR.txt')} ; V17_5_RAPPORT.json {sha(OUT/'V17_5_RAPPORT.json')} ; "
  f"V17_5_JOURNAL.txt (en cours d'ecriture).\n", encoding='utf-8')
note('Fichiers : V17_5_SCRIPT.py | V17_5_RAPPORT.json | V17_5_THEOREME_NON_CONFINEMENT.txt | V17_5_CERTIFICAT_MINEUR.txt | V17_5_JOURNAL.txt')
```
