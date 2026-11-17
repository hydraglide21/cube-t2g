# CUBE-T2G — Géométrie cubique et dynamique

## Présentation

Ce dépôt rassemble les travaux consacrés à l'étude d'une structure géométrique cubique, de ses symétries et de ses propriétés spectrales et dynamiques.

Le programme de recherche suit un fil directeur :

**Cube → théorème QXY32_A (Will Hunting) → icosaèdre → dodécaèdre → dynamique du cube élastique.**

L'objectif est d'étudier les relations entre géométrie, couplages mécaniques, observabilité et dynamique des systèmes.

## Dossier V17.5 — Étude de non-confinement

V17.5 porte sur un modèle linéaire réduit à 13 degrés de liberté.

La partition étudiée comprend :
- 10 coordonnées initiales, notées N ;
- 3 coordonnées observées, notées L.

La matrice d'évolution est définie par :

\[
A=M^{-1}K,\qquad K=K_0+\frac12K_d.
\]

L'étude porte sur la possibilité qu'une condition initiale non nulle dans N ne produise jamais de composante observable dans L.

## Résultats annoncés dans les artefacts V17.5

- Matrice d'observabilité : 39 × 10.
- Rang exact annoncé : 10.
- Noyau annoncé : réduit à zéro.
- Mineur carré : 10 × 10.
- Déterminant exact annoncé :

\[
-108252777380873936258305659923321246663647119739453440000
\]

- Vérification de Cayley-Hamilton annoncée.

Ces résultats sont rapportés dans les livrables V17.5 et devront être reproduits à partir des sources et scripts correspondants avant d'être considérés comme indépendamment vérifiés.

## Livrables V17.5

1. `V17_5_SCRIPT.py` — script principal de calcul.
2. `V17_5_THEOREME_NON_CONFINEMENT.txt` — énoncé et conclusion.
3. `V17_5_CERTIFICAT_MINEUR.txt` — certificat du mineur.
4. `V17_5_RAPPORT.json` — rapport structuré.
5. `V17_5_JOURNAL.txt` — journal des opérations.
6. `V17_5_MANIFESTE_SHA256.json` — empreintes de contrôle.

## Limites scientifiques

La conclusion de non-confinement concerne le modèle linéaire réduit étudié.

Elle ne constitue pas, à elle seule, une démonstration concernant :
- le cube physique complet ;
- les effets non linéaires ;
- les phénomènes d'amortissement ;
- toutes les conditions aux limites possibles ;
- l'existence ou l'absence de solitons et de breather.

## Méthode de validation

Les statuts seront distingués explicitement :

- **Démontré** : preuve mathématique complète.
- **Reproduit** : calcul exécuté indépendamment avec résultat concordant.
- **Rapporté** : résultat présent dans un artefact antérieur.
- **Conditionnel** : résultat dépendant d'hypothèses ou de sources non encore vérifiées.
- **Ouvert** : contrôle restant à effectuer.

Les versions antérieures doivent être conservées. Aucun livrable historique ne doit être remplacé silencieusement.

---

Projet de recherche : CUBE-T2G  
Version documentaire initiale.
