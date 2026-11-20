# Modèle cubique exact — construction à 59 degrés de liberté

> **Statut documentaire.** Cette page rassemble les résultats exacts obtenus et rapportés dans les audits V14–V16 du projet. Elle distingue les démonstrations algébriques déjà réalisées de leur reproductibilité depuis ce dépôt. Les sources historiques V16.1 ne sont pas encore toutes archivées ici ; leur réintégration est nécessaire pour refaire l'ensemble de la chaîne à partir de zéro.

## 1. Objet du modèle

Le modèle étudie une construction vectorielle sur le cube de référence
\[
\Omega=[-\tfrac12,\tfrac12]^3,
\]
en tenant compte de ses symétries cubiques. Les champs sont construits à partir de composantes polynomiales homogènes de degrés 1, 3 et 5. L'objectif est de construire des espaces de degrés de liberté compatibles avec les représentations irréductibles cubiques, puis d'étudier exactement les matrices de masse, de raideur et les couplages dynamiques associés.

La symétrie considérée est celle du groupe octaédrique complet \(O_h\), représenté par les 48 matrices de permutations signées des coordonnées.

## 2. Décomposition des 59 degrés de liberté

L'espace retenu se décompose en deux secteurs :

| Secteur | Dimension | Décomposition rapportée |
|---|---:|---|
| \(E_g\) | 20 | 14 degrés de liberté à divergence nulle et 6 relèvements scalaires |
| \(T_{2g}\) | 39 | 30 degrés de liberté à divergence nulle et 9 relèvements scalaires |
| **Total** | **59** | \(20+39\) |

Les bases cubiques et les matrices associées ont été construites sur les rationnels. Les résultats d'audit rapportent une invariance exacte sous les 48 opérations de \(O_h\).

Cette page décrit le modèle cubique étudié ; elle ne prétend pas que ces 59 degrés de liberté constituent à eux seuls une discrétisation universelle de tout cube élastique continu. Leur interprétation dépend de la construction variationnelle et des espaces de champs effectivement retenus.

## 3. Matrices exactes

Les audits rapportent la construction exacte de matrices rationnelles :
- \(M\), matrice de masse ;
- \(K_0\), première contribution de raideur ;
- \(K_d\), contribution de couplage ou de déformation ;
- \(C\), opérateur d'observation centrale.

Les contrôles de symétrie, de rationalité et les assemblages par secteurs ont été effectués dans les scripts de calcul historiques. Les matrices complètes sont des pièces essentielles de la preuve reproductible : elles devront être déposées avec leurs scripts générateurs, leurs conventions d'ordre des bases et leurs empreintes SHA-256.

## 4. Entrelaceurs rationnels et réduction par symétrie

Des entrelaceurs rationnels explicites ont été construits entre les espaces complets et les représentations irréductibles :

- secteur \(E_g\) : 10 entrelaceurs indépendants de taille \(20\times2\), assemblés en une matrice carrée \(S_E\) inversible ; déterminant rapporté : \(-1/8\) ;
- secteur \(T_{2g}\) : 13 entrelaceurs indépendants de taille \(39\times3\), assemblés en une matrice carrée \(S_T\) inversible ; déterminant rapporté : \(3/128\).

Ces constructions établissent une séparation exacte entre la dimension de représentation et la multiplicité.

Pour le secteur \(T_{2g}\), la réduction V16.1 rapportée est :
\[
S_T^{\mathsf T} A S_T=I_3\otimes A_{T,13},
\]
avec une matrice réduite de multiplicité de taille \(13\times13\).

Pour le secteur \(E_g\), la forme rapportée est :
\[
S_E^{\mathsf T} A S_E=H_E\otimes A_{E,10},
\qquad H_E=\operatorname{diag}(1,3),
\]
avec une matrice de multiplicité de taille \(10\times10\).

Ces identités sont des résultats algébriques exacts rapportés par les audits ; elles ne sont pas des approximations numériques.

## 5. Raccords exacts entre versions

Plusieurs raccords de versions ont également été rapportés comme réussis :

- raccord V14-H2 vers les bases V16.1, avec égalité exacte des coefficients dans le corps algébrique utilisé : 6 591 coefficients sur 6 591 pour le bloc principal et 5 616 sur 5 616 pour le bloc associé ;
- raccord V16.1 vers V14-G, avec zéro différence sur les 144 coefficients comparés dans \(\mathbb{Q}(\beta)\).

Ces contrôles vérifient que des constructions issues de versions différentes décrivent les mêmes objets algébriques dans les bases raccordées. Les journaux et scripts d'origine doivent accompagner les matrices pour rendre ces comparaisons entièrement reproductibles depuis le dépôt.

## 6. Observabilité centrale : résultats et limites

Les audits d'observabilité rapportent l'étude des mesures centrales dans les deux secteurs. Des polynômes exacts de conditions exceptionnelles ont été construits :

- \(P_E\), de degré 24, pour le secteur \(E_g\) ;
- \(P_T\), de degré 33, pour le secteur \(T_{2g}\).

Les racines réelles positives rapportées sont approximativement :

| Secteur | Racines approximatives rapportées |
|---|---|
| \(E_g\) | 0,65312369 ; 5,59153181 ; 23,18621228 |
| \(T_{2g}\) | 5,03703936 ; 11,27304326 |

Les audits rapportent des racines simples et des noyaux invisibles de dimension un par copie aux paramètres exceptionnels, soit des corangs d'état de 2 pour le système du second ordre dans les cas décrits.

Il faut distinguer :
1. le rang géométrique des opérateurs de mesure ;
2. le rang de Krylov ou le critère PBH pour la dynamique ;
3. l'observabilité de l'état complet du système du second ordre.

Ces notions ne sont pas interchangeables. Les résultats de rang des profils d'arêtes et des intégrales ne constituent pas, à eux seuls, une preuve d'observabilité dynamique complète.

Les trois mesures proposées sur une arête de référence et leur test PBH/GCD restent un point de travail distinct ; ils ne sont pas présentés ici comme un résultat fermé.

## 7. Lien avec V17.5

Le dépôt contient déjà les artefacts V17.5 consacrés à un théorème de non-confinement pour un modèle linéaire réduit à 13 degrés de liberté. Ce résultat porte sur la matrice réduite et sur les conditions initiales explicitement définies dans l'énoncé V17.5.

Il ne remplace pas la démonstration de la réduction des 59 degrés de liberté vers les secteurs de multiplicité, et ne doit pas être interprété comme une preuve directe sur le cube continu complet.

## 8. État de validation

| Élément | Statut |
|---|---|
| Construction des secteurs 20 et 39 DDL | Résultat exact rapporté dans les audits |
| Invariance cubique sous les 48 opérations de \(O_h\) | Vérification exacte rapportée |
| Entrelaceurs rationnels \(S_E\) et \(S_T\) | Constructions exactes rapportées |
| Réduction \(T_{2g}:39\to13\) | Identité exacte rapportée en V16.1 |
| Réduction \(E_g:20\to10\) | Identité exacte rapportée en V16.1 |
| Raccords V14-H2/V16.1 et V16.1/V14-G | Égalités exactes rapportées, avec nombres de coefficients indiqués |
| Analyse d'observabilité centrale | Résultats algébriques rapportés ; périmètre précisé ci-dessus |
| Reproduction indépendante de toute la chaîne depuis le dépôt | **À faire** : sources V16.1 et données de base à archiver |

« Exact rapporté » signifie que le résultat est annoncé comme obtenu par calcul symbolique exact dans les audits du projet. Cela ne signifie pas qu'il a été recalculé indépendamment pendant la rédaction de cette page.

## 9. Fichiers à archiver pour la reproduction complète

La publication reproductible devra contenir, sans écraser les versions antérieures :

1. les bases exactes \(E_g\) et \(T_{2g}\) ;
2. les matrices complètes \(M,K_0,K_d,C\) et leur convention d'indexation ;
3. les entrelaceurs \(S_E,S_T\) ;
4. les matrices réduites V16.1 et les scripts de construction ;
5. les scripts de vérification des invariances, des identités de réduction et des raccords de versions ;
6. les journaux d'exécution et un manifeste SHA-256.

Les fichiers V16.1 historiques identifiés comme manquants dans la branche principale devront être récupérés depuis leur source d'origine avant de pouvoir annoncer une reproduction intégrale.

## 10. Principe de publication

Le dépôt conserve les versions historiques. Toute nouvelle vérification doit être ajoutée comme un nouvel artefact, avec :
- les hypothèses et conventions ;
- le script exact exécuté ;
- le résultat brut ;
- les empreintes des fichiers d'entrée et de sortie ;
- les contrôles négatifs pertinents ;
- un statut explicite : démontré, reproduit, rapporté, conditionnel ou ouvert.

L'ambition est de rendre la chaîne géométrique et algébrique vérifiable par un tiers, sans confondre résultat mathématique, exécution informatique et interprétation physique.
