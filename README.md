# ST2 (Théorie des Jeux) - EI Algorithmique Génétique

![CentraleSupelec Logo](https://www.centralesupelec.fr/sites/all/themes/cs_theme/medias/common/images/intro/logo_nouveau.jpg)

## MEMBRES DU GROUPE

Alexis VRIELYNCK
Tonino DIONISI
William AUROUX
Erwan BERNE
Timothée BOURRET

## CONTEXTE
Les nucléotides caractérisent ainsi l’information génétique par leur disposition sur un brin d’ADN. Au sein
d’une cellule, l’ADN permet alors de gérer les interactions de la cellule avec son environnement, où le positionnement dans l’espace du brin est important. Cela motive alors l’étude de la trajectoire tri-dimensionnelle d’une
molécule d’ADN. En 1993, des biophysiciens source à citer ont établi un modèle de conformation 3D permettant
de transformer une séquence ADN en une trajectoire tri-dimensionnelle.
On observe régulièrement dans la nature que les brins d’ADN sont circulaires, c’est-à-dire que les deux
extrémités d’un brin sont reliées. On souhaite alors pouvoir trouver les paramètres adaptés sur le modèle 3D
permettant d’obtenir une trajectoire tri-dimensionnelle circulaire, autrement dit où les nucléotides des extrémités
sont en position exacte pour qu’elles soient juxtaposées.
Ainsi, la problématique est la suivante : une table de rotation associe à une paire de nucléotide juxtaposées un triplet d’angles décrivant une rotation en coordonnées sphériques. Étant donné une séquence ADN, comment choisir une table de rotation permettant d’obtenir un brin d’ADN associé circulaire ? 

## DESCRIPTION
Ce programme propose une solution au problème considéré via deux approches différentes : une utilisant le principe du recuit simulé et l’autre consistant en un algorithme génétique.


## RESSOURCES
Sont fournis :

- Le fichier <tt>Traj3D.py</tt> implémentant le moteur de calcul d’une trajectoire 3D
- Le fichier <tt>Rot_Table.py</tt> contenant la table de rotation (avec les écarts-types) nécessaires au calcul d’une trajectoire 3D. Ce fichier contient également l’ensemble des outils permettant de traiter ces tables de rotations et les populations (ensemble de tables de rotations)
- Le fichier <tt>Main.py</tt> permettant d’exécuter les algorithmes recuit simulé et algorithme génétique
- Le fichier <tt>recuit_simulé.py</tt> implémentant l’algorithme du recuit simulé
- Le fichier <tt>genetique.py</tt> implémentant l’algorithme génétique
- L’ensemble de fichiers tests <tt>tests</tt>
Pour lancer l’ensemble des tests et effectuer un coverage, exécuter : <code>coverage run -m pytest tests</code>. 
Pour obtenir le bilan du coverage, exécuter <code>coverage html</code> puis ouvrir le fichier <tt>index.html</tt> dans le dossier htmlcov.


## EXECUTION
- Exécution : <code>python -m dna --help</code>
- L'exécution du paramètre <code>--pace</code> ou <code>-p</code> permet une exécution plus rapide et utilisant moins les capacités du processeur. Sans ce paramètre, le programme produit un meilleur résultat mais est plus lent et utilise toutes les capacités du processeur. Il n'est pas recommandé d'utiliser le mode lent pour le plasmide 180k (temps d'exécution long)
- Pour obtenir un résultat, utiliser les paramètres <code>-r</code> ou <code>-g</code> pour utiliser le recuit simulé ou l'algorithme génétique.

Après exécution, le programme affiche la table de rotation trouvée et la trajectoire 3D de la séquence d'ADN.