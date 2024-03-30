from .RotTable import RotTable, population, import_RotTable, degre_lib, completion_dico
from random import random


def algorithme_génétique(nb_indiv, seq_dna, nombre_coeurs, nb_iter=100):
    """ Returns the result of the genetic algorithm

    Args:
        nb_indiv (int): Number of individuals  
        seq_dna (str): DNA sequence
        nombre_coeurs (int): Number of cores
        nb_iter (int): Number of iterations
    Returns:
        (RotTable): The approximated RotTable
    """
    pop = population(nb_indiv)
    pop.creation_liste()  # initialisation de la population
    for iterencours in range(nb_iter):
        pop.tournoi(seq_dna, nombre_coeurs)  # selection des individus
        new_liste = []
        for indiv in pop.liste_individus:  # mutation des individus
            if random() < 0.005:  # 0.5% de chance de muter
                indiv = indiv.mutation(pop.individu_ref, iterencours)
            new_liste.append(indiv)
        pop.liste_individus = new_liste
        pop.croisement()  # croisement des individus

    score = [(pop.liste_individus[i].energy(seq_dna), i)
             # calcul de l'énergie de chaque individu
             for i in range(len(pop.liste_individus))]
    Best = min(score)  # sélection le meilleur individu
    return pop.liste_individus[Best[1]]
