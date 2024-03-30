from .RotTable import RotTable, import_RotTable
from math import exp
from random import random
from .Traj3D import Traj3D


def recuit_simulé(rot_table_init, seq_dna, nb_iter, T_init, précision):
    """ Returns the result of the simulated annealing algorithm

    Args:
        rot_table_init (RotTable): Initial RotTable
        seq_dna (str): DNA sequence
        nb_iter (int): Number of iteration
        T_init (float): Temperature
        précision : Determines how precise the result has to be in comparison to the optimal RotTable

    Returns:
        (RotTable): The approximated RotTable
        (int): The best energy found 
        (list): The list of energies
        traj_point_final (Traj3D): trajectory of the last point trough the simulation
    """
    traj_fin = [] #cette liste stocke les coordonnées du dernier point à chaque itération
    ref_table = import_RotTable()
    s = rot_table_init
    T = T_init #témpérature
    decrease_rate = 0.99 #taux de décroissance pour la température
    e = s.energy(seq_dna) #calcule de l'energie initiale
    energy_list = [e] #liste des energies successives 
    liste_k = [0]
    k = 0
    while k < nb_iter and e > précision: #s'arrête après nb_iter itérations, ou si l'energie obtenue est excellente
        s_n = s.neighbors(ref_table)  #s_n est in voisin de s
        e_n = s_n.energy(seq_dna) #energie de s_n
        energy_list.append(e_n)
        liste_k.append(k)
        if e_n < e or random() < exp((e-e_n)/T): #on compare l'energie de s et de son voisin
            s = s_n #maj de l'individu gardé
            e = e_n #maj de l'energie
            traj = Traj3D
            Traj3D.compute(traj, seq_dna, s_n)
            traj_fin.append(Traj3D.getTraj(traj)[-1])
        T = decrease_rate * T #diminituion de la température
        k += 1 
    traj_point_final = Traj3D
    Traj3D.set_traj(traj_point_final, traj_fin)
    return s, energy_list, liste_k, traj_point_final 
