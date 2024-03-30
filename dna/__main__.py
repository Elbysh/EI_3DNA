from .RotTable import RotTable, import_RotTable
from .Traj3D import Traj3D
from .recuit_simulé import recuit_simulé
from .genetique import algorithme_génétique
import time
import argparse
import matplotlib.pyplot as plt
import time

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="input filename of DNA sequence")
# Ajout d'argument pour lancer les algorithmes de recuit ou génétique
recuit_bool = parser.add_argument(
    '-r', '--recuit', action='store_true', help="launch the recuit algorithm")
ag_bool = parser.add_argument(
    '-g', '--genetique', action='store_true', help="launch the genetic algorithm")
# Ajout d'argument pour lancer les algorithmes de recuit ou génétique de manière rapide ou optimale (par défaut)
pace_fast = parser.add_argument('-rapide', '--pace',
                                action='store_true', help="entrer -rapide pour avoir un résultat rapide")
parser.parse_args()
args = parser.parse_args()


def main():

    # Read file
    lineList = [line.rstrip('\n') for line in open(args.filename)]
    # Formatting
    seq = ''.join(lineList[1:])

    if args.pace:
        if args.recuit:
            rot_table = import_RotTable()
            rot_table, energy_list, k_list, traj_fin = recuit_simulé(
                rot_table, seq, 1000, 300, 0.01)
            traj = Traj3D()
            traj.compute(seq, rot_table)

            print(rot_table.getTable())  # Affiche la table de rotation
            print(rot_table.energy(seq))  # Affiche l'énergie de la table de rotation

            traj.draw()
            traj.write(args.filename+".png")

        if args.genetique:
            rot_table = algorithme_génétique(500, seq, 2, 50)
            traj = Traj3D()

            print(rot_table.getTable())  # Affiche la table de rotation
            print(rot_table.energy(seq))  # Affiche l'énergie de la table de rotation

            traj.compute(seq, rot_table)
            traj.draw()
            traj.write(args.filename+".png")
    else:
        if args.recuit:
            rot_table = import_RotTable()
            rot_table, energy_list, k_list, traj_fin = recuit_simulé(
                rot_table, seq, 50000, 300, 0.01)
            traj = Traj3D()
            traj.compute(seq, rot_table)

            print(rot_table.getTable())
            print(rot_table.energy(seq))

            traj.draw()
            traj.write(args.filename+".png")

        if args.genetique:
            rot_table = algorithme_génétique(10000, seq, None, 50)
            traj = Traj3D()

            print(rot_table.getTable())
            print(rot_table.energy(seq))

            traj.compute(seq, rot_table)
            traj.draw()
            traj.write(args.filename+".png")


def benchmark():
    # Initialize seq and other necessary parameters
    lineList = [line.rstrip('\n') for line in open(args.filename)]
    seq = ''.join(lineList[1:])
    rot_table = import_RotTable()
    itérations_recuit = [10, 100, 1000, 10000, 50000]
    itérations_AG = [10, 50, 100, 150, 200]
    taille_pop = [10, 100, 500, 1000, 5000]
    temperature = 200
    precision = 0.01

    # Benchmark for recuit_simulé
    temps_exec_recuit = {i: 0 for i in itérations_recuit}
    res_exec_recuit = {i: 0 for i in itérations_recuit}
    for i in itérations_recuit:
        start_time = time.time() # On lance le chrono
        rot_table, energy_list, k_list, traj = recuit_simulé(
            rot_table, seq, i, temperature, precision)
        end_time = time.time() # On arrête le chrono

        temps_exec_recuit[i] = end_time-start_time
        res_exec_recuit[i] = rot_table.energy(seq), rot_table.getTable()

    # Benchmark for genetique
    temps_exec_AG = {(i, j): 0 for i in itérations_AG for j in taille_pop}
    res_exec_AG = {i: 0 for i in itérations_AG}
    for iter in itérations_AG:
        for taille in taille_pop:
            start_time = time.time() # On lance le chrono
            rot_table = algorithme_génétique(taille, seq, None, iter)
            end_time = time.time() # On arrête le chrono
            temps_exec_AG[(iter, taille)] = end_time-start_time
            res_exec_AG[(iter, taille)] = rot_table.energy(
                seq), rot_table.getTable()
            
    print(temps_exec_recuit)
    print(" ")
    print(res_exec_recuit)
    print(" ")
    print(temps_exec_AG)
    print(" ")
    print(res_exec_AG)


if __name__ == "__main__":
    main()
