# For computing
from mathutils import (
    Matrix,
    Vector
)
import math

# For drawing
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from .RotTable import RotTable


class Traj3D:
    """Represents a 3D trajectory"""

    # Vertical translation (elevation) between two di-nucleotides
    __MATRIX_T = Matrix.Translation((0.0, 0.0, 3.38/2.0, 1.0))

    def __init__(self):
        self.__Traj3D = {}
        self.fig = plt.figure()
        self.ax = plt.axes(projection='3d')

    def getTraj(self) -> dict:
        return self.__Traj3D

    def compute(self, dna_seq: str, rot_table: RotTable):

        # Matrice cumulant l'ensemble des transformations géométriques engendrées par la séquence d'ADN
        total_matrix = Matrix()

        # On enregistre la position du premier nucléotide
        self.__Traj3D = [Vector((0.0, 0.0, 0.0, 1.0))]

        matrices_Rz = {}
        matrices_Q = {}
        # On parcourt la sequence, nucléotide par nucléotide
        for i in range(1, len(dna_seq)):
            # On lit le dinucleotide courant
            dinucleotide = dna_seq[i-1]+dna_seq[i]
            # On remplit au fur et à mesure les matrices de rotation
            if dinucleotide not in matrices_Rz:
                matrices_Rz[dinucleotide] = Matrix.Rotation(
                    math.radians(rot_table.getTwist(dinucleotide)/2),
                    4,
                    'Z'
                )
                matrices_Q[dinucleotide] = \
                    Matrix.Rotation(
                        math.radians(
                            (rot_table.getDirection(dinucleotide)-90)),
                        4,
                        'Z'
                ) \
                    @ Matrix.Rotation(
                        math.radians((-rot_table.getWedge(dinucleotide))),
                        4,
                        'X'
                ) \
                    @ Matrix.Rotation(
                        math.radians(
                            (90-rot_table.getDirection(dinucleotide))),
                        4,
                        'Z'
                )

            # On calcule les transformations géométriques
            # selon le dinucleotide courant,
            # et on les ajoute à la matrice totale
            total_matrix @= \
                self.__MATRIX_T \
                @ matrices_Rz[dinucleotide] \
                @ matrices_Q[dinucleotide] \
                @ matrices_Rz[dinucleotide] \
                @ self.__MATRIX_T

            # On calcule la position du nucléotide courant
            # en appliquant toutes les transformations géométriques
            # à la position du premier nucléotide
            self.__Traj3D.append(total_matrix @ self.__Traj3D[0])

    def draw(self):
        xyz = np.array(self.__Traj3D)
        x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
        self.ax.plot(x, y, z)
        plt.show()

    def write(self, filename: str):
        self.fig.savefig(filename)

    def set_traj(self, traj):
        self.__Traj3D = traj

    def draw2(self, self2):
        """Cette fonction permet de dessiner 2 trajectoires 3D sur le même graphique,
        pour afficher à la fois la forme finale du brin d'ADN et la trajectoire du 
        point final de la chaîne au fil des itérations

        Args:
            self2 (_type_): _description_
        """
        # On récupère les trajectoires 2 trajectoires et on les convertit en tableau numpy
        xyz = np.array(self.__Traj3D)
        xyz2 = np.array(Traj3D.getTraj(self2))

        # On récupère les coordonnées x, y et z des deux trajectoires
        x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
        x2, y2, z2 = xyz2[:, 0], xyz2[:, 1], xyz2[:, 2]

        # On dessine les deux trajectoires sur le même graphique
        self.ax.plot(x, y, z)
        self.ax.plot(x2, y2, z2)
        plt.show()
