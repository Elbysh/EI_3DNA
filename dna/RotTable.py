from json import load as json_load
from os import path as os_path
from mathutils import (
    Matrix,
    Vector,
)
import math
from random import uniform, choice, randint, random, shuffle, randrange

from multiprocessing import Pool


here = os_path.abspath(os_path.dirname(__file__))


class RotTable:
    """Represents a rotation table"""

    # 3 first values: 3 angle values
    # 3 last values: SD values

    def __init__(self, d={"AA": [0 for i in range(6)], "AC": [0 for i in range(6)], "AG": [0 for i in range(6)], "AT": [0 for i in range(6)], "CA": [0 for i in range(6)], "CC": [0 for i in range(6)], "CG": [0 for i in range(6)], "CT": [0 for i in range(6)], "GA": [0 for i in range(6)], "GC": [0 for i in range(6)], "GG": [0 for i in range(6)], "GT": [0 for i in range(6)], "TA": [0 for i in range(6)], "TC": [0 for i in range(6)], "TG": [0 for i in range(6)], "TT": [0 for i in range(6)]}):
        self.rot_table = d

    ###################
    # WRITING METHODS #
    ###################
    def setTwist(self, dinucleotide: str, value: float):
        self.rot_table[dinucleotide][0] = value

    def setWedge(self, dinucleotide: str, value: float):
        self.rot_table[dinucleotide][1] = value

    def setDirection(self, dinucleotide: str, value: float):
        self.rot_table[dinucleotide][2] = value

    ###################
    # READING METHODS #
    ###################
    def getTwist(self, dinucleotide: str) -> float:
        return self.getTable()[dinucleotide][0]

    def getWedge(self, dinucleotide: str) -> float:
        return self.getTable()[dinucleotide][1]

    def getDirection(self, dinucleotide: str) -> float:
        return self.getTable()[dinucleotide][2]

    def getTable(self) -> dict:
        return self.rot_table

    ###################

    def compute(self, dna_seq: str):
        """ returns the 3d trajectory of a DNA sequence, given a rotation table

        Args:
            dna_seq : DNA sequences, string

        Returns:
            traj : trajectory 3d of the DNA sequence, vectors list
        """

        MATRIX_T = Matrix.Translation((0.0, 0.0, 3.38/2.0, 1.0))
        # Matrice cumulant l'ensemble des transformations géométriques engendrées par la séquence d'ADN
        total_matrix = Matrix()

        # On enregistre la position du premier nucléotide
        traj = [Vector((0.0, 0.0, 0.0, 1.0))]

        matrices_Rz = {}
        matrices_Q = {}
        # On parcourt la sequence, nucléotide par nucléotide
        for i in range(1, len(dna_seq)):
            # On lit le dinucleotide courant
            dinucleotide = dna_seq[i-1]+dna_seq[i]
            # On remplit au fur et à mesure les matrices de rotation
            if dinucleotide not in matrices_Rz:
                matrices_Rz[dinucleotide] = Matrix.Rotation(
                    math.radians(self.getTwist(dinucleotide)/2),
                    4,
                    'Z'
                )
                matrices_Q[dinucleotide] = \
                    Matrix.Rotation(
                        math.radians((self.getDirection(dinucleotide)-90)),
                        4,
                        'Z'
                ) \
                    @ Matrix.Rotation(
                        math.radians((-self.getWedge(dinucleotide))),
                        4,
                        'X'
                ) \
                    @ Matrix.Rotation(
                        math.radians((90-self.getDirection(dinucleotide))),
                        4,
                        'Z'
                )

            # On calcule les transformations géométriques
            # selon le dinucleotide courant,
            # et on les ajoute à la matrice totale
            total_matrix @= \
                MATRIX_T \
                @ matrices_Rz[dinucleotide] \
                @ matrices_Q[dinucleotide] \
                @ matrices_Rz[dinucleotide] \
                @ MATRIX_T

            # On calcule la position du nucléotide courant
            # en appliquant toutes les transformations géométriques
            # à la position du premier nucléotide
            traj.append(total_matrix @ traj[0])
        return traj

    def energy(self, dna_seq):
        """ Computes how close are the first and the last segment of the constructed DNA strand to be linked. The DNA strand is constructed from dna_seq with respect of the given RotTable

        Args:
            dna_seq (str): The DNA sequence

        Returns:
            int: The distance which describes how close are the first and the last segment to be linked correctly
        """
        traj = self.compute(dna_seq + dna_seq[0:2]) #calcule la trajectoire augmentée en ajoutant les deux premiers tout à la fin
        dernier_point, image_premier_deuxieme = traj[-2], traj[-1] #point n et n+1
        deuxieme_point = traj[1]
        premier_point = traj[0]
        vector_predict = image_premier_deuxieme - dernier_point #vecteur entre le point n et n+1
        vector_data = traj[1] - traj[0] #vecteur entre le point 1 et 2
        scalar_product = vector_predict.dot(vector_data) 
        product_norm = vector_predict.magnitude * vector_data.magnitude
        cos_angle = scalar_product / product_norm #cos de l'angle entre ces deux vecteurs
        cos_angle_norm = (1 + cos_angle)/2
        return math.sqrt(((dernier_point-premier_point).magnitude)**2 + ((deuxieme_point-image_premier_deuxieme).magnitude)**2)/cos_angle_norm

    def neighbors(self, reference_table):
        """ Returns a random RotTable in the neighborhood of the initial RotTable, without violating the standard devitions relative to the means precised in the reference RotTable

        Args:
            reference_table (RotTable): The RotTable wich describe the domain where the returned neighbor has to be in

        Returns:
            RotTable: The neighbor, it is determined by a uniform distribution on a smaller neighborhood of the initial RotTable, which depends on the fixed constant k
        """

        def completion(dict):
            """ Auxiliar function for treating constraints between values in a RotTable """

            nouv_dict = {}
            for elem in dict:
                if elem == "AA":
                    nouv_dict["TT"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                if elem == "AC":
                    nouv_dict["GT"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                if elem == "AG":
                    nouv_dict["CT"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                if elem == "CA":
                    nouv_dict["TG"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                if elem == "CC":
                    nouv_dict["GG"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                if elem == "GA":
                    nouv_dict["TC"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                nouv_dict[elem] = dict[elem]
            return nouv_dict 
        
        initial = self.rot_table
        ref = reference_table.rot_table
        k = 30 #plus k est grand, plus le voisin généré sera proche

        # Nucleotide pairs to treat. The others are redundant

        list_dinucleotide = ["AA", "AC", "AG", "AT",
                             "CA", "CC", "CG", "GA", "GC", "TA"]

        next = {}

        for r in list_dinucleotide:

            # Standard deviations
            next[r] = [0, 0, 0, ref[r][3], ref[r][4], ref[r][5]] #les écarts-type ne sont pas changés

            # Twist, Wedge and Direction
            for i in range(3):
                a = max(ref[r][i]-ref[r][i+3], initial[r][i]-ref[r][i+3]/k)
                b = min(ref[r][i]+ref[r][i+3], initial[r][i]+ref[r][i+3]/k)
                next[r][i] = uniform(a, b)

        return RotTable(completion(next))

    def mutation(self, reference_table, iterencours):
        """ Returns a randomly mutated RotTable of the initial RotTable. It modifies a random parameter of the RotTable, on the same principle as the neighbors function

        Args:
            reference_table (RotTable): The RotTable wich describe the domain where the returned mutated has to be in

        Returns:
            RotTable: The mutated RotTable, on the same principle as neighbors
        """
        def completion(dict):
            """ Auxiliar function for treating constraints between values in a RotTable """

            nouv_dict = {}
            for elem in dict:
                if elem == "AA":
                    nouv_dict["TT"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                if elem == "AC":
                    nouv_dict["GT"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                if elem == "AG":
                    nouv_dict["CT"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                if elem == "CA":
                    nouv_dict["TG"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                if elem == "CC":
                    nouv_dict["GG"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                if elem == "GA":
                    nouv_dict["TC"] = [dict[elem][0], dict[elem][1], -dict[elem]
                                       [2], dict[elem][3], dict[elem][4], dict[elem][5]]
                nouv_dict[elem] = dict[elem]
            return nouv_dict
        initial = self.rot_table
        ref = reference_table.rot_table
        k = 1

        # Nucleotide pairs to treat. The others are redundant
        list_dinucleotide = ["AA", "AC", "AG", "AT",
                             "CA", "CC", "CG", "GA", "GC", "TA"]

        next = {}

        # Random choice of the parameter
        c = choice(list_dinucleotide)
        v = randint(0, 2)

        for r in list_dinucleotide:

            # Standard deviations
            next[r] = [0, 0, 0, ref[r][3], ref[r][4], ref[r][5]]

            # Twist, Wedge and Direction

            for i in range(3):
                if r == c and i == v:
                    a = max(ref[r][i]-ref[r][i+3], initial[r]
                            [i]-ref[r][i+3]/(k+iterencours/50))
                    b = min(ref[r][i]+ref[r][i+3], initial[r]
                            [i]+ref[r][i+3]/(k+iterencours/50))
                    next[r][i] = uniform(a, b)
                else:
                    next[r][i] = initial[r][i]

        return RotTable(completion(next))

    def is_valid(self):
        """ Teste si la table de rotation est valide, 
        c'est à dire les valeurs restent dans la plage donnée par les écarts-types et les conditions de 
        symétrie par rapport aux complémentaires sont respectées"""
        ref_table = import_RotTable()
        list_dinucleotide = ["AA", "AC", "AG", "AT",
                             "CA", "CC", "CG", "GA", "GC", "TA", "GT", "CT", "TG", "GG", "TC"]
        ref_table = ref_table.getTable()
        rot_table_to_test = self.getTable()
        for dinucleotide in list_dinucleotide: 
            for i in range(3):
                if rot_table_to_test[dinucleotide][i] < ref_table[dinucleotide][i] - ref_table[dinucleotide][i+3] or rot_table_to_test[dinucleotide][i] > ref_table[dinucleotide][i] + ref_table[dinucleotide][i+3]:
                    return False #vérifie les conditions sur les écarts-types
        return True and rot_table_to_test == completion_dico(rot_table_to_test) #vérifie les conditions géométriques


def import_RotTable(filename: str = 'table.json'):
    """ Allows a RotTable to be imported from a json file

    Args:
        filename (str, optional): The filename of the json file in the current directory. Defaults to "table.json".

    Returns:
        RotTable: The exctracted RotTable from the specified json file. 
    """
    return RotTable(json_load(open(os_path.join(here, filename)))) 


def completion_dico(dict):
    """ Auxiliar function for treating constraints between values in a RotTable """
    nouv_dict = {}

    for elem in dict:
        if elem == "AA":
            nouv_dict["TT"] = [dict[elem][0], dict[elem][1], -dict[elem]
                               [2], dict[elem][3], dict[elem][4], dict[elem][5]]
        if elem == "AC":
            nouv_dict["GT"] = [dict[elem][0], dict[elem][1], -dict[elem]
                               [2], dict[elem][3], dict[elem][4], dict[elem][5]]
        if elem == "AG":
            nouv_dict["CT"] = [dict[elem][0], dict[elem][1], -dict[elem]
                               [2], dict[elem][3], dict[elem][4], dict[elem][5]]
        if elem == "CA":
            nouv_dict["TG"] = [dict[elem][0], dict[elem][1], -dict[elem]
                               [2], dict[elem][3], dict[elem][4], dict[elem][5]]
        if elem == "CC":
            nouv_dict["GG"] = [dict[elem][0], dict[elem][1], -dict[elem]
                               [2], dict[elem][3], dict[elem][4], dict[elem][5]]
        if elem == "GA":
            nouv_dict["TC"] = [dict[elem][0], dict[elem][1], -dict[elem]
                               [2], dict[elem][3], dict[elem][4], dict[elem][5]]
            
        nouv_dict[elem] = dict[elem]
    return nouv_dict


def degre_lib():
    """renvoit la liste des 10 di-nucléotides qui n'ont pas de liens géométriques entre eux"""
    return ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "GA", "GC", "TA"]


class population:

    """Une population regroupe un ensemble d'individus ROtTable """
    def __init__(self, taille=100, indiv_ref=import_RotTable()):

        self.liste_individus = [] #list de RotTable
        self.taille = taille #taille de la population
        self.individu_ref = indiv_ref #RotTable de référence

    def creation_liste(self):

        """Initialise les individus de la population de manière aléatoire, tout en respectant les conditions de validité"""
        dico_indiv_initial = self.individu_ref.getTable()
        liste_dinucleotides = degre_lib()
        ecarts_types = {
            dinucleo: dico_indiv_initial[dinucleo][3:] for dinucleo in liste_dinucleotides}

        for _ in range(self.taille):
            distrib = {}
            for nucl in liste_dinucleotides:
                tirage = [2*random()*ecarts_types[nucl][i]+dico_indiv_initial[nucl]
                          [i]-ecarts_types[nucl][i] for i in range(3)]
                distrib[nucl] = tirage

            indiv = RotTable(completion_dico(
                {dinucleo: distrib[dinucleo]+ecarts_types[dinucleo] for dinucleo in liste_dinucleotides}))
            self.liste_individus.append(indiv)

    def tournoi(self, dna_seq: str, nombre_coeurs):
        
        """Modifie la liste des individus de la population pour ne garder que la meilleur moitié
        Cette opération divise donc par deux la taille de la population"""
        
        n = self.taille
        with Pool(nombre_coeurs) as p: #pour paralléliser les calculs d'energy
            results = p.starmap(energieparallele, [
                                (self.liste_individus[i], dna_seq, i) for i in range(n)])


        Best = min(results)
        results.remove(Best) #le meilleur ne participe pas au tournoi 
        list_indices = [Best[1]] #liste des indices des gagnants
        if n % 2 == 0:
            Worst = max(results)
            results.remove(Worst) #si il reste un nombre impair après avoir retiré le meilleur, ou supprime le dernier
        shuffle(results)
        n = len(results)
        k = 500
        for i in range(n//2):

            if 0 == randint(0, k-1): #probabilité de 1/k que ce soit le moins bon du duel qui est conservé
                if results[i] >= results[i+n//2]:
                    list_indices.append(results[i][1])
                else:
                    list_indices.append(results[i+n//2][1])
            else:
                if results[i] < results[i+n//2]: #test pour savoir quelle ROtTable est meilleur
                    list_indices.append(results[i][1]) #ajout du gagnant du duel à la liste des indices des gagnants
                else:
                    list_indices.append(results[i+n//2][1]) ##ajout du gagnant du duel à la liste des indices des gagnants

        self.taille = len(list_indices)
        self.liste_individus = [self.liste_individus[i] for i in list_indices] 

    def croisement(self):
        
        """Modifie la liste des individus d'une population : chaque paire d'individu donne naissance à deux nouveaux individus
        Cette opération double donc la taille de la population"""
        
        lb = degre_lib()
        shuffle(lb)
        shuffle(self.liste_individus) #mélange les individus pour créer des paires aléatoires
        individus_croisés = self.liste_individus.copy()

        nbre_indi = self.taille

        for k in range(nbre_indi//2):

            # croisement entre l'individu 2k et 2k+1 :

            ind_1 = self.liste_individus[2*k].rot_table #parent 1
            ind_2 = self.liste_individus[2*k+1].rot_table #parent 2

            ind_croisé_1 = {} #fils 1
            ind_croisé_2 = {} #fils 2

            p = [random() for _ in range(10)]

            for i in range(10):
                #la ligne i du fils 1 est la somme de p*parent_1 et (1-p)*parent_2
                ind_croisé_1[lb[i]] = [(p[i]*ind_1[lb[i]][k]+(1-p[i])*ind_2[lb[i]][k])
                                       for k in range(2)] + [ind_1[lb[i]][k] for k in range(2, 6)]
                #la ligne i du fils 2 est la somme de (p-1)*parent_1 et p*parent_2
                ind_croisé_2[lb[i]] = [(p[i]*ind_2[lb[i]][k]+(1-p[i])*ind_1[lb[i]][k])
                                       for k in range(2)] + [ind_1[lb[i]][k] for k in range(2, 6)]

            individus_croisés.append(RotTable(completion_dico(ind_croisé_1)))
            individus_croisés.append(RotTable(completion_dico(ind_croisé_2)))
        self.taille = len(individus_croisés)
        self.liste_individus = individus_croisés


def energieparallele(elem, dna_seq, numero):
    """Même fonction qu'energy, mais est appelée de manière parallélisée"""

    traj = RotTable.compute(elem, dna_seq + dna_seq[0:2])
    dernier_point, image_premier_deuxieme = traj[-2], traj[-1] #point n et n+1
    deuxieme_point = traj[1]
    premier_point = traj[0]
    
    vector_predict = image_premier_deuxieme - dernier_point #vecteur entre le point n+1 et n+2
    vector_data = traj[1] - traj[0] #vecteur entre le point 1 et 2
    scalar_product = vector_predict.dot(vector_data)
    product_norm = vector_predict.magnitude * vector_data.magnitude
    cos_angle = scalar_product / product_norm #cos de l'angle entre ces deux vecteurs
    cos_angle_normed = (1 + cos_angle)/2

    return (math.sqrt(((dernier_point - premier_point).magnitude)**2 + ((deuxieme_point-image_premier_deuxieme).magnitude)**2)/cos_angle_normed, numero)
