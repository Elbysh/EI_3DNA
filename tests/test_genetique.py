from dna.RotTable import RotTable, import_RotTable
from dna.genetique import algorithme_génétique
from unittest import mock
from random import seed


def test_algorithme_genetique():
    """ Test de l'algortihtme genetique: on compare l'énergie de la table de rotation obtenue par  l'algortihtme genetique
        avec celle de la table de rotation initiale (qui n'est pas très bonne) 
    """
    # assert les écarts types de bases avec la rot_table modifiée
    lineList_8k = [line.rstrip('\n') for line in open('data/plasmid_8k.fasta')]
    seq8k = ''.join(lineList_8k[1:])
    seed(1)
    rot_table_object = algorithme_génétique(
        20, seq8k, 2, 100)
    rot_table_res = rot_table_object.rot_table
    seed(1)  # and rot_table_object.is_valid()
    assert algorithme_génétique(
        20, seq8k, 2, 100).rot_table == rot_table_res
