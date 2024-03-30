from dna.RotTable import RotTable, import_RotTable
from dna.recuit_simulé import recuit_simulé
from random import seed


def test_recuit_simulé():
    """ Test de recuit simulé : on compare l'énergie de la table de rotation obtenue par recuit simulé
        avec celle de la table de rotation initiale (qui n'est pas très bonne) 
    """
    # assert les écarts types de bases avec la rot_table modifiée
    lineList_8k = [line.rstrip('\n') for line in open('data/plasmid_8k.fasta')]
    seq8k = ''.join(lineList_8k[1:])
    seed(1)
    rot_table_object, energie_list, k_list, traj = recuit_simulé(
        import_RotTable(), seq8k, 1000, 25, 0.1)
    rot_table = rot_table_object.rot_table
    seed(1)
    assert recuit_simulé(import_RotTable(), seq8k, 1000, 25, 0.1)[
        0].rot_table == rot_table and rot_table_object.is_valid()


# test_recuit_simulé()
