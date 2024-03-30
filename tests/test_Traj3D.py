from dna.Traj3D import Traj3D
from dna.RotTable import RotTable, import_RotTable
from matplotlib import pyplot as plt
from dna.genetique import algorithme_génétique


def test_Traj3D():
    # Read file
    lineList = [line.rstrip('\n') for line in open('data/plasmid_8k.fasta')]
    # Formatting
    seq = ''.join(lineList[1:])
    rot_table = import_RotTable()
    rot_table2 = algorithme_génétique(10, seq, None, 5)
    traj = Traj3D()
    traj.compute(seq, rot_table)
    traj2 = Traj3D()
    traj2.compute(seq, rot_table2)
    traj.draw()
    traj.draw2(traj2)
    # print(traj.getTraj())
    assert traj.fig.savefig('data/plasmid_8k.fasta' +
                            ".png") == traj.write('data/plasmid_8k.fasta'+".png")
