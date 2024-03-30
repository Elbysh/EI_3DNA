from dna.RotTable import RotTable, import_RotTable
from dna.Traj3D import Traj3D
import numpy as np


def test_compute():
    r = import_RotTable()
    lineList_8k = [line.rstrip('\n') for line in open('data/plasmid_8k.fasta')]
    seq8k = ''.join(lineList_8k[1:])
    traj8k = Traj3D()

    lineList_180k = [line.rstrip('\n')
                     for line in open('data/plasmid_180k.fasta')]
    seq180k = ''.join(lineList_180k[1:])
    traj180k = Traj3D()

    traj8k.compute(seq8k, r)
    traj180k.compute(seq180k, r)

    assert (r.compute(seq8k) == traj8k.getTraj()
            and r.compute(seq180k) == traj180k.getTraj())


# test_compute()
