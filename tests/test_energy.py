from dna.RotTable import RotTable, import_RotTable,energieparallele
from dna.Traj3D import Traj3D
import math


def test_energy():
    lineList = [line.rstrip('\n') for line in open("data/plasmid_8k.fasta")]
    seq = ''.join(lineList[1:])
    r = import_RotTable()
    traj = Traj3D
    Traj3D.compute(traj, seq + seq[0:2], r)
    traj = Traj3D.getTraj(traj)
    dernier_point, image_premier_deuxieme = traj[-2], traj[-1]
    vector_predict = image_premier_deuxieme - dernier_point
    vector_data = traj[1] - traj[0]
    scalar_product = vector_predict.dot(vector_data)
    product_norm = vector_predict.magnitude * vector_data.magnitude
    cos_angle = scalar_product / product_norm
    cos_angle_normed = (1 + cos_angle)/2
    assert r.energy(seq) == math.sqrt(
        (traj[-2]-traj[0]).magnitude**2 + (traj[-1]-traj[1]).magnitude**2)/cos_angle_normed


def test_energieparallel():
    lineList = [line.rstrip('\n') for line in open("data/plasmid_8k.fasta")]
    seq = ''.join(lineList[1:])
    r = import_RotTable()
    traj = Traj3D
    Traj3D.compute(traj, seq + seq[0:2], r)
    traj = Traj3D.getTraj(traj)
    dernier_point, image_premier_deuxieme = traj[-2], traj[-1]
    vector_predict = image_premier_deuxieme - dernier_point
    vector_data = traj[1] - traj[0]
    scalar_product = vector_predict.dot(vector_data)
    product_norm = vector_predict.magnitude * vector_data.magnitude
    cos_angle = scalar_product / product_norm
    cos_angle_normed = (1 + cos_angle)/2
    assert energieparallele(r, seq, 2) == (math.sqrt(
        (traj[-2]-traj[0]).magnitude**2 + (traj[-1]-traj[1]).magnitude**2)/cos_angle_normed, 2)
