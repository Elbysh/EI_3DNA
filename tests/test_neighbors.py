from dna.RotTable import RotTable
from dna.RotTable import import_RotTable

ref = import_RotTable()
list_dinucleotide = ["AA", "AC", "AG", "AT",
                     "CA", "CC", "CG", "GA", "GC", "TA"]


def test_neighbors():
    k = 30
    for j in range(30):
        t = ref.neighbors(ref)
        u = t.neighbors(ref)
        for dinucleotide in list_dinucleotide:
            for i in range(3):
                assert abs(u.rot_table[dinucleotide][i] - t.rot_table[dinucleotide]
                           [i]) <= ref.rot_table[dinucleotide][i+3]/k
                assert abs(u.rot_table[dinucleotide][i] - ref.rot_table[dinucleotide]
                           [i]) <= ref.rot_table[dinucleotide][i+3]
