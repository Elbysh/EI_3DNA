from dna.RotTable import RotTable, import_RotTable, population

def test_croisement() :
    pop = population(taille=100)
    pop.creation_liste()
    pop.croisement()
    assert pop.taille == 200
    for indi in pop.liste_individus :
        assert indi.is_valid()
    
