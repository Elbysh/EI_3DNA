from dna.RotTable import RotTable, import_RotTable


def test_getTwist():
    rot_table = import_RotTable()
    assert rot_table.getTwist("AA") == 35.62


def test_getWedge():
    rot_table = import_RotTable()
    assert rot_table.getWedge("AA") == 7.2


def test_getDirection():
    rot_table = import_RotTable()
    assert rot_table.getDirection("AA") == -154


def test_is_valid():
    rot_Table = import_RotTable()
    rot_Table_bis = import_RotTable()
    rot_Table_bis.rot_table["AA"][1] = 0.61  # juste en dehors de l'écart-type
    assert rot_Table.is_valid() == True and rot_Table_bis.is_valid() == False
