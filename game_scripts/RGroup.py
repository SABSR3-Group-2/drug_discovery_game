# from rdkit.Chem.rdchem import Mol
# from rdkit import Chem
# from rdkit.Chem import Draw
#
# class GameMol(Mol):
#     """Wrapper class of `rdkit.Chem.rdchem.Mol` to add additional functionality for molecules used in the Drug
#     Discovery Game."""
#
#     def __init__(self, mol):
#         super().__init__(mol)
#
#
#
# x = Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]')
# print(x.GetNumAtoms())
# y = GameMol(x)
# print(y.GetNumAtoms())
# print(y.MakePNG())