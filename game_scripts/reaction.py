'''
Script that takes a molecular scaffold and an R group, both with placeholders, and returns a molecule with the two correctly attached
'''

import rdkit 
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

def join_fragments(r_group, scaffold, r_group_vector, scaffold_vector):

    print(scaffold_vector)
    print(r_group)

    #Define the reaction
    rxn = AllChem.ReactionFromSmarts("[Rb][*:2].[Rb][*:3]>>[*:2][*:3]")
    #Run the reaction with the provided information
    product = rxn.RunReactants( [scaffold, r_group] )

    return product

#test = join_fragments(Chem.MolFromSmiles('Oc1ccc(C[*:1])cc1'), Chem.MolFromSmiles('O=C(O)C(NS(=O)(=O)c1ccc([*:2])cc1)[*:1]'),
                       # Chem.MolFromSmiles('[*:1]'), Chem.MolFromSmiles('[*:1]'))

#print(Chem.MolToSmiles(test))

test = join_fragments(Chem.MolFromSmiles("[Rb]c1ccccc1"), Chem.MolFromSmiles("[Rb]OC[Cs]"), 'a', 'b')

print(Chem.MolToSmiles(test[0][0]))

smi = 'c1cc[*:1]ccc1.0[*:1]C'
mol = Chem.MolFromSmiles(smi)
mol