import numpy as np

import grid

"""This module defines a class for scoring docked PPI's."""

class Score(object):
    """Base class for handling grid files and scoring docked poses."""

    # atom selections for each grid type (refine later to select atoms)
    grid_types = {'positive' : '(resname ARG HIS LYS) and sidechain and element N', 
                  'negative' : '(resname ASP GLU) and sidechain and element O', 
                  'polar' : '(resname SER THR ASN GLN) and sidechain and (element N O)',
                  'hydrophobic' : '(resname ALA VAL ILE LEU MET PHE TYR TRP) and sidechain and (element C S N)'}

    def __init__(self, **kwargs):
        """Receptor grid files are taken as keyword arguments at instantiation. 
	    
	    Ex: myScore = Score(polar = 'defaults_IPRO.dx')

	    """
        # should check here to make sure keywords match standard grid types
        self.grids = kwargs
        
    def scorePose(self, ligand, selection=None):
        """Method for scoring a docked pose.

        Ligand atom selection should already be aligned to the receptor. It
        should also contain only those atoms that are within the receptor 
        grid.

        Ex: ligand = ligand.select('within 4 of thing', thing = receptor)

        Optional input selection string is combined with each grid's default
        selection string.

        The final score of a docking is a sum of values in the voxels that are 
        inhabited by the selected atoms on the ligand.

        """

        score = 0

        for grid_tag in self.grids:
            
            # get coordinates of appropriate atoms in ligand
            selstr = Score.grid_types[grid_tag]
            lig_atoms = ligand.select(selstr)
            if selection:
                lig_atoms = lig_atoms.select(selection)
            
            if lig_atoms:
                lig_coords = lig_atoms.getCoords()

                # sum values in corresponding voxels
                gridDX = grid.OpenDX(self.grids[grid_tag])
                score += sum(gridDX[lig_coords])
                print grid_tag, sum(gridDX[lig_coords])

        return score
