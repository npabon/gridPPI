Synopsis

GridPPI is a Python module for comparative scoring of predicted protein docking poses. Docked poses are scored via mapping of the protein ligand interface residues to the receptor's probe enrichment grids. These enrichment grids are generated from a druggability simulation of the unbound receptor with probe molecules. Poses that place receptor residues in within the enriched regions of their correspondiong probe grids receive higher scores.


USAGE

The 'score.py' module defines the base Score class that stores probe enrichment grids for a particular receptor and scores predicted poses of a ligand docking to that receptor. Score objects are instantiated with the receptor's probe grids that are to be used in scoring.

Ex: :
  
  myScore = Score(polar='receptor_polar.dx', hydrophobic='receptor_hydrophobic.dx')


The 'scorePose' method is then used to score a ProDy atom selection representing the ligand in the docked pose. The ligand must be aligned to the receptor grids prior to scoring. The ligand atom selection must also only contain atoms that are within the bounds of the receptor enrichment grids.

Ex::
  
  r = parsePDB('receptor_heavyatoms.pdb')
  p = parsePDB('dock_1.pdb')
  p = matchAlign(p, r)[0]
  m = matchChains(p, r, subset='all')
  pm = m[0][0]
  l = ~pm
  l = l.select('within 4 of thing', thing=r)
  testScore = myScore.scorePose(l)
  


