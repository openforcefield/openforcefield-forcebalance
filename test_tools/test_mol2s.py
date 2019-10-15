import os

from openforcefield.topology import Molecule as Off_Molecule
from openforcefield.topology import Topology as Off_Topology
from openforcefield.typing.engines.smirnoff import ForceField

test_ff = ForceField("../../forcefield/param_valence.offxml", allow_cosmetic_attributes=True)


for f in os.listdir('.'):
    if f.endswith('mol2'):
        print(f)
        off_molecule = Off_Molecule.from_file(f)
        off_topology = Off_Topology.from_molecules(off_molecule)
        test_ff.create_openmm_system(off_topology)