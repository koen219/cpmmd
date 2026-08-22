from abc import ABC, abstractmethod
from typing import Optional, List, Tuple, Dict, Set
from .network import Network
from .network import BOND, BONDTYPE, BEADID, FIBREID, BONDID, Network

import logging

_logger = logging.getLogger(__name__)

class CrosslinkDistributer(ABC):
    def distribute_crosslinkers(self, network: Network):
        bonds_to_add = list()
        types_to_add = list()

        beads_with_crosslinker = set()

        selected_bonds_and_types = self.select_bonds(network)

        for bond, bond_typ in selected_bonds_and_types:
            sorted_bond: List[int] = sorted(bond)

            if (
                sorted_bond in bonds_to_add
                or sorted_bond[0] in beads_with_crosslinker
                or sorted_bond[1] in beads_with_crosslinker
            ):
                continue

            bonds_to_add.append(bond)
            types_to_add.append(bond_typ)
            beads_with_crosslinker.add(bond[0])
            beads_with_crosslinker.add(bond[1])

        network.bonds_groups.extend(bonds_to_add)
        network.bonds_types.extend(types_to_add)

    @abstractmethod
    def select_bonds(self, network: Network) -> List[Tuple[BOND, BONDTYPE]]:
        pass

class TipToTailCrosslinkDistributer(CrosslinkDistributer):
    def __init__(self, number_of_beads_per_strand, number_of_strands):
        self._num_beads_per_strand = number_of_beads_per_strand
        self._num_strands = number_of_strands
        
    def select_bonds(self, network: Network): 
        selected_bonds = []
        crosslink_type = "crosslinker"
        
        network.details_of_bondtypes[crosslink_type] = {'r0': 0, 'k': 1}

        for i in range(self._num_strands-1):
            b0 = self._num_beads_per_strand * i + self._num_beads_per_strand - 1
            b1 = b0 + 1
            selected_bonds.append( ((b0, b1), crosslink_type))
            
        return selected_bonds

