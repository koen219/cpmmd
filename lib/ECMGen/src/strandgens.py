from .network import Network, BEADTYPE, BONDTYPE
from .stranddistributions import StrandDistribution

import numpy as np
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple, List, Optional

import logging

_logger = logging.getLogger(__name__)


class StrandGenerator(ABC):
    @abstractmethod
    def build_strands(self, network: Network):
        pass

    @abstractmethod
    def fix_boundaries(self, network: Network):
        pass


from .parameters import RandomStrandGeneratorParameters


class RandomStrandGenerator(StrandGenerator):
    def __init__(
        self,
        par: RandomStrandGeneratorParameters,
        strand_distribution: StrandDistribution,
    ):
        self._strand_distribution = strand_distribution
        self._par = par

    def build_strands(self, network: Network):
        particlepos, types = self._pos_gen()
        bondsgroup, bondstypes = self._bond_gen()
        anglegroup, angletypes = self._angle_gen()

        network.beads_positions.extend(particlepos)
        network.beads_types.extend(types)

        network.bonds_groups.extend(bondsgroup)
        network.bonds_types.extend(bondstypes)

        network.angle_groups.extend(anglegroup)
        network.angle_types.extend(angletypes)

        return network

    def fix_boundaries(self, network: Network):
        pos = np.array(network.beads_positions)
        typeid = np.array(network.beads_types, dtype=object)

        sizex = network.domain.sizex
        sizey = network.domain.sizey
        sizez = network.domain.sizez

        if network.domain.fix_boundary or network.domain.fix_boundary_north:
            boundary_particles = pos[:, 1] > sizey
            typeid[boundary_particles] = "boundary"

        if network.domain.fix_boundary or network.domain.fix_boundary_south:
            boundary_particles = pos[:, 1] < 0
            typeid[boundary_particles] = "boundary"

        if network.domain.fix_boundary or network.domain.fix_boundary_east:
            boundary_particles = abs(pos[:, 0]) > sizex
            typeid[boundary_particles] = "boundary"

        if network.domain.fix_boundary or network.domain.fix_boundary_west:
            boundary_particles = pos[:, 0] < 0
            typeid[boundary_particles] = "boundary"

        if network.domain.fix_boundary:
            boundary_particles = abs(pos[:, 2]) > sizey
            typeid[boundary_particles] = "boundary"

        if network.domain.fix_boundary:
            boundary_particles = pos[:, 2] < 0
            typeid[boundary_particles] = "boundary"

        _logger.debug("Fixed %s boundary particles" % np.sum(typeid == "boundary"))

        network.beads_types = typeid.tolist()

    def _pos_gen(self):
        num_particles = (
            self._par.number_of_strands * self._par.number_of_beads_per_strand
        )
        num_strands = self._par.number_of_strands
        num_beads = self._par.number_of_beads_per_strand
        contour_length = self._par.contour_length_of_strand

        dist = self._strand_distribution

        middle_of_strand = num_beads // 2

        pos = np.zeros((num_strands, num_beads, 3))
        pos[:, middle_of_strand, 0] = dist.pos_x_dist(num_strands)
        pos[:, middle_of_strand, 1] = dist.pos_y_dist(num_strands)
        pos[:, middle_of_strand, 2] = dist.pos_z_dist(num_strands)

        h = contour_length / (num_beads - 1)

        angles = dist.angle_dist(num_strands)
        phis = dist.angle_dist(num_strands)
        for bead in range(num_beads):
            if bead == middle_of_strand:
                continue
            xs = np.cos(angles) * np.sin(phis)
            ys = np.sin(angles) * np.sin(phis)
            zs = np.cos(phis)
            v = h * np.column_stack([xs, ys, zs])
            pos[:, bead, :] = pos[:, middle_of_strand, :] + v * (
                middle_of_strand - bead
            )
        pos = pos.reshape((num_particles, 3))
        typeid = np.array(["free"] * num_particles, dtype=object)

        typeid = typeid.tolist()

        return pos, typeid

    def _bond_gen(self):
        num_strands = self._par.number_of_strands
        num_beads = self._par.number_of_beads_per_strand
        bondsgroup = np.empty(shape=(num_strands * (num_beads - 1), 2), dtype=int)
        bonds_indices = np.repeat(
            num_beads * np.arange(0, num_strands, dtype=int), num_beads - 1
        )
        bondsgroup[:, 0] = bonds_indices + np.tile(
            np.arange(0, num_beads - 1, 1, dtype=int), num_strands
        )
        bondsgroup[:, 1] = bonds_indices + np.tile(
            np.arange(1, num_beads, 1, dtype=int), num_strands
        )
        bondsgroup = bondsgroup.tolist()
        bondstype = ["polymer"] * len(bondsgroup)

        return bondsgroup, bondstype

    def _angle_gen(self):
        num_strands = self._par.number_of_strands
        num_beads = self._par.number_of_beads_per_strand
        angles_group = np.empty(shape=(num_strands * (num_beads - 2), 3), dtype=int)
        angles_indices = np.repeat(
            num_beads * np.arange(0, num_strands, dtype=int), num_beads - 2
        )
        angles_group[:, 0] = angles_indices + np.tile(
            np.arange(0, num_beads - 2, 1, dtype=int), num_strands
        )
        angles_group[:, 1] = angles_indices + np.tile(
            np.arange(1, num_beads - 1, 1, dtype=int), num_strands
        )
        angles_group[:, 2] = angles_indices + np.tile(
            np.arange(2, num_beads, 1, dtype=int), num_strands
        )
        anglegroup = angles_group.tolist()
        types = ["polymer_bend"] * len(anglegroup)
        return anglegroup, types
