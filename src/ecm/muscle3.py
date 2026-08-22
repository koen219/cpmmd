"""Helper functions for connecting to MUSCLE3

This module contains some functions for getting information from
MUSCLE3 into our custom data types, and vice versa for taking our
data and sending it to other components via MUSCLE3.
"""

from dataclasses import asdict, fields
import logging
from typing import Any, Optional, Type, TypeVar

import numpy as np

from .cell_ecm_interactions import (
    AddAdhesionParticles,
    CellECMInteractions,
    ChangeTypeInArea,
    MoveAdhesionParticles,
    RemoveAdhesionParticles,
)
from .ecm import (
    AngleCsts,
    AngleCstTypes,
    Bonds,
    BondTypes,
    MDState,
    Particles,
    ParticleType,
)
from .ecm_boundary_state import ECMBoundaryState

# from tissue_simulation_toolkit.ecm.muscle3_mpi_wrapper import Instance


_logger = logging.getLogger(__name__)


T = TypeVar("T")

# def from_settings(ParType: Type[T], instance: Instance) -> T:
#     """Create object with parameters from MUSCLE3 settings.
#
#     Note that the first argument must be a type, specifically a
#     dataclass, not an object. The result will be of this type.
#
#     Example:
#         par = from_settings(CreationParameters, instance)
#
#     Args:
#         par_type: The type of the parameter object to create
#         instance: MUSCLE3 instance to get settings from.
#     """
#     # This uses introspection, which confuses mypy, so we tell it to ignore
#     # the perceived problem.
#     return ParType(**{
#         field.name: instance.get_setting(
#             field.name, field.type.__name__)    # type: ignore [call-overload]
#         for field in fields(ParType)})          # type: ignore [arg-type]


def encode_mdstate(state: Optional[MDState]) -> Any:
    """Encode an MDState into a MUSCLE3 message-compatible data object

    This does not copy the data, so the state object passed in must not
    be changed or the result will change too.

    Args:
        state: The MDState to encode
    """
    if state is None:
        return None

    return {
        "particles": {
            "positions": state.particles.positions,
            "types": state.particles.type_ids.astype(np.int32),
        },
        "bond_types": {"r0": state.bond_types.r0, "k": state.bond_types.k},
        "bonds": {
            "groups": state.bonds.particle_groups.astype(np.int32),
            "types": state.bonds.typ.astype(np.int32),
        },
        "angle_cst_types": {
            "t0": state.angle_cst_types.t0,
            "k": state.angle_cst_types.k,
        },
        "angle_csts": {
            "groups": state.angle_csts.particle_groups.astype(np.int32),
            "types": state.angle_csts.typ.astype(np.int32),
        },
    }


def decode_mdstate(data: Any) -> MDState:
    """ "Decode an MDState from a MUSCLE3 data object.

    Args:
        data: The received message data to decode
    """
    particles = Particles(
        data["particles"]["positions"].copy(),
        data["particles"]["types"].copy(),
    )
    bond_types = BondTypes(
        data["bond_types"]["r0"].copy(), data["bond_types"]["k"].copy()
    )
    bonds = Bonds(data["bonds"]["groups"].copy(), data["bonds"]["types"].copy())
    angle_cst_types = AngleCstTypes(
        data["angle_cst_types"]["t0"].copy(),
        data["angle_cst_types"]["k"].copy(),
    )
    angle_csts = AngleCsts(
        data["angle_csts"]["groups"].copy(),
        data["angle_csts"]["types"].copy(),
    )

    return MDState(particles, bond_types, bonds, angle_cst_types, angle_csts)


def decode_cell_ecm_interactions(data: Any) -> CellECMInteractions:
    """Decode a CellECMInteractions from received data

    Args:
        data: The received data
    """
    ctia = data["change_type_in_area"]
    if ctia:
        change_type_in_area = ChangeTypeInArea(
            change_area=ctia["change_area"],
            num_particles=ctia["num_particles"],
            from_type=ParticleType(ctia["from_type"]),
            to_type=ParticleType(ctia["to_type"]),
        )
    else:
        change_type_in_area = ChangeTypeInArea()

    aap = data["add_adhesion_particles"]
    if aap:
        add_adhesion_particles = AddAdhesionParticles(
            new_pos=aap["new_pos"],
            bond_attempt_radius=aap["bond_attempt_radius"],
        )
    else:
        add_adhesion_particles = AddAdhesionParticles()

    maps = data["move_adhesion_particles"]
    if maps:
        move_adhesion_particles = MoveAdhesionParticles(
            par_id=np.array(maps["par_id"]), new_pos=np.array(maps["new_pos"])
        )
    else:
        move_adhesion_particles = MoveAdhesionParticles()

    rap = data["remove_adhesion_particles"]
    if rap:
        remove_adhesion_particles = RemoveAdhesionParticles(par_id=np.array(rap["par_id"], dtype=np.int64))
    else:
        remove_adhesion_particles = RemoveAdhesionParticles()

    return CellECMInteractions(
        change_type_in_area=change_type_in_area,
        add_adhesion_particles=add_adhesion_particles,
        move_adhesion_particles=move_adhesion_particles,
        remove_adhesion_particles=remove_adhesion_particles,
    )


def convert_arrays_to_lists(d):
    if isinstance(d, dict):
        return {k: convert_arrays_to_lists(v) for k, v in d.items()}
    elif isinstance(d, np.ndarray):
        return d.tolist()
    else:
        return d


def encode_ecm_boundary_state(boundary: Optional[ECMBoundaryState]) -> Any:
    """Encode the state of the ECM boundary.

    Args:
        boundary: The current state of the boundary to encode
    """
    if boundary is None:
        return None

    return {
        "particles": {
            "par_ids": boundary.particles.par_ids.tolist(),
            "positions": boundary.particles.positions.flatten().tolist(),
            "types": boundary.particles.types.tolist(),
        },
        "bonds": {
            "bond_ids": boundary.bonds.bond_ids.tolist(),
            "group": boundary.bonds.particle_groups.tolist(),
            "types": boundary.bonds.types.tolist(),
        },
        "bonds_types": {
            "bond_type_ids": boundary.bond_types.bond_type_ids.tolist(),
            "k": boundary.bond_types.k.tolist(),
            "r0": boundary.bond_types.r0.tolist(),
        },
        "angles": {
            "angle_cst_ids": boundary.angle_csts.angle_cst_ids.tolist(),
            "group": boundary.angle_csts.particle_groups.tolist(),
            "types": boundary.angle_csts.types.tolist(),
        },
        "angles_types": {
            "angle_type_ids": boundary.angle_cst_types.angle_cst_type_ids.tolist(),
            "k": boundary.angle_cst_types.k.tolist(),
            "t0": boundary.angle_cst_types.t0.tolist(),
        },
    }
    return convert_arrays_to_lists(asdict(boundary))
