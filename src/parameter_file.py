from CPMMD.ecm.parameters import (
    EvolutionParameters,
    GenerationParameters,
    InputOutputParameters,
    # CPMParameters,
)

import yaml


from pathlib import Path
from dataclasses import dataclass


@dataclass
class Parameters:
    input_output: InputOutputParameters
    evolution: EvolutionParameters
    generation: GenerationParameters
    cpm: str # as specified in config_spec.yaml in cpmmd


def copy_parfile(destination: Path, file: Path):
    destination = Path(destination) / "configuration.yaml"
    with open(file, "r") as f:
        contents = yaml.safe_load(f)
    with open(destination, "w") as f:
        yaml.safe_dump(contents, f)


def load_file(file: Path):
    with open(file, "r") as f:
        raw_par = list(yaml.safe_load_all(f))[0]

    gen = GenerationParameters(
        box_size_x=raw_par["cpm"]["grid_sizex"],
        box_size_y=raw_par["cpm"]["grid_sizey"],
        Lx=raw_par["cpm"]["grid_sizex"] // 2,
        Ly=raw_par["cpm"]["grid_sizey"] // 2,
        fixed_boundary=raw_par["gen_ecm"]["fixed_boundary"],
        bottom_fixed=raw_par["gen_ecm"]["bottom_fixed"],
        top_fixed=raw_par["gen_ecm"]["top_fixed"],
        contour_length=raw_par["gen_ecm"]["contour_length_fiber"],
        strands=raw_par["gen_ecm"]["number_of_strands"],
        beads=raw_par["gen_ecm"]["number_of_beads_per_strand"],
        spring_r0=raw_par["gen_ecm"]["spring_r0"],
        spring_k=raw_par["gen_ecm"]["spring_k"],
        crosslink_k=raw_par["gen_ecm"]["crosslink_k"],
        helix_angle=raw_par["gen_ecm"]["helix_angle"],
        bend_t0=raw_par["gen_ecm"]["bend_t0"],
        bend_k=raw_par["gen_ecm"]["bend_k"],
        num_init_crosslinks=raw_par["gen_ecm"]["num_init_crosslinks"],
        crosslink_max_r=raw_par["gen_ecm"]["crosslink_max_r"],
        crosslink_quant_step=raw_par["gen_ecm"]["crosslink_quant_step"],
        crosslink_bin_size=raw_par["gen_ecm"]["crosslink_bin_size"],
        network_seed=raw_par["gen_ecm"]["network_seed"],
    )

    io = InputOutputParameters(
        output_folder=raw_par["storage"]["output_folder"],
        output_interval=raw_par["storage"]["stride"],
    )

    evol = EvolutionParameters(
        box_size_x=raw_par["cpm"]["grid_sizex"],
        box_size_y=raw_par["cpm"]["grid_sizey"],
        box_size_z=raw_par["cpm"]["grid_sizez"],
        contour_length=raw_par["gen_ecm"]["contour_length_fiber"],
        md_use_gpu=raw_par["md"]["gpu"],
        md_seed=raw_par["md"]["seed"],
        md_dt=raw_par["md"]["dt"],
        md_its=raw_par["md"]["its"],
        overdamped=raw_par["md"]["overdamped"],
        md_kT=raw_par["md"]["kT"],
        viscosity=raw_par["md"]["viscosity"],
    )
    cpm = yaml.dump(raw_par["cpm"])
#    cpm = CPMParameters(
#        grid_sizex=raw_par["cpm"]["grid_size_x"],
#        grid_sizey=raw_par["cpm"]["grid_size_y"],
#        grid_sizez=raw_par["cpm"]["grid_size_z"],
#        init_cell_size=raw_par["cpm"]["init_cell_size"],
#        init_cell_number_divisions=raw_par["cpm"]["init_cell_number_divisions"],
#        temperature=raw_par["cpm"]["temperature"],
#        target_area_lambda=raw_par["cpm"]["target_area_lambda"],
#        target_area=raw_par["cpm"]["target_area"],
#        j_matrix=raw_par["cpm"]["j_matrix"],
#        num_init_adhesions=raw_par["cpm"]["num_init_adhesions"],
#        init_adhesion_radius=raw_par["cpm"]["init_adhesion_radius"],
#        adhesion_annihilation_penalty=raw_par["cpm"]["adhesion_annihilation_penalty"],
#        adhesion_overflow_number=raw_par["cpm"]["adhesion_overflow_number"],
#        adhesion_overflow_penalty=raw_par["cpm"]["adhesion_overflow_penalty"],
#        check_connectivity=raw_par["cpm"]["check_connectivity"],
#        mcs=raw_par["cpm"]["mcs"],
#    )

    return Parameters(generation=gen, input_output=io, evolution=evol, cpm=cpm)


if __name__ == "__main__":
    import sys

    file = load_file(sys.argv[1])
    print(file)
