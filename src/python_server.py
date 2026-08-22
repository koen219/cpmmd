import socket
import sys
import json
import pandas as pd
import struct
from ecmgen import Network, random_network, single_strand
from CPMMD.ecm.ecm_boundary_state import ECMBoundaryState
from CPMMD.ecm.make_ecm import encode_net_as_dict
from CPMMD.ecm.ecm import MDState, Particles, BondTypes, Bonds, AngleCsts, AngleCstTypes
from CPMMD.ecm.muscle3 import (
    encode_mdstate,
    decode_mdstate,
    decode_cell_ecm_interactions,
    encode_ecm_boundary_state,
)
import numpy as np
from CPMMD.ecm.simulation import Simulation, EvolutionParameters
from dataclasses import dataclass, field
from parameter_file import load_file
from pathlib import Path


def check_boundary(boundary: ECMBoundaryState, sigma):
    pass
#    for pos, typ in zip(boundary.particles.positions, boundary.particles.types):
#        print(pos, typ)
#        if typ == 2:
#            pos = list(map(int, map(np.floor, pos)))
#            print(f"Adh at {pos} heeft spin ")
#            print(np.array(sigma)[*pos])

def test_setup(net_dim, ecmgen_par) -> MDState:
    net =single_strand(
        net_dim[0],
        net_dim[1],
        net_dim[2],
        net_dim[0] / 2-5,
        net_dim[1]/2,
        net_dim[2]/2,
        0.0,
        ecmgen_par.beads,
        ecmgen_par.contour_length
    )
    mdstate = decode_mdstate(encode_net_as_dict(ecmgen_par, net))
    return mdstate


def generate_network(net_dim, ecmgen_par) -> MDState:
    net = random_network(
        net_dim[0],
        net_dim[1],
        ecmgen_par.beads,
        ecmgen_par.strands,
        ecmgen_par.contour_length,
        ecmgen_par.crosslink_max_r,
        ecmgen_par.num_init_crosslinks,
        ecmgen_par.crosslink_bin_size,
        seed=ecmgen_par.network_seed,
        fix_boundary=ecmgen_par.fixed_boundary,
        sizez=net_dim[2],
    )

    mdstate = decode_mdstate(encode_net_as_dict(ecmgen_par, net))
    return mdstate


def receive_large_data(conn):
    # Read the 4-byte size header
    size_data = conn.recv(4)
    size = struct.unpack(">I", size_data)[0]  # ">I" means big-endian unsigned int

    # Read the data in chunks
    data = bytearray()
    while len(data) < size:
        packet = conn.recv(min(size - len(data), 4096))  # Read in chunks of 4096 bytes
        if not packet:
            raise ConnectionError("Connection closed while receiving data")
        data.extend(packet)
    return data


def send_large_data(conn, data):
    # Send the size of the data as a 4-byte header
    conn.sendall(struct.pack(">I", len(data)))
    # Send the actual data
    conn.sendall(data.encode("utf-8"))


def save_state(state: dict, cpm, time, folder):
    if time == 0:
        np.savez_compressed(
            folder + "/state_" + str(time).zfill(7),
            cpm=cpm,
            particle_positions=state["particles"]["positions"],
            particle_types=state["particles"]["types"],
            bonds_groups=state["bonds"]["groups"],
            bonds_types=state["bonds"]["types"],
            allow_pickle=False,
        )
    else:
        np.savez_compressed(
            folder + "/state_" + str(time).zfill(7),
            cpm=cpm,
            particle_positions=state["particles"]["positions"],
            particle_types=state["particles"]["types"],
            # bonds_groups=state["bonds"]["groups"],
            # bonds_types=state["bonds"]["types"],
            allow_pickle=False,
        )


def main(host, port, parfile):
    inout_par = parfile.input_output
    evol_par = parfile.evolution
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as server:
        server.bind((host, port))
        server.listen()
        print("Python server listening on", (host, port))
        ecm = generate_network(
            (evol_par.box_size_x, evol_par.box_size_y, evol_par.box_size_z),
            parfile.generation,
        )
#        ecm = test_setup(
#            (evol_par.box_size_x, evol_par.box_size_y, evol_par.box_size_z),
#            parfile.generation,
#        )

        simulation = Simulation(evol_par, ecm)

        conn, addr = server.accept()
        with conn:
            print("Connected by", addr)

            while True:
                # Receive a message from Rust

                try:
                    data = receive_large_data(conn)  # conn.recv(1024).decode('utf-8')
                except ConnectionError as e:
                    print(f"Connection error: {e}")
                    data = None
                if not data:
                    print("Connection closed by peer. Exiting.")
                    break

                message = json.loads(data)

                if message["type"] == "INIT":
                    data = message["data"]
                    data2 = message["data2"]
                    cell = np.array(message["grid"]).reshape(
                        (
                            evol_par.box_size_x,
                            evol_par.box_size_y,
                            evol_par.box_size_z,
                        )
                    )
#                    with open("ecm_adhesion_zone_dump.data", "w") as file:
#                        np.savetxt(
#                            file, np.array(data["change_type_in_area"]["change_area"])
#                        )

                    interactions = decode_cell_ecm_interactions(data)
                    interactions2 = decode_cell_ecm_interactions(data2)
                    simulation.apply_interactions(interactions)
                    simulation.apply_interactions(interactions2)
                    boundary = simulation.get_boundary_state()

                    # check_boundary(boundary, cell)

                    msg = json.dumps(encode_ecm_boundary_state(boundary))
                    cpm = message["grid"]
                    send_large_data(conn, msg)
                    save_state(
                        encode_mdstate(simulation.get_state()),
                        cpm,
                        0,
                        inout_par.output_folder,
                    )

                    print("Init done")

                elif message["type"] == "STEP":
                    print("Processing step data:")
                    cell = message["grid"]
                    time = message["time"] + 1
                    interactions = decode_cell_ecm_interactions(message["data"])
                    simulation.apply_interactions(interactions)
                    simulation.run()

                    if time % inout_par.output_interval == 0:
                        save_state(
                            encode_mdstate(simulation.get_state()),
                            cell,
                            time,
                            inout_par.output_folder,
                        )
                    boundary = simulation.get_boundary_state()
                    msg = json.dumps(encode_ecm_boundary_state(boundary))
                    send_large_data(conn, msg)

                elif message["type"] == "DONE":
                    print("Simulation completed!")
                    break


if __name__ == "__main__":
    HOST = "127.0.0.1"  # Localhost
    port = int(sys.argv[1])  # Port to listen on
    parfile = load_file(Path(sys.argv[2]))

    main(HOST, port, parfile)
