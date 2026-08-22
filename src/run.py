import subprocess
import time
import socket
import sys
from pathlib import Path
from .parameter_file import load_file, copy_parfile
import click

def wait_for_server(host, port, timeout=30):
    """Wait until a server is available on a given host and port.

    Args:
        host (str): The hostname or IP address.
        port (int): The port number.
        timeout (int): Maximum number of seconds to wait.

    Returns:
        bool: True if the server is available, False if timed out.
    """
    start_time = time.time()
    while True:
        try:
            with socket.create_connection((host, port), timeout=2):
                print(f"Server on {host}:{port} is ready!")
                return True
        except (ConnectionRefusedError, socket.timeout):
            if time.time() - start_time > timeout:
                print(f"Timeout reached: Server on {host}:{port} did not start within {timeout} seconds.")
                return False
            time.sleep(0.5)  # Wait a short period before retrying


def find_free_port():
    with socket.socket() as s:
        s.bind(("", 0))  # Bind to a free port provided by the host.
        return s.getsockname()[1]


@click.command()
@click.argument("parfile")
def main(parfile):
    output_folder = load_file(parfile).input_output.output_folder
    Path(output_folder).mkdir(exist_ok=True)
    copy_parfile(output_folder, parfile)

    port = find_free_port()
    ecm_code_output = open(f"{output_folder}/ecm.txt", "w")
    ecm_code = subprocess.Popen(
        [
            (Path(__file__).resolve().parent / "../.venv/bin/python3").__str__(),
            (Path(__file__).resolve().parent / "python_server.py").__str__(),
            str(port),
            parfile,
        ],
        stdout=ecm_code_output,
        stderr=ecm_code_output,
    )
    print("Start ecm!")
    time.sleep(5)

#    if not wait_for_server("127.0.0.1", port, timeout=30):
#        print("Python server did not start in time. Exiting.")
#        sys.exit(1)


    connection_adress = f"127.0.0.1:{port}"
    print(port, connection_adress)
    cpm_code_output = open(f"{output_folder}/cpm.txt", "w")
    cpm_code = subprocess.Popen(
        [
            ## "../venv/bin/python3",
            (Path(__file__).resolve().parent / "../.venv/bin/python3").__str__(),
            (Path(__file__).resolve().parent / "rust_server.py").__str__(),
            parfile,
        ],
        env=dict(
            CONNADRESS=connection_adress,
            RUST_BACKTRACE="1",
        ),
        stdout=cpm_code_output,
        stderr=cpm_code_output,
    )

    cpm_code.wait()
    cpm_code_output.flush()
    ecm_code_output.flush()
    cpm_code_output.close()
    ecm_code_output.close()
    print("Simulation complete!")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        parfile = Path(sys.argv[1])
        if not parfile:
            raise FileNotFoundError(f"Parameterfile {sys.argv[1]} not found.")
    else:
        raise RuntimeError(
            'No parameterfile supplied. Please run as "python run.py NAMEOFPARFILE"'
        )

    main(sys.argv[1])
