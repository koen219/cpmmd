import sys
import numpy as np
import vispy
from vispy import scene, app
import vispy.geometry
from vispy.scene import visuals
from PIL import Image
import click
import yaml
from pathlib import Path


def select_folder() -> Path | None:
    """Opens a folder selection dialog using Qt and returns the selected folder as a Path object."""
    try:
        from PyQt5.QtWidgets import QApplication, QFileDialog
    except ImportError:
        try:
            from PySide6.QtWidgets import QApplication, QFileDialog
        except ImportError:
            raise ImportError(
                "Neither PyQt5 nor PySide6 is installed. Install one using pip."
            )

    # Ensure a QApplication exists (Qt requires one)
    app = QApplication.instance() or QApplication(sys.argv)

    # Show folder selection dialog
    folder = QFileDialog.getExistingDirectory(None, "Select Folder")

    return Path(folder) if folder else None


@click.group()
def main():
    pass


def plot_boxes_at_positions(scatter_data, color):
    box_vertices, box_faces, _ = vispy.geometry.create_box(width=1, height=1, depth=1)
    all_vertices, all_faces = list(), list()
    for i, pos in enumerate(scatter_data):
        # Translate only the 'position' field of the vertices
        translated_vertices = box_vertices["position"] + pos

        # Add other vertex attributes (e.g., normal, color) to the structured array
        translated_vertex_data = np.zeros_like(box_vertices)
        translated_vertex_data["position"] = translated_vertices
        translated_vertex_data["normal"] = box_vertices["normal"]
        translated_vertex_data["texcoord"] = box_vertices["texcoord"]
        translated_vertex_data["color"] = box_vertices["color"]

        all_vertices.append(translated_vertex_data)

        # Adjust the face indices for the current box
        translated_faces = box_faces + i * len(box_vertices)
        all_faces.append(translated_faces)

    # Concatenate all vertices and faces
    all_vertices = np.concatenate(all_vertices)
    all_faces = np.concatenate(all_faces)

    # Create the mesh
    mesh = scene.visuals.Mesh(
        vertices=all_vertices["position"],
        faces=all_faces,
        color=color,
        # shading='smooth',
        # parent=view.scene,
    )
    mesh.set_gl_state(depth_test=True, blend=True, cull_face=False)
    # mesh.set_gl_state(cull_face=False)
    return mesh


def preprocess_ecm(positions, group, particle_types, bonds_types):
    flat_group = group.reshape(-1)
    frame_x = positions[flat_group, 0]
    frame_y = positions[flat_group, 1]
    frame_z = positions[flat_group, 2]
    lines = np.vstack([frame_x, frame_y, frame_z]).T
    return lines, (0.3, 0.3, 0.3, 0.1)


def make_lines(view, positions, group, color):

    flat_group = group.reshape(-1)
    frame_x = positions[flat_group, 0]
    frame_y = positions[flat_group, 1]
    frame_z = positions[flat_group, 2]
    lines = np.vstack([frame_x, frame_y, frame_z]).T

    # lines, colors = make_lines(positions, lines, color)
    poly_visual = visuals.Line(pos=lines, color=color, width=2, connect="segments")
    view.add(poly_visual)


def set_data(view, data, draw_ecm=True):
    cpm, ecm = data
    # Add ECM
    print(set(ecm["bond_types"]))
    polymer = ecm["group"][ecm["bond_types"] == 0]
    crosslink = ecm["group"][ecm["bond_types"] > 0]

    print(f"Number of crosslinks = {crosslink.shape}")

    pos = ecm["positions"]
    if draw_ecm:
        make_lines(view, pos, polymer, (1.0, 1.0, 1.0, 0.7))
        make_lines(view, pos, crosslink, (0.0, 0.7, 0.0, 1.0))
    #    lines, colors = make_lines(pos, polymer, (0.3, 0.3, 0.3, 0.1))
    #    poly_visual = visuals.Line(pos=lines, color=colors, width=2, connect="segments")
    #    view.add(poly_visual)
    #
    #    cross_lines, cross_colors = make_lines(pos, crosslink, (0.0, 0.7, 0.0, 1.0))
    #    cross_visual = visuals.Line(
    #        pos=cross_lines, color=cross_colors, width=2, connect="segments"
    #    )
    #    view.add(cross_visual)

    # Add boundary points
    boundary_pos = ecm["positions"][ecm["particle_types"] == 1]
    boundary = visuals.Markers()
    boundary.set_data(boundary_pos, face_color="yellow", edge_color="yellow", size=4)
    view.add(boundary)

    # Add adhesions
    adhesions_pos = ecm["positions"][ecm["particle_types"] == 2]
    print(adhesions_pos.shape)
    adhesions = visuals.Markers()
    adhesions.set_data(adhesions_pos, face_color="blue", edge_color="blue", size=8)
    view.add(adhesions)

    # Add CPM cell
    scatter_data = np.array(np.nonzero(cpm)).T
    if scatter_data.shape[0] > 0:
        mesh = plot_boxes_at_positions(scatter_data, (1, 0, 0, 0.5))
        view.add(mesh)

    # Add 3D axes
    axes = visuals.XYZAxis(parent=view.scene)
    # Scale the axes to make them visible
    axes.transform = scene.transforms.STTransform(
        scale=(200, 200, 200)
    )  # Adjust scale as needed


def plot_3d_grid(data, time=0, save_as=None):
    cpm, ecm = data
    canvas = scene.SceneCanvas(
        keys="interactive" if save_as is None else None,
        bgcolor="gray",
        show=True if save_as is None else False,
        size=(800, 600),
    )

    # Create a 3D view
    view = canvas.central_widget.add_view()

    # Use an arcball camera for 3D interaction
    camera = scene.cameras.TurntableCamera(fov=60)
    camera.center = (cpm.shape[0] / 2, cpm.shape[1] / 2, cpm.shape[2] / 2)
    camera.distance = 50  # np.sqrt(3 * (max(cpm.shape) / 2) ** 2) * 2
    camera.clip = (0, max(cpm.shape) * 10)
    camera.up = "z"
    camera.azimuth = time * 10.0 / (2 * 3.1415)
    view.camera = camera

    set_data(view, data)

    if save_as:
        # Render the scene to an image array
        image_data = canvas.render()

        # Save the image using PIL
        image = Image.fromarray(image_data)
        image.save(save_as)  # Save as PNG
    else:
        app.run()


@main.command()
@click.argument(
    "filename",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
@click.option("--bgcolor", default="black", type=str)
@click.option("--ecm/--no-ecm", default=True, type=bool)
def interactive(filename, ecm, bgcolor):
    draw_ecm = ecm
    filename = Path(filename)
    data = load_data(filename, 0)
    cpm, ecm = data
    canvas = scene.SceneCanvas(
        keys="interactive",
        bgcolor=bgcolor,
        show=True,
        size=(800, 600),
    )

    # Create a 3D view
    view = canvas.central_widget.add_view()

    # Use an arcball camera for 3D interaction
    camera = scene.cameras.TurntableCamera(fov=60)
    camera.center = (cpm.shape[0] / 2, cpm.shape[1] / 2, cpm.shape[2] / 2)
    camera.distance = 50  # np.sqrt(3 * (max(cpm.shape) / 2) ** 2) * 2
    camera.clip = (0, max(cpm.shape) * 10)
    camera.up = "z"
    # camera.azimuth = time * 10.0 / (2 * 3.1415)
    view.camera = camera
    set_data(view, data, draw_ecm)
    current_time = 0

    def update_scene(event):
        nonlocal current_time
        print(f"Key pressed {event.key}")
        if event.key == "Right":
            current_time += 1
        elif event.key == "Left":
            current_time -= 1
            current_time = max(0, current_time)
        elif event.key == "W":
            folder = select_folder()

            for time in which_times_are_there(filename):
                print(f"Rendering time {time}")
                data = load_data(filename, time)
                for child in list(
                    view.scene.children[2:]
                ):  # Use a copy of the children list
                    child.parent = None
                set_data(view, data, draw_ecm)
                image_data = canvas.render()
                # Save the image using PIL
                image = Image.fromarray(image_data)
                image.save(folder / ("image_" + str(time).zfill(7) + ".png"))
        elif event.key == "Q":
            exit()
        else:
            return

        # for child in list(view.scene.children[2:]):  # Use a copy of the children list
        #     child.parent = None
        try:
            data = load_data(filename, current_time)
        except FileNotFoundError as err:
            print(err)
        else:
            for child in list(
                view.scene.children[2:]
            ):  # Use a copy of the children list
                child.parent = None
            set_data(view, data, draw_ecm)

    canvas.events.key_press.connect(update_scene)

    app.run()


def which_times_are_there(filename: Path):
    filenames = filename.glob("state_*.npz")
    return [
        int(name.stem.removeprefix("state_").removesuffix(".npz")) for name in filenames
    ]


def load_data(filename, time):
    if isinstance(time, int):
        time = str(time).zfill(7)
    with open(filename / "configuration.yaml", "r") as f:
        config = yaml.safe_load(f)
        sizex = config["cpm"]["grid_sizex"]
        sizey = config["cpm"]["grid_sizey"]
        sizez = config["cpm"]["grid_sizez"]
    state = np.load(filename / ("state_" + time + ".npz"))
    cpm = state["cpm"].reshape((sizex, sizey, sizez))
    # Assume data is stored in a NumPy .npy file for simplicity
    # Replace with your data loading logic if needed
    ecm = {
        "group": state["bonds_groups"],
        "positions": state["particle_positions"],
        "particle_types": state["particle_types"],
        "bond_types": state["bonds_types"],
    }
    # cpm = np.load(filename / ("./cpm_" + time + ".npy")).reshape((200, 200, 200))
    return cpm, ecm


@main.command()
@click.argument(
    "folder",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
def video(folder):
    folder = Path(folder)
    for k, data_file in enumerate(Path(folder).glob("*npz")):
        time = data_file.stem.removesuffix(".npz").split("_")[1]
        print("Making picture for time = {time}")
        # Load data
        data = load_data(data_file.parents[0], time)
        time = int(data_file.stem.split("_")[1].split(".")[0])
        # Plot the 3D grid
        plot_3d_grid(
            data, time, save_as=f"{folder}/3d_simulation_" + str(time).zfill(7) + ".png"
        )

    # plot_3d_grid(data=load_data(data_file))


if __name__ == "__main__":

    # interactive(Path(sys.argv[1]))
    # plot_3d_grid(data=load_data(Path(sys.argv[1]), time=str(100).zfill(7)))
    main()
