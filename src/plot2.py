import sys
import click
import json
import yaml
import numpy as np
from pathlib import Path
from PIL import Image
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import vispy
from vispy import scene, app
from vispy.scene import visuals
import vispy.geometry
from scipy.ndimage import convolve

# ----------------------------------------------------------------------
# Utility Functions
# ----------------------------------------------------------------------
def select_folder(default_folder) -> Path | None:
    """Opens a folder selection dialog using Qt/PySide and returns the selected folder as a Path object."""
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
    qt_app = QApplication.instance() or QApplication(sys.argv)

    # Show folder selection dialog
    folder = QFileDialog.getExistingDirectory(None, "Select Folder", directory= str(default_folder))
    return Path(folder) if folder else None

def select_time(minValue,maxValue) -> Path | None:
    from PyQt5 import QtWidgets

    # qt_app = QApplication.instance() or QApplication(sys.argv)

    time, ok = QtWidgets.QInputDialog().getInt(None, "Time to jump to:",  f"{minValue=}\n{maxValue=}", min=minValue, max=maxValue)

    if ok:
        return time

    return None

_open_editor = None
def ask_for_drawing_options(what_to_draw):
    global _open_editor
    from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QCheckBox, QPushButton, QLabel, QSlider
    from PyQt5.QtCore import Qt
    class OptionsEditor(QWidget):
        def __init__(self, d):
            super().__init__()
            self.setWindowTitle("Options")
            self.d = d
            self.checkboxes = dict()
            self.sliders = dict()

            layout = QVBoxLayout(self)
            for key, value in d.items():
                if isinstance(value, bool):
                    name = key.removeprefix("draw_").capitalize()
                    cb = QCheckBox(name)
                    cb.setChecked(value)
                    layout.addWidget(cb)
                    self.checkboxes[key] = cb
                if isinstance(value, dict) and value['type'] == 'slider':
                    name = key
                    label = QLabel(f"{name}: {value.get('value', 0)}")
                    slider = QSlider(Qt.Horizontal)
                    slider.setMinimum(value.get('min', 0))
                    slider.setMaximum(value.get('max', 100))
                    slider.setValue(value.get('value', 0))
                    slider.valueChanged.connect(lambda val, l=label, n=name: l.setText(f"{n}: {val}"))
    
                    slider_layout = QVBoxLayout()
                    slider_layout.addWidget(label)
                    slider_layout.addWidget(slider)
                    layout.addLayout(slider_layout)
                    self.sliders[key] = slider


            ok_button = QPushButton("OK")
            ok_button.clicked.connect(self.update_dict)
            layout.addWidget(ok_button)

        def update_dict(self):
            for key, cb in self.checkboxes.items():
                self.d[key] = cb.isChecked()
            for key, slider in self.sliders.items():
                self.d[key]['value'] = slider.value()
            self.close()

    _open_editor = OptionsEditor(what_to_draw)
    _open_editor.show()

def select_file(default_folder) -> Path | None:
    """Opens a file selection dialog using Qt/PySide and returns the selected file as a Path object."""
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
    qt_app = QApplication.instance() or QApplication(sys.argv)

    # Show folder selection dialog
    file, _ = QFileDialog.getSaveFileName(None, "Select File", directory= str(default_folder), filter="Images (*.png *.xpm *.jpg)")
    return Path(file) if file else None

def ask_for_analysis(folder):
    from PyQt5.QtWidgets import QInputDialog, QMessageBox
    commands = {
        "cellshape/Volume",
        "cellshape/Surface",
        "count-beads/Adhesion",
        "network/perculation"
    }

    choice, ok = QInputDialog.getItem(None, "Select Command", "Subcommand:", commands, 0, False)
    if not ok:
        return None
    command, what = choice.split('/')

    print(["cpmmdplot2", "analyse", folder, "--what", command, "--yvalue", what, "--force"])
    subprocess.Popen(["cpmmdplot2", "analyse", folder, "--what", command, "--yvalue", what, "--force"])
    
    return command, ok

def which_times_are_there(folder: Path):
    """Returns a list of time steps found in the given folder (based on 'state_*.npz' files)."""
    filenames = folder.glob("state_*.npz")
    return [
        int(name.stem.removeprefix("state_").removesuffix(".npz")) for name in filenames
    ]


def load_data(folder: Path, time):
    """
    Loads the CPM data and ECM data from 'state_<time>.npz' in the specified folder.
    Also reads grid sizes from 'configuration.yaml'.
    """
    if isinstance(time, int):
        time = str(time).zfill(7)
    with open(folder / "configuration.yaml", "r") as f:
        config = yaml.safe_load(f)
        sizex = config["cpm"]["grid_sizex"]
        sizey = config["cpm"]["grid_sizey"]
        sizez = config["cpm"]["grid_sizez"]
    state0 = np.load(folder / f"state_{'0'.zfill(7)}.npz")
    state = np.load(folder / f"state_{time}.npz")
    cpm = state["cpm"].reshape((sizex, sizey, sizez))
    ecm = {
        "group": state0["bonds_groups"],
        "positions": state["particle_positions"],
        "particle_types": state["particle_types"],
        "bond_types": state0["bonds_types"],
    }
    metadata = config
   #  dict(
   #     shape=(sizex,sizey,sizez),
   #     time_stride=time_stride,
   #     )
    return cpm, ecm, metadata


# ----------------------------------------------------------------------
# Artist class for drawing 3D elements
# ----------------------------------------------------------------------
class Artist:
    """
    Encapsulates the methods needed to draw various 3D elements such as:
      - Boxes for the CPM cells
      - Lines for ECM (polymer vs. crosslinks)
      - Markers for boundary and adhesion particles
      - Axes
    """
    def __init__(self, view):
        """
        :param view: A vispy.scene.widgets.viewbox.ViewBox object,
                     typically obtained from `canvas.central_widget.add_view()`.
        """
        self.view = view

    def draw_cpm_cells(self, cpm_data, color=(1, 0, 0, 0.5)):
        """
        Draw boxes at positions given by nonzero entries in the cpm_data.
        :param cpm_data: A 3D numpy array representing the CPM grid.
        :param color: A tuple (r, g, b, a) for color.
        """
        # Find all non-zero indices
        scatter_data = np.array(np.nonzero(cpm_data)).T
        if scatter_data.shape[0] == 0:
            return  # nothing to draw

        # Create and add a mesh of boxes
        mesh = self._plot_boxes_at_positions(scatter_data, color)
        self.view.add(mesh)

    def draw_ecm_lines(
        self,
        positions: np.ndarray,
        group: np.ndarray,
        bond_types: np.ndarray,
        polymer_color=(1.0, 1.0, 1.0, 0.7),
        crosslink_color=(0.0, 0.7, 0.0, 1.0),
    ):
        """
        Draw lines representing ECM polymer vs. crosslinks.
        :param positions: (N,3) array of all particle positions.
        :param group: (N,2) or (N,...) array grouping indices for line connectivity.
        :param bond_types: array indicating type of bond.
        :param polymer_color: RGBA for polymer lines.
        :param crosslink_color: RGBA for crosslink lines.
        """
        polymer = group[bond_types == 0]
        crosslink = group[bond_types > 0]

        # polymer lines
        self._make_lines(positions, polymer, polymer_color)

        # crosslink lines
        self._make_lines(positions, crosslink, crosslink_color)

    def draw_boundary_points(
        self, positions: np.ndarray, particle_types: np.ndarray, color="yellow"
    ):
        """
        Draw boundary points (where `particle_types == 1`) as markers.
        :param positions: (N,3) array of all particle positions.
        :param particle_types: array of types for each particle.
        :param color: color for boundary markers.
        """
        boundary_pos = positions[particle_types == 1]
        if boundary_pos.shape[0] > 0:
            boundary = visuals.Markers()
            boundary.set_data(boundary_pos, face_color=color, edge_color=color, size=4)
            self.view.add(boundary)

    def draw_adhesion_points(
        self, positions: np.ndarray, particle_types: np.ndarray, color="blue"
    ):
        """
        Draw adhesion points (where `particle_types == 2`) as markers.
        :param positions: (N,3) array of all particle positions.
        :param particle_types: array of types for each particle.
        :param color: color for adhesion markers.
        """
        adhesions_pos = positions[particle_types == 2]
        if adhesions_pos.shape[0] > 0:
            adhesions = visuals.Markers()
            adhesions.set_data(
                adhesions_pos, face_color=color, edge_color=color, size=8
            )
            self.view.add(adhesions)

    def draw_axes(self, scale=200):
        """
        Draw 3D axes to help visualize orientation.
        :param scale: Scaling factor for the axes length.
        """
        axes = visuals.XYZAxis(parent=self.view.scene)
        # Scale them to make them visible
        axes.transform = scene.transforms.STTransform(scale=(scale, scale, scale))

    # ----------------------------------------------------------------------
    # Internal helper methods
    # ----------------------------------------------------------------------
    def _plot_boxes_at_positions(self, scatter_data, color):
        box_vertices, box_faces, _ = vispy.geometry.create_cube() # vispy.geometry.create_box(width=1, height=1, depth=1)
        all_vertices, all_faces = [], []

        for i, pos in enumerate(scatter_data):
            # Translate only the 'position' field of the vertices
            translated_vertices = box_vertices["position"] + pos

            # Add other attributes (e.g., normal, color) to the structured array
            translated_vertex_data = np.zeros_like(box_vertices)
            translated_vertex_data["position"] = translated_vertices
            translated_vertex_data["normal"] = box_vertices["normal"]
            translated_vertex_data["texcoord"] = box_vertices["texcoord"]
            translated_vertex_data["color"] = box_vertices["color"]

            all_vertices.append(translated_vertex_data)

            # Adjust face indices for the current box
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
            shading='flat',
        )
        mesh.set_gl_state(depth_test=True, blend=True, cull_face=False)
        return mesh

    def _make_lines(self, positions, group, color):
        """
        Creates line segments between pairs of points specified in 'group'.
        :param positions: (N,3) array of all particle positions.
        :param group: Indices that define start-end pairs for lines.
        :param color: RGBA or named color.
        """
        if group.size == 0:
            return
        flat_group = group.reshape(-1)
        frame_x = positions[flat_group, 0]
        frame_y = positions[flat_group, 1]
        frame_z = positions[flat_group, 2]
        lines = np.vstack([frame_x, frame_y, frame_z]).T

        line_visual = visuals.Line(pos=lines, color=color, width=2, connect="segments")
        self.view.add(line_visual)


# ----------------------------------------------------------------------
# SceneController class: manages time stepping and key events
# ----------------------------------------------------------------------
class SceneController:
    """
    Manages the 3D scene, including:
      - current time index
      - responding to key presses
      - loading & drawing data via an Artist
    """
    def __init__(
        self,
        folder: Path,
        canvas: scene.SceneCanvas,
        draw_ecm=True,
        draw_adh=True,
        draw_boundary = True,
        draw_cpm=True,
        cpm_color=(1, 0, 0, 0.5),
        ecm_polymer_color=(1.0, 1.0, 1.0, 0.7),
        ecm_crosslink_color=(0.0, 0.7, 0.0, 1.0),
        boundary_color="yellow",
        adhesion_color="blue",
        start_time=0,
        time_stride=None, # None uses the time from configuration.yaml
        only_adh_bonds=False,
    ):
        self.folder = folder
        self.canvas = canvas
        self.view = canvas.central_widget.add_view()
        self.what_to_draw = dict(
            draw_ecm = draw_ecm,
            draw_cpm = draw_cpm,
            draw_adh = draw_adh,
            draw_boundary = draw_boundary,
            only_adh_bonds=only_adh_bonds,
            draw_distance = dict(value=0,min=0,max=100, type='slider'),
        )
        self.cpm_color = cpm_color
        self.ecm_polymer_color = ecm_polymer_color
        self.ecm_crosslink_color = ecm_crosslink_color
        self.boundary_color = boundary_color
        self.adhesion_color = adhesion_color
        # self.only_adh_bonds = only_adh_bonds
        self._bonds_to_draw = None
        self._draw_axis = True

        if start_time < 0:
            start_time = max(self._times())
        self.current_time = start_time
        _, _, meta_data = load_data(self.folder, 0)
        self._meta_data = meta_data
        if time_stride == None:
            self.time_stride = meta_data['storage']['stride']
        else:
            self.time_stride = time_stride
        self._configure_camera()
        self.artist = Artist(self.view)

        self._saved_camera_state = None

        # Initially draw the scene
        self.draw_scene()
    
    def _times(self):
            return which_times_are_there(self.folder)

    def _configure_camera(self):
        """Set up the camera for 3D interaction."""
        # We start by loading time 0 just for shape references
        cpm_data, _, _  = load_data(self.folder, 0)
        camera = scene.cameras.TurntableCamera(fov=60)
        camera.center = (
            cpm_data.shape[0] / 2,
            cpm_data.shape[1] / 2,
            cpm_data.shape[2] / 2,
        )
        camera.distance = 50
        camera.clip = (0, max(cpm_data.shape) * 10)
        camera.up = "z"
        self.view.camera = camera

    def on_key_press(self, event):
        """
        Respond to keyboard events:
          Right/Left arrows to step forward/backward in time,
          W to export images for all available times,
          Q to quit.
        """
        key = event.key.name if hasattr(event.key, "name") else event.key
        if key == "Right":
            # self.current_time += 1
            self.current_time = min(max(which_times_are_there(self.folder)), self.current_time + self.time_stride)
            print(f"Time = {self.current_time}")
            self.draw_scene()
        elif key == "Left":
            self.current_time = max(0, self.current_time - self.time_stride)
            print(f"Time = {self.current_time}")
            self.draw_scene()
        elif key == 'R':
            self._bonds_to_draw = None
            self.draw_scene()
        elif key == "T":
            time = select_time( min(self._times()), max(self._times()))

            if time is not None:
                self.current_time = time
                self.draw_scene()
        elif key == "A":
            ask_for_analysis(self.folder)
        elif key == "O":
            ask_for_drawing_options(self.what_to_draw)
        elif key == 'S':
            save_name = select_file(self.folder)
            if save_name is None:
                print("No file selected. Aborting.")
                return 
            save_name = save_name if str(save_name).endswith('.png') else str(save_name) + '.png'
            image_data = self.canvas.render()
            image = Image.fromarray(image_data)
            image.save(save_name) 
            print(f"Image saved as {save_name}")
        elif key == "W":
            if self._saved_camera_state is None:
                self._saved_camera_state = self._get_camera_state()
                print(f"Saved Camera State")
                return
            else:
                which_times = tuple(map(int, input("Time interval of movie like this: \"T1/T2\". Use T2 = -1 for the final frame. Use T1/T2/frames to set number of frames per timepoint > 1.").split('/')))
                assert len(which_times)>=2,f"Input parsed to {which_times}. Please correct input"
                tmin = which_times[0]
                tmax = which_times[1] if which_times[1] > 0 else max(which_times_are_there(self.folder))
                if len(which_times) > 2:
                    number_of_frames_per_time = which_times[2]
                else:
                    number_of_frames_per_time = 1

                folder_out = select_folder(self.folder)
                if folder_out is None:
                    print("No folder selected. Aborting animation.")
                    self._saved_camera_state = None 
                    return 
            self._export_all_images(number_of_frames_per_time, tmin,tmax, folder_out, self._saved_camera_state, self._get_camera_state())
            self._saved_camera_state = None 
        elif key == "Q":
            exit()

    def _get_camera_state(self):
        """Extract the current camera state as a dict."""
        return self.view.camera.get_state()

    def _identify_fibers_without_excluded(self, ecm_data, possible_bonds):
        excluded,  = np.where(ecm_data['particle_types'] == 3)
        bonds_indices_not_touching_excluded, = np.where(np.logical_and(
            np.isin(possible_bonds[:, 0], excluded, invert=True),
            np.isin(possible_bonds[:, 1], excluded, invert=True)
        ))
        b = self._meta_data['gen_ecm']['number_of_beads_per_strand']
        return bonds_indices_not_touching_excluded

    def _identify_fibers_with_adhesion(self, ecm_data):
        adhesions,  = np.where(ecm_data['particle_types'] == 2)
        strands = self._meta_data['gen_ecm']['number_of_strands']
        beads = self._meta_data['gen_ecm']['number_of_beads_per_strand']

        bonds_indices_touching_adhesion, = np.where(np.logical_or(
            np.isin(ecm_data['group'][:, 0], adhesions),
            np.isin(ecm_data['group'][:, 1], adhesions)
        ))
        fiber_indices =np.unique( np.floor(bonds_indices_touching_adhesion / (beads - 1)).astype(int) )
        bonds = np.unique( ((beads - 1) * fiber_indices[:, None] + np.arange(0, beads)).flatten().astype(int) )
        
        adh_bonds = ecm_data['group'][bonds, :]

        temp = np.unique(ecm_data['bond_types'][:(beads-1) * strands])
        assert len(temp) == 1, temp

        return bonds_indices_touching_adhesion
    
    def _select_selection_of_fibers(self, ecm_data, percent):
        print(ecm_data['group'].shape)
        print(self._meta_data['gen_ecm']['number_of_beads_per_strand'] , self._meta_data['gen_ecm']['number_of_strands'] )
        print( (self._meta_data['gen_ecm']['number_of_beads_per_strand']-1) * self._meta_data['gen_ecm']['number_of_strands'] )
        bonds_per_fiber = self._meta_data['gen_ecm']['number_of_beads_per_strand'] - 1
        total_fibers = self._meta_data['gen_ecm']['number_of_strands']
        number_to_select = max(1, int(total_fibers * percent))

        selected_fibers = np.random.choice(total_fibers, number_to_select, replace=False)
        selected_fibers.sort()

        bonds_indices = (selected_fibers[:, None] * bonds_per_fiber + np.arange(bonds_per_fiber)).flatten()
        print(np.max(bonds_indices))

        return ecm_data['group'][bonds_indices]

    def _select_bonds_that_moved(self, ecm_data, draw_distance):
        from TST_analyse import RSimulation, displacement_from_start
        sim = RSimulation(self.folder)
        df = displacement_from_start(sim).query(f"time == {self.current_time}")
        bead_index = df.query(f"distance >= {draw_distance}")['bead_index'].astype(int)
        return np.where(np.logical_or(
            np.isin(ecm_data['group'][:, 0], bead_index),
            np.isin(ecm_data['group'][:, 1], bead_index)))[0]

    def draw_scene(self):
        """Clear the existing scene objects and redraw them for the current time."""

        print(self._get_camera_state())

        # Remove old children (except camera, which is child index 0 or 1)
        for child in list(self.view.scene.children[2:]):
            child.parent = None

        try:
            cpm_data, ecm_data, _ = load_data(self.folder, self.current_time)
        except FileNotFoundError as e:
            print(e)
            return


        # Draw ECM lines
        if self.what_to_draw['draw_ecm']:
            if self.what_to_draw['only_adh_bonds']:
                self._bonds_to_draw = self._identify_fibers_with_adhesion(ecm_data)
            else:
                if self._bonds_to_draw is None:
                    self._bonds_to_draw = self._select_bonds_that_moved(ecm_data, self.what_to_draw['draw_distance']['value'])
                    print(self._bonds_to_draw)
                    # self._bonds_to_draw = np.arange(ecm_data['group'].shape[0])
            beads = self._meta_data['gen_ecm']['number_of_beads_per_strand']
            strands = self._meta_data['gen_ecm']['number_of_strands']
            # print(ecm_data['group'][:(beads-1)*strands].tolist())
            self.artist.draw_ecm_lines(
                positions=ecm_data["positions"],
                group=ecm_data["group"][self._bonds_to_draw, :],
                bond_types=ecm_data["bond_types"][self._bonds_to_draw],
                polymer_color=self.ecm_polymer_color,
                crosslink_color=self.ecm_crosslink_color,
            )
        if self.what_to_draw['draw_boundary']:
            self.artist.draw_boundary_points(
                ecm_data["positions"], ecm_data["particle_types"], self.boundary_color
            )
        if self.what_to_draw['draw_adh']:
            self.artist.draw_adhesion_points(
                ecm_data["positions"], ecm_data["particle_types"], self.adhesion_color
            )

        # Draw CPM cells
        if self.what_to_draw['draw_cpm']:
            cpm = cpm_data
            kernel = np.ones((3, 3, 3), dtype=np.uint8)
            kernel[1, 1, 1] = 0
            zero_mask = (cpm == 0).astype(np.uint8)
            zero_neighbour_count = convolve(zero_mask, kernel, mode="constant", cval=0)
            border = np.logical_and(cpm == 1, zero_neighbour_count > 0)
            self.artist.draw_cpm_cells(border, color=(1.0, 0.0, 0.0, 0.5) )# self.cpm_color)

        # Draw axes
        if self._draw_axis:
            self.artist.draw_axes(scale=200)
        

    def _export_all_images(self, frames_per_time, tmin, tmax, folder_out, start_state, end_state):
        """Open a folder dialog and save images for all time steps."""
        tmin /= self.time_stride
        tmax /= self.time_stride
        state = start_state.copy()
        times = sorted(which_times_are_there(self.folder))
        for t,_ in enumerate(times):
            if t < tmin or t > tmax:
                continue
            for j in range(frames_per_time):
                print(f"Rendering time {t * self.time_stride} frame {j+1}/{frames_per_time}")
                frame = frames_per_time * t + j
                i = frame / (frames_per_time * (tmax - tmin))
                new_scale_factor = (1 - i) * start_state['scale_factor'] + i * end_state['scale_factor']
                new_elevation = (1 - i) * start_state['elevation'] + i * end_state['elevation']
                new_azimuth = (1 - i) * start_state['azimuth'] + i * end_state['azimuth']
                state['scale_factor'] = new_scale_factor
                state['elevation'] = new_elevation
                state['azimuth'] = new_azimuth


                # Update the camera state.
                camera = self.view.camera
                camera.set_state(state)
#                camera.scale_factor = new_scale_factor
#                camera.elevation = new_elevation
#                camera.azimuth = new_azimuth

                self.current_time = t * self.time_stride
                self.draw_scene()

                image_data = self.canvas.render()
                image = Image.fromarray(image_data)
                image.save(folder_out / f"image_{str(t).zfill(7)}_frame_{str(frame).zfill(7)}.png")


# ----------------------------------------------------------------------
# High-level function for single-snapshot rendering
# ----------------------------------------------------------------------
def plot_3d_grid(data, time=0, save_as=None,
                 draw_ecm=True, draw_cpm=True,
                 cpm_color=(1,0,0,0.5),
                 ecm_polymer_color=(1.0, 1.0, 1.0, 0.7),
                 ecm_crosslink_color=(0.0, 0.7, 0.0, 1.0),
                 boundary_color="yellow",
                 adhesion_color="blue"):
    """
    Creates a 3D scene, adds the data, and optionally saves the final render as an image.
    :param data: A tuple (cpm, ecm) containing the 3D data.
    :param time: Used to set camera orientation or for naming output images.
    :param save_as: File path to save the rendered image, or None for interactive mode.
    :param draw_ecm: Whether to draw ECM lines.
    :param draw_cpm: Whether to draw CPM boxes.
    :param cpm_color, ecm_polymer_color, ecm_crosslink_color, boundary_color, adhesion_color:
                    Various color parameters.
    """
    cpm_data, ecm_data = data

    # Create a SceneCanvas
    canvas = scene.SceneCanvas(
        keys="interactive" if save_as is None else None,
        bgcolor="gray",
        show=(save_as is None),
        size=(800, 600),
    )

    view = canvas.central_widget.add_view()
    camera = scene.cameras.TurntableCamera(fov=60)
    camera.center = (
        cpm_data.shape[0] / 2,
        cpm_data.shape[1] / 2,
        cpm_data.shape[2] / 2,
    )
    camera.distance = 50
    camera.clip = (0, max(cpm_data.shape) * 10)
    camera.up = "z"
    camera.azimuth = time * 10.0 / (2 * 3.1415)
    view.camera = camera

    # Use our Artist class
    artist = Artist(view)

    if draw_ecm:
        artist.draw_ecm_lines(
            positions=ecm_data["positions"],
            group=ecm_data["group"],
            bond_types=ecm_data["bond_types"],
            polymer_color=ecm_polymer_color,
            crosslink_color=ecm_crosslink_color,
        )
        artist.draw_boundary_points(
            ecm_data["positions"], ecm_data["particle_types"], color=boundary_color
        )
        artist.draw_adhesion_points(
            ecm_data["positions"], ecm_data["particle_types"], color=adhesion_color
        )

    if draw_cpm:
        artist.draw_cpm_cells(cpm_data, color=cpm_color)

    artist.draw_axes(scale=200)

    if save_as:
        # Render the scene to an image array
        image_data = canvas.render()
        image = Image.fromarray(image_data)
        image.save(save_as)
    else:
        app.run()


# ----------------------------------------------------------------------
# CLI Definition
# ----------------------------------------------------------------------
@click.group()
def cli():
    pass


@cli.command()
@click.argument(
    "folder",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
@click.option("--bgcolor", default="black", type=str, help="Background color for the scene.")
@click.option("--draw-ecm/--no-draw-ecm", default=True, help="Toggle drawing ECM lines.")
@click.option("--draw-adh/--no-draw-adh", default=True, help="Toggle drawing adhesion markers.")
@click.option("--draw-boundary/--no-draw-boundary", default=True, help="Toggle drawing boundary markers.")
@click.option("--draw-cpm/--no-draw-cpm", default=True, help="Toggle drawing CPM boxes.")
@click.option("--cpm-color", default="red", help="Color for CPM boxes (RGBA or name).")
@click.option("--ecm-polymer-color", default="white", help="Color for ECM polymer lines.")
@click.option("--ecm-crosslink-color", default="green", help="Color for ECM crosslink lines.")
@click.option("--boundary-color", default="yellow", help="Color for boundary markers.")
@click.option("--adhesion-color", default="blue", help="Color for adhesion markers.")
@click.option("--start-time", "--draw-time", "--time", default =0, type=int, help="Time point to draw")
@click.option("--time-stride", default=None, type=int, help="Stride between timepoints.")
@click.option("--only-adh-bonds", default =False, is_flag=True, help="Draw only bonds connected to an adhesion")
def interactive(
    folder,
    bgcolor,
    draw_ecm,
    draw_adh,
    draw_boundary,
    draw_cpm,
    cpm_color,
    ecm_polymer_color,
    ecm_crosslink_color,
    boundary_color,
    adhesion_color,
    start_time,
    time_stride,
    only_adh_bonds,
):
    """
    Launch an interactive visualization of the 3D data in the given FOLDER.
    Use arrow keys to move in time, 'W' to write out images, 'Q' to quit.
    """
    folder = Path(folder)

    canvas = scene.SceneCanvas(
        keys="interactive",
        bgcolor=bgcolor,
        show=True,
        size=(800, 600),
    )

    # Create the controller (which also draws the initial scene)
    controller = SceneController(
        folder=folder,
        canvas=canvas,
        draw_ecm=draw_ecm,
        draw_cpm=draw_cpm,
        draw_adh=draw_adh,
        draw_boundary = draw_boundary,
        cpm_color=cpm_color,
        ecm_polymer_color=ecm_polymer_color,
        ecm_crosslink_color=ecm_crosslink_color,
        boundary_color=boundary_color,
        adhesion_color=adhesion_color,
        start_time=start_time,
        only_adh_bonds=only_adh_bonds,
        time_stride=time_stride,
    )

    # Connect the controller to the canvas's key_press event
    canvas.events.key_press.connect(controller.on_key_press)
    app.run()


@cli.command()
@click.argument(
    "folder",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True
    ),
)
@click.option("--draw-ecm/--no-draw-ecm", default=True, help="Toggle drawing ECM lines.")
@click.option("--draw-cpm/--no-draw-cpm", default=True, help="Toggle drawing CPM boxes.")
def video(folder, draw_ecm, draw_cpm):
    """
    Renders all states in FOLDER to a series of PNG images (for later compilation into a video).
    """
    folder = Path(folder)
    for data_file in sorted(folder.glob("state_*.npz")):
        time_str = data_file.stem.split("_")[1]  # e.g. "0000000"
        print(f"Making picture for time = {time_str}")
        cpm, ecm, _ = load_data(folder, time_str)
        out_name = folder / f"3d_simulation_{time_str}.png"

        plot_3d_grid(
            (cpm,ecm),
            time=int(time_str),
            save_as=str(out_name),
            draw_ecm=draw_ecm,
            draw_cpm=draw_cpm,
            cpm_color=(1, 0, 0, 0.5),
            ecm_polymer_color=(1.0, 1.0, 1.0, 0.7),
            ecm_crosslink_color=(0.0, 0.7, 0.0, 1.0),
        )

@cli.command()
@click.argument("folder", nargs=-1,type=click.Path(file_okay=False,dir_okay=True,exists=True, readable=True))
@click.option("--par", "-p", multiple=True, type=str, help="Parameters which are taken for the simulation")
@click.option("--what", "-w", multiple=False, type=str, help="Which tool to use.")
@click.option("--yvalue", "-y", multiple=False, type=str, help="Which yvalue to plot (depending on the --what parameter)")
@click.option("--hue", "-h", multiple=False, type=str, help="Which hue is used for plotting")
@click.option("--plot/--no-plot", is_flag=True, default=True, help="To create a plot")
@click.option("--show/--no-show", is_flag=True, default=True, help="If plot is created, toggle to show plot.")
@click.option("--force", '-f', is_flag=True, default=False, help="Overwrite existing analyse files in datadir folder")
@click.option("--savename", multiple=False, type=click.Path(file_okay=True, dir_okay=False, exists=False), help="Savename of plot")
@click.option( "--axispar", multiple=False, default=dict(), help='Axis parameters as JSON dicts, e.g. \'{"xlabel": "Time"}\'')
@click.option( "--snskwargs", multiple=False, default=dict(), help='Seaborn parameters as JSON dicts, e.g. \'{"xlabel": "Time"}\'')
@click.pass_context
def analyse(ctx, folder, par, what, **kwargs): #  yvalue, hue, plot, show, force, savename, axispar):
    """Used to analyse simulations using TST_analyse package (if installed)"""
    from TST_analyse import RSimulation, Experiment

    Experiment.Simulation = RSimulation

    if len(par) == 0:
        par = None
    kwargs['axispar'] = json.loads(kwargs['axispar'])
    kwargs['snskwargs'] = json.loads(kwargs['snskwargs'])
    exper = Experiment.fromIter(folder, parameters=par)

    functions = {
            'network': lambda: network(exper, **kwargs),
            'cellshape': lambda: cellshape(exper, **kwargs),
            'count-beads': lambda: count_beads(exper, **kwargs)
            }
    if what in functions.keys():
        functions[what]()
    else:
        print("Error: what parameter not valid: ", what, "\nOptions are: ", list(functions.keys()))
#    if what == 'network':
#        network(exper, **kwargs)
#    elif what == 'cellshape':
#        cellshape(exper, **kwargs) # yvalue, hue, plot, show, force, savename, axispar, snskwargs)
#    elif what == "count-beads":
#        count_beads(exper, **kwargs) # yvalue, hue, plot, show, force, savename, axispar, snskw)
#    else:
        #print("Error: what parameter not valid: ", what)


def plot_data(df, x, y, hue, plot, savename, show, axispar, snskwargs=dict()):
    if plot:
        ax = sns.relplot(
            data=df,
            x=x,
            y=y,
            hue=hue,
            **snskwargs
        )
        if axispar:
            ax.set(**axispar)
        if savename is not None:
            plt.savefig(savename, dpi=300)
        if show:
            plt.show()

def cellshape(exper, yvalue, hue, plot, show, force, savename, axispar, **kwargs):
    from TST_analyse import perimiter

    if force:
        perimiter.overwrite = True
    df = exper.GetSummary(perimiter, save=not plot)

    y = {
        'volume': 'area',
        'surface': 'perimiter',
    }.get(yvalue.lower())

    plot_data(df, 'time', y, hue, plot, savename, show, axispar, **kwargs)

def count_beads(exper, yvalue, hue, plot, show, force, savename, axispar, **kwargs):
 
    from TST_analyse import number_of_beads
 
    number_of_beads.kwargs = dict(which_type={
             'free': 0,
             'boundary': 1,
             'adhesion': 2,
             'excluded': 3,
    }.get(yvalue.lower()))
    number_of_beads.overwrite = True
    df = exper.GetSummary(number_of_beads, save=False)
    df[yvalue] = df['number']
    df.drop('number', axis=1)
    plot_data(df, 'time', yvalue, hue, plot, savename, show, axispar, **kwargs)

 
def network(exper, yvalue, hue, plot, show, force, savename, axispar, **kwargs):

    from TST_analyse import perculation, displacement_from_start
    from TST_analyse.inputoutput.muscle_faker import MakeAFake

    if 'perculation' in yvalue:
        exper._simulations = [MakeAFake(sim) for sim in exper._simulations ]
        if force:
            perculation.overwrite = True
        df = exper.GetSummary(perculation, save=True)
        print("-"*10 + "Network Perculation" + 10*"-")
        print(df)
        print("-"*39)
    elif 'displacement' == yvalue:
        if force:
            displacement_from_start.overwrite = True
        df = exper.GetSummary(displacement_from_start)
        print(df)
if __name__ == "__cli__":
    cli()

"""


def plot_series(x, y, axis_params=None, label=None, title=None, save_path=None):
    fig, ax = plt.subplots(1, 1)
    ax.plot(x, y, label=label)

    if axis_params:
        ax.set(**axis_params)
    if title:
        ax.set_title(title)
    if label:
        ax.legend()

    if save_path:
        fig.savefig(Path(save_path), dpi=300)
    return fig, ax

def run_plot_series_all(
    x_data, y_data_dict, y_keys,
    axispar_list, savename_list,
    do_plot=True, do_show=True
):
    if do_plot:
        for i, ykey in enumerate(y_keys):
            y = y_data_dict.get(ykey)
            if y is None:
                print(f"Skipping invalid y-value: {ykey}")
                continue

            axis_kwargs = {}
            if i < len(axispar_list):
                axis_kwargs = json.loads(axispar_list[i])

            save_path = None
            if savename_list and i < len(savename_list):
                save_path = savename_list[i]

            label = ykey
            title = f"{ykey} vs Time"

            plot_series(x_data, y, axis_params=axis_kwargs, label=label, title=title, save_path=save_path)

        if do_show:
            plt.show()

"""
