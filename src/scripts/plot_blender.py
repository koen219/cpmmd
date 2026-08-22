import bpy
import numpy as np
import click
from pathlib import Path

import bpy
import mathutils
import numpy as np


def clear_scene():
    """Delete all objects in the Blender scene."""
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()


def create_camera():
    """Create and position a camera in the scene."""
    cam = bpy.data.objects.new("Camera", bpy.data.cameras.new("Camera"))
    bpy.context.collection.objects.link(cam)
    bpy.context.scene.camera = cam
    cam.location = (150, 150, 150)

    direction = mathutils.Vector((-200, -200, -200))
    rot_quat = direction.to_track_quat("-Z", "Y")
    # cam.angle = (90.0/180.0) * np.pi

    cam.data.lens = 20
    # assume we're using euler rotation
    cam.rotation_euler = rot_quat.to_euler()


# def create_camera():
#     # create the first camera object
#     cam1 = bpy.data.cameras.new("Camera 1")
#     cam1.lens = 18
#
#     cam_obj1 = bpy.data.objects.new("Camera 1", cam1)
#     cam_obj1.location = (100, 100, 100)
#     cam_obj1.rotation_euler = (0.0, 0.0, 0.0)
#     bpy.context.scene.collection.objects.link(cam_obj1)
#     bpy.context.scene.camera = cam_obj1


#     """Create and position a camera in the scene, focusing on (100, 100, 100) from a distance of 100."""
#     focus_point = np.array([100, 100, 100])
#     distance = 100
#
#     # Choose a good camera position: directly above and offset
#     camera_position = focus_point + np.array(
#         [0, -distance, distance]
#     )  # Adjust as needed
#
#     cam = bpy.data.objects.new("Camera", bpy.data.cameras.new("Camera"))
#     bpy.context.collection.objects.link(cam)
#     bpy.context.scene.camera = cam
#     cam.location = camera_position
#
#     cam.rotation_euler = (np.radians(60), 0, 0)  # rot_quat.to_euler()


def direction_to_rotation(direction):
    """Convert a direction vector to a rotation quaternion so the camera looks at a target."""
    direction = direction / np.linalg.norm(direction)  # Normalize the vector
    z_axis = np.array([0, 0, 1])  # Default up direction
    rotation_axis = np.cross(z_axis, direction)
    rotation_angle = np.arccos(np.dot(z_axis, direction))

    if np.linalg.norm(rotation_axis) < 1e-6:
        return mathutils.Quaternion((1, 0, 0, 0))  # No rotation needed

    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
    return mathutils.Quaternion(rotation_axis.tolist() + [rotation_angle])


def create_material(color, emission_strength=5.0):
    """Create an emissive material to enhance visibility."""
    mat = bpy.data.materials.new(name="EmissiveMaterial")
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links

    # Clear default nodes
    for node in nodes:
        nodes.remove(node)

    # Create emission material
    output_node = nodes.new(type="ShaderNodeOutputMaterial")
    emission_node = nodes.new(type="ShaderNodeEmission")

    emission_node.inputs[0].default_value = color  # RGBA color
    emission_node.inputs[1].default_value = emission_strength
    links.new(emission_node.outputs[0], output_node.inputs[0])

    return mat


def create_lines(vertices, edges, color=(1, 1, 1, 1)):
    """Create a single mesh object containing multiple lines."""
    mesh = bpy.data.meshes.new(name="LinesMesh")
    obj = bpy.data.objects.new("LinesObject", mesh)
    bpy.context.collection.objects.link(obj)

    mesh.from_pydata(vertices, edges, [])
    mesh.update()

    mat = create_material(color, emission_strength=5.0)
    obj.data.materials.append(mat)


def create_large_curve(vertices, edges, thickness=0.02, color=(1, 1, 1, 1)):
    """Create a single curve object contaifning all lines efficiently."""

    # Create one Curve object
    curve_data = bpy.data.curves.new(name="LargeCurve", type="CURVE")
    curve_data.dimensions = "3D"

    # Convert edges into multiple splines inside the same curve
    for edge in edges:
        spline = curve_data.splines.new(type="POLY")
        spline.points.add(1)
        spline.points[0].co = (*vertices[edge[0]], 1)
        spline.points[1].co = (*vertices[edge[1]], 1)

    # Set bevel depth to make lines visible
    curve_data.bevel_depth = thickness  # Controls line thickness

    # Create the single object to store all lines
    curve_obj = bpy.data.objects.new("LargeLinesObject", curve_data)
    bpy.context.collection.objects.link(curve_obj)

    # Optional: Assign a material with emission
    mat = bpy.data.materials.new(name="LineMaterial")
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        mat.node_tree.nodes.remove(bsdf)
    emission = mat.node_tree.nodes.new("ShaderNodeEmission")
    emission.inputs["Color"].default_value = color  # Line color
    emission.inputs["Strength"].default_value = 5.0  # Glow effect
    mat.node_tree.links.new(
        emission.outputs["Emission"],
        mat.node_tree.nodes["Material Output"].inputs["Surface"],
    )

    curve_obj.data.materials.append(mat)


def create_markers(positions, color=(1, 1, 0, 1), size=0.1):
    """Create a single mesh object containing all markers as vertices."""
    mesh = bpy.data.meshes.new(name="MarkersMesh")
    obj = bpy.data.objects.new("MarkersObject", mesh)
    bpy.context.collection.objects.link(obj)

    vertices = [tuple(p) for p in positions]
    edges = []  # No edges, only points
    mesh.from_pydata(vertices, edges, [])
    mesh.update()

    mat = create_material(color, emission_strength=5.0)
    obj.data.materials.append(mat)

    obj.show_instancer_for_viewport = True
    obj.show_instancer_for_render = True


def enable_viewport_shading():
    """Ensure Blender is in rendered shading mode."""
    for area in bpy.context.screen.areas:
        if area.type == "VIEW_3D":
            for space in area.spaces:
                if space.type == "VIEW_3D":
                    space.shading.type = "MATERIAL"


def visualize_data(cpm, ecm):
    """Visualize ECM and CPM data in Blender efficiently."""
    clear_scene()
    create_camera()
    enable_viewport_shading()

    pos = ecm["positions"]
    polymer = ecm["group"][ecm["bond_types"] == 0]
    crosslink = ecm["group"][ecm["bond_types"] > 0]

    vertices = [tuple(p) for p in pos]
    polymer_edges = [(polymer[i, 0], polymer[i, 1]) for i in range(len(polymer))]
    crosslink_edges = [
        (crosslink[i, 0], crosslink[i, 1]) for i in range(len(crosslink))
    ]

    # create_lines(vertices, polymer_edges, color=(0.3, 0.3, 0.3, 1))
    create_large_curve(vertices, polymer_edges, color=(0.3, 0.3, 0.3, 1))
    # create_lines(vertices, crosslink_edges, color=(0.0, 0.7, 0.0, 1))
    create_large_curve(vertices, crosslink_edges, color=(0.0, 0.7, 0.0, 1))
    boundary_pos = pos[ecm["particle_types"] == 1]
    create_markers(boundary_pos, color=(1, 1, 0, 1), size=0.15)

    create_pixels(cpm, box_size=1.0, color=(1, 0, 0, 1))

    adhesion_pos = pos[ecm["particle_types"] == 2]
    create_markers(adhesion_pos, color=(0, 0, 1, 1), size=0.2)

    print("Visualization complete.")


import bmesh


def create_pixels(grid, box_size=1.0, color=(1, 0, 0, 0.8)):
    """Create a single mesh object containing all red boxes where grid == 1."""

    mesh = bpy.data.meshes.new(name="PixelsMesh")
    obj = bpy.data.objects.new("PixelsObject", mesh)
    bpy.context.collection.objects.link(obj)

    bm = bmesh.new()

    # Iterate over grid and add cubes where grid[x, y, z] == 1
    dx, dy, dz = grid.shape
    for x in range(dx):
        for y in range(dy):
            for z in range(dz):
                if grid[x, y, z] == 1:
                    # Create a cube at (x, y, z)
                    bmesh.ops.create_cube(
                        bm,
                        size=box_size,
                        matrix=mathutils.Matrix.Translation((x, y, z)),
                    )
                # bmesh.ops.create_cube(bm, size=box_size, matrix=np.eye(4) @ np.array([
                #     [1, 0, 0, x],
                #     [0, 1, 0, y],
                #     [0, 0, 1, z],
                #     [0, 0, 0, 1]
                # ]))

    # Convert BMesh to Blender mesh
    bm.to_mesh(mesh)
    bm.free()

    # Create material
    mat = bpy.data.materials.new(name="PixelMaterial")
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf:
        bsdf.inputs["Base Color"].default_value = color  # Red color
    obj.data.materials.append(mat)


def render_image(output_path):
    """Render and save the image."""
    bpy.context.scene.render.filepath = output_path
    bpy.ops.render.render(write_still=True)
    print(f"Image saved to {output_path}")


def load_data(filename, time):
    if isinstance(time, int):
        time = str(time).zfill(7)

    state = np.load(filename / ("state_" + time + ".npz"))
    cpm = state["cpm"].reshape((200, 200, 200))
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


@click.group
def main():
    pass


# Example usage (assuming data is loaded as numpy arrays)
@main.command
@click.argument("folder")
def video(folder):
    folder = Path(folder)
    create_camera()
    for k, data_file in enumerate(Path(folder).glob("*npz")):
        print(data_file)
        time = data_file.stem.removesuffix(".npz").split("_")[1]
        # Load data
        cpm, ecm = load_data(data_file.parents[0], time)
        time = int(data_file.stem.split("_")[1].split(".")[0])
        # Plot the 3D grid
        visualize_data(cpm, ecm)
        render_image(f"{folder}/3d_simulation_{str(time).zfill(7)}.png")


@main.command
def test():
    # Generate random 3D positions for ECM
    num_points = 100
    positions = np.random.rand(num_points, 3) * 200  # Scale positions for visibility

    # Create random connectivity groups for polymer and crosslinks
    num_edges = 50
    group_polymer = np.random.randint(0, num_points, size=(num_edges, 2))
    group_crosslink = np.random.randint(0, num_points, size=(num_edges, 2))

    # Assign random bond types (0 = polymer, 1+ = crosslink)
    bond_types = np.zeros(num_edges, dtype=int)
    bond_types[: num_edges // 2] = 1  # Half are crosslinks

    # Assign particle types (1 = boundary, 2 = adhesion)
    particle_types = np.zeros(num_points, dtype=int)
    particle_types[np.random.choice(num_points, num_points // 10, replace=False)] = (
        1  # Some boundaries
    )
    particle_types[np.random.choice(num_points, num_points // 10, replace=False)] = (
        2  # Some adhesions
    )

    # Create CPM grid (for now, an empty array)
    cpm = np.zeros((10, 10, 10), dtype=int)

    # Construct ECM dictionary
    ecm = {
        "positions": positions,
        "group": np.vstack([group_polymer, group_crosslink]),
        "bond_types": np.concatenate(
            [np.zeros(len(group_polymer)), np.ones(len(group_crosslink))]
        ),
        "particle_types": particle_types,
    }

    # Run visualization with test data
    visualize_data(cpm, ecm)
    render_image(f"test.png")


if __name__ == "__main__":
    main()
