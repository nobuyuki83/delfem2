import math
import bpy, mathutils


def set_tex_environment(world: bpy.types.World, path: str):
    # shader node tree
    world.use_nodes = True
    bpysnt = world.node_tree # ShaderNodeTree
    bpysn_texenv = bpysnt.nodes.new(type="ShaderNodeTexEnvironment") # ShaderNodeTexEnvironment
    bpysn_texenv.image = bpy.data.images.load(path)
    bpysn_texcrd = bpysnt.nodes.new(type="ShaderNodeTexCoord")
    bpysnt.links.new(bpysn_texcrd.outputs["Generated"], bpysn_texenv.inputs["Vector"])
    bpysnt.links.new(bpysn_texenv.outputs["Color"], bpysnt.nodes["Background"].inputs["Color"])


def set_camera_ydirection():
    # camera (enviroment tex only work for perspective projection)
    bpyobj_cam = bpy.context.scene.camera
    bpyobj_cam.data.sensor_fit = 'HORIZONTAL'
    bpyobj_cam.data.sensor_width = 36.0
    bpyobj_cam.data.sensor_height = 24.0
    bpyobj_cam.data.lens = 20  # focus length
    bpyobj_cam.data.dof.aperture_fstop = 1.4 # f-stop
    bpyobj_cam.rotation_euler = mathutils.Euler((math.radians(90.0), 0.0, 0.0), 'XYZ') # look in the direction of y-axis
    bpyobj_cam.location = (0.0, -6.0, 0.0)
