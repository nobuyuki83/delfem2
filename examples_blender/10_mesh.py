import os, sys, math
import bpy, bmesh, mathutils
import PyDelFEM2 as dfm2

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import myutil # my functions

def main():
	# initialize
    bpy.data.objects.remove(bpy.data.objects["Cube"])   

    dfm2_msh = dfm2.Mesh()
    dfm2_msh.read("../test_inputs/bunny_2k.ply")
    dfm2_msh.normalize(3.0)

    bpy_msh = bpy.data.meshes.new(name="cubemesh")
    bpy_msh.from_pydata(dfm2_msh.np_pos.tolist(), [], dfm2_msh.np_elm.tolist())
    bpy_msh.update()
    bpyobj_msh = bpy.data.objects.new(name="cube", object_data=bpy_msh)
    bpy.context.scene.collection.objects.link(bpyobj_msh)

    myutil.set_tex_environment(
            bpy.context.scene.world,
            '../test_inputs/epping_forest_01_1k.hdr')

    myutil.set_camera_ydirection()

    # render
    bpy.context.scene.render.resolution_percentage = 50
    bpy.context.scene.cycles.samples = 60
    bpy.context.scene.render.engine = 'CYCLES'
    bpy.ops.render.render()
    bpy.data.images['Render Result'].save_render(filepath = os.path.dirname(__file__)+'/out/10_out.png')

if __name__ == "__main__":    
	main()   
