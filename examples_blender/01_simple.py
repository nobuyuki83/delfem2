import bpy
import os, sys

def main():	
    # render
    bpy.context.scene.render.resolution_percentage = 50
    bpy.context.scene.cycles.samples = 60
    bpy.context.scene.render.engine = 'CYCLES'
    bpy.ops.render.render()
    bpy.data.images['Render Result'].save_render(filepath = os.path.dirname(__file__)+'/out/01_out.png')

if __name__ == "__main__":    
	main()   
