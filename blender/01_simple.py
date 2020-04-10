import bpy
import os
import sys

if __name__ == "__main__":
    
    output_file_path = "out/01_"
    resolution_percentage = 50
    num_samples = 64

    camera_object = bpy.data.objects["Camera"]
    bpy.context.scene.render.resolution_percentage = resolution_percentage
    bpy.context.scene.render.filepath = output_file_path

    bpy.context.scene.render.image_settings.file_format = 'PNG'
    bpy.context.scene.render.engine = 'CYCLES'

    bpy.context.scene.cycles.samples = num_samples
