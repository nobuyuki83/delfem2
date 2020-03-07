import numpy as np
import pickle
import json
import sys


def export_numpy(path_in, fname_out):

  with open(path_in, 'rb') as f:
    src_data = pickle.load(f, encoding="latin1")

  vertices_template = np.array(src_data['v_template'])
  face_indices = np.array(src_data['f'] + 1)  # starts from 1
  weights = np.array(src_data['weights'])
  shape_blend_shapes = np.array(src_data['shapedirs'])
  pose_blend_shapes = np.array(src_data['posedirs'])
  joint_regressor = np.array(src_data['J_regressor'].toarray())
  kinematic_tree = np.array(src_data['kintree_table'])    

  model_data_np = {
    'vertices_template': vertices_template,
    'face_indices': face_indices,
    'weights': weights,
    'shape_blend_shapes': shape_blend_shapes,
    'pose_blend_shapes': pose_blend_shapes,
    'joint_regressor': joint_regressor,
    'kinematic_tree': kinematic_tree
  }

  np.savez(fname_out, **model_data_np)


if __name__ == '__main__':
  export_numpy('basicModel_f_lbs_10_207_0_v1.0.0.pkl', "smpl_model_f.npz")
  export_numpy('basicModel_m_lbs_10_207_0_v1.0.0.pkl', "smpl_model_m.npz")  
