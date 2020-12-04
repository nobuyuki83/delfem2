# DelFEM2 C++ Examples using SMPL data





## SMPL Model Preparation

1. Downlaod zip model file```SMPL_python_v.1.0.0.zip``` from https://smpl.is.tue.mpg.de/downloads
2. Put ```basicModel_f_lbs_10_207_0_v1.0.0.pkl``` and ```basicmodel_m_lbs_10_207_0_v1.0.0.pkl``` to ```delfem2/test_inputs```.
3. Install ```chumpy``` with the command ```pip3 install chmpy```
4. Run ```delfem2/test_inputs/smpl_preprocess.py```
5. Then ```smpl_model_f.npz```and ```smpl_model_m.npz``` will apper under ```delfem2/test_inputs/```



## Note 

These demos use OpenGL version 2.1 and GLSL shaer version 1.2 which are depricated in many environment. But still it is convenient to use legacy functions such as glBegin(), glEnd(). We will eventually consider porting these demo into newer OpenGL >= 3.3 in the examples_glfwnew folder.





## Demos

### [00_PoseBone](00_PoseBone)
<img src="00_PoseBone/thumbnail.png" width=200>

### [01_RigTransfer](01_RigTransfer)
<img src="01_RigTransfer/thumbnail.png" width=200>

### [02_Ik](02_Ik)
<img src="02_Ik/thumbnail.png" width=200>

### [03_IkArap](03_IkArap)
<img src="03_IkArap/thumbnail.png" width=200>

### [04_IkImage](04_IkImage)
<img src="04_IkImage/thumbnail.png" width=200>

### [05_Ui](05_Ui)
<img src="05_Ui/thumbnail.png" width=200>

### [10_Cloth](10_Cloth)
<img src="10_Cloth/thumbnail.png" width=200>

### [11_ClothPose](11_ClothPose)
<img src="11_ClothPose/thumbnail.png" width=200>

### [12_IkImageCloth](12_IkImageCloth)

<img src="12_IkImageCloth/thumbnail.png" width=200>

### [13_ClothPoseRig](13_ClothPoseRig)

<img src="13_ClothPoseRig/thumbnail.png" width=200>



### [14_ClothPoseTexture](14_ClothPoseTexture)

<img src="14_ClothPoseTexture/thumbnail.png" width=200>



### [20_PoseBoneBlendshape](20_PoseBoneBlendshape)

<img src="20_PoseBoneBlendshape/thumbnail.png" width=200>




