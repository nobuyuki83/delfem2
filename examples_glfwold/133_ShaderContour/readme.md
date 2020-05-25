# 132_ShaderContour
![](thumbnail.png)

## About

Contour drawing of a 3D oejct using a Non-Photorealistic Rendering (NPR) technique. We first draw normal of the mesh on a texture using the fragment shader. Then the edge is extracted by computing Laplacian of the imge-space normal.


## Technical detail

please take look at the survey for more detail:
```
"Line drawings from 3D models: a tutorial. Foundations and Trends"
Pierre Bénard, Aaron Hertzmann.
in Computer Graphics and Vision, Now Publishers, 2019, 11 (1-2), pp.159.
```

