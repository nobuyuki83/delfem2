# 622_FemRodHair
![](../../docs/imgs/glfwold_622_FemRodHairStatic.png)



## About
This demo shows the minimization of the deformation energy of the discrete elastic rod. 
Starting from random configuration, 
the minimization of the deformation energy leads to the initial un-deformed configuration.

[613_FemRod](../613_FemRod/readme.md) solves the similar problem, 
but the formulation of the system matrix is specialized for the discretization of hairs.



## Technical detail
The implementation of the rod is based on following two papers:

```    
"Discrete Viscous Threads" 
Miklós Bergou, Basile Audoly, Etienne Vouga, Max Wardetzky, Eitan Grinspun
ACM Transactions on Graphics (SIGGRAPH) 2010
http://www.cs.columbia.edu/cg/threads/
```
   
```   
"Efficient Yarn-based Cloth with Adaptive Contact Linearization"
Jonathan M. Kaldor, Doug L. James, Steve Marschner
ACM SIGGRAPH 2010
https://research.cs.cornell.edu/YarnCloth/   
```
 
