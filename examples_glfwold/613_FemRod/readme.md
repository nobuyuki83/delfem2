# 613_FemRod
![](../../docs/imgs/glfwold_613_FemRod.png)


## About

This program demonstrates the minimization of the elastic energy of discrete elastic rod.
Starting from the random configuration, the minimization make the rod return to the initial un-deformed configuration.   

The 3x3 block sparse matrix data structure is defined such that it allows branching rods or coupling with other 3D deformation.    

see also
[622_FemRodStatic](../622_FemRodHairStatic) and [623_FemRodStatic](../622_FemRodHairDynamic).




## Technical detail
This implementation of the discrete elastic rod is based on the following two papers.

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

