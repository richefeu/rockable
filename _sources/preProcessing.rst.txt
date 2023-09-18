.. _prePro:

Pre-processing commands
=======================

* ``stickVerticesInClusters`` (*double*) **Epsilon** 
  This command will add glued interfaces between bodies having the same cluster identifier. 
  Only bonds between vertices (spheres) are added when the distance is less than **Epsilon**.

* ``stickClusters`` (*double*) **Epsilon**  
  This command will add glued interfaces between bodies having different cluster identifier. 
  Bonds are added when the distance is less than **Epsilon**.
  
* ``setAllVelocities`` (*vec3r*) **velocity**
  Set the velocity vector of all particles (that are not driven) to the prescribed value **velocity**.




