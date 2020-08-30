<p align="center">
<img src="https://github.com/andrewkpeterson/surface-sketcher/blob/master/examples/results.png"
</p>
  
# surface-sketcher

This is an implementation of the main surface-solving algorithm from [BendSketch](https://haopan.github.io/bendsketch.html) from SIGGRAPH 2017 by Li et al. This program takes in a set of boundary, contour, and bending (i.e. principal curvature) strokes and creates a surface described by the strokes. 

<p align="center">
<img src="https://github.com/andrewkpeterson/surface-sketcher/blob/master/examples/surface_sketcher.png" height="500">
</p>

Note that this is not a full implementation of the BendSketch program. Only the main steps of the algorithm presented in the BendSketch paper are implemented here. In particular, flat strokes, ridges, valleys, sharp features, and Laplacian deformation are not implemented.

### Overview of Pipeline

Boundary strokes describe the boundaries of the surface being created. Contour lines describe the boundaries of the surface and additionally constrain the surface normals so that they are as orthogonal as possible to the domain bounadry in the drawing plane. Convex/concave bending strokes describe convex/concave principal curvature lines inside the domain.

<p align="center">
<img src="https://github.com/andrewkpeterson/surface-sketcher/blob/master/examples/curvature_directions.png" width="500" height="500">
</p>

### References

These are the books and papers that I found helpful when working on this project.

