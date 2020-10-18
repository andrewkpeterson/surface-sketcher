<p align="center">
<img src="https://github.com/andrewkpeterson/surface-sketcher/blob/master/examples/results.png"
</p>
  
# surface-sketcher

This is an independent sketch-based modeling project. In this project, I implemented the main surface-solving algorithm from [BendSketch](https://haopan.github.io/bendsketch.html) by Li et al. This program takes in a set of boundary, contour, and bending (i.e. principal curvature) strokes and creates a surface described by the strokes. 

The program is shown below. The user can pick which kind of stroke they would like to create, and then draw the stroke using a mouse or tablet in the window. Black lines represent boundary strokes or contour strokes, red lines represent convex bending strokes, and blue lines represent concave bending strokes. The "Start Height", "End Height", and "Radius" controls are for defining the height curves of contour and boundary strokes.
 
<p align="center">
<img src="https://github.com/andrewkpeterson/surface-sketcher/blob/master/examples/surface-sketcher.png" height="500">
</p>

Note that this is not a full implementation of the BendSketch program. Only the main steps of the surface-solving algorithm presented in the BendSketch paper are implemented here. In particular, automatic concave/convex stroke parsing, flat strokes, ridges, valleys, sharp features, Laplacian deformation, and multi-view modeling are not implemented.

### Overview of Pipeline

Boundary strokes describe the boundaries of the surface being created. Contour lines describe the boundaries of the surface and additionally constrain the surface normals so that they are as orthogonal as possible to the domain boundary in the drawing plane. Convex/concave bending strokes describe convex/concave principal curvature lines inside the domain.

The final surface created by the program is a triangulation of the two-dimensional domain described by the boundary and contour strokes. The z-coordinates of the vertices in this two-dimensional domain are optimized so that the surface has the principal curvatures described by the bending strokes, and matches the height curves of the domain boundary.

#### Domain Triangulation

To triangulate the domain defined by the boundary and contour strokes, I used the Delaunay triangulation implementation from CGAL.

#### Creating the Principal Curvature Direction Field

Next, the principal curvature directions are calculated at every triangle. The principal curvature direction field is described by four directions at each triangle. The curvature directions are paired, so that the two directions in a pair are related by a rotation of pi radians. Note that the principal curvature field is not orthogonal, as the orthongonality of principal curvature directions is not preserved through projection. 

The principal curvature directions are calculated by optimizing a PolyVector representation. This involves optimizing the coefficients of a complex polynomial defined at each triangle whose roots are the desired principal curvature directions at that triangle. We optimize these coefficients so that they are smooth across the mesh, encouraged to be orthogonal, and obey the constraints defined by the bending lines (i.e. the direction defined by the bending line is a root of the polynomial of the triangle that the bending line intersects.) The energies that are optimized to find the PolyVector representation are all quadratic functions of the coefficients, so optimizing the PolyVector representation reduces to solving a system of linear equations.

Below is the principal curvature direction field of the teapot example.

<p align="center">
<img src="https://github.com/andrewkpeterson/surface-sketcher/blob/master/examples/curvature_directions.png" width="500" height="500">
</p>

Note that the BendSketch paper optimizes the BendField energy of the principal curvature direction field to further refine the curvature directions. This step is implemented in this project, but it is not used in any of the examples in this README.

#### Jointly Optimizing the Surface and Principal Curvature Magnitudes

Before the we can jointly optimize the surface (i.e. the vertices' z-coordinates) and the principal curvature magnitudes defined at each triangle, we must intiailize the surface. To do this, we first assign the same arbitrary magnitude to each of the bending strokes. Then, the curvature magnitudes are smoothed out through the mesh under the constraints of the bending strokes. Next, the surface is optimized to have the known curvature values and directions at each triangle.

Now, the surface and principal curvature magnitudes can be jointly optimized. This process is the same as the initialization step, except the curvature magnitudes along the bending strokes of the current surface are estimated, rather than given arbitrary values. Estimating the curvature magnitudes involves fitting a circle to the surface along the bending line, and then taking the reciprocal of the circle's radius. Fitting the circle is a nonlinear least squares problem, so I used the Ceres Solver for this step. Smoothing the curvature magnitudes and optimizing the surface to have the known curvature magnitudes and directions proceeds in the same way as above. This process continues until the surface converges, which takes 5-8 iterations.

The image below shows the curvature directions and mangitudes of the teapot surface in the middle of optimization. Red arrows represent convex principal curvature, and blue arrows represent concave principal curvature.

<p align="center">
<img src="https://github.com/andrewkpeterson/surface-sketcher/blob/master/examples/curvature_magnitudes.png" height="500">
</p>

Finally, we have the solved surface.

<p align="center">
<img src="https://github.com/andrewkpeterson/surface-sketcher/blob/master/examples/teapot_image.png" height="500">
</p>

### Dependecies

* Eigen 3.3.7
* CGAL 5.0.2
* Ceres-Solver 1.14.0

### References

These are the books and papers that I found helpful when working on this project.

* BendSketch: Modeling Freeform Surfaces Through 2D Sketching by Li et al. (This was the main motivation for this project.)
* Polygon Mesh Processing by Botsch et al. (This book was useful for understanding how a gradient of a function defined at each vertex in a triangle mesh is calculated.)
* Differential Geometry of Curves and Surfaces by Tapp (This book is a great introduction to differential geometry.)
* Designing N-PolyVector Fields with Complex Polynomials by Diamanti et al. (This paper helped me understand the PolyVector representation, which is used to find the curvature direction field.)
* BendFields: Regularized Curvature Fields from Rough Concept Sketches by Iarussi et al. (This paper was useful for understanding the BendField energy used in the BendSketch paper.)
* Rotational Symmetry Field Design on Surfaces by Palacios et al. (This paper helped me learn more about N-RoSy's, or N-way rotational symmetries, which were generalized by the PolyVector representation.)

