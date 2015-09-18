# poly3d

This repository reproduces code for `poly3d`, a boundary element method program first released in the M.S. thesis Thomas, 1993.

## poly3d_parallel

A C++ patch that parallelizes parts of `poly3d` is available for beta testing. This version, called `poly3d_parallel`, can currently handle complex fault models with up to 25K elements in around 18 hours using two nodes with dual 2.2GHz, 16-core CPUs. Please contact [rmsare@stanford.edu](rmsare@stanford.edu) if you are interested in testing `poly3d_parallel` as part of an academic project.

## TODO

`poly3d` needs
* Unit tests
* Further parallelization
* Object-oriented framework

Traditionally, `poly3d` has been used in regional studies employing relatively coarse fault meshes. If you use the original or parallel version successfully with a complex model with intersecting faults or a fine mesh, please consider sharing your input file as a test for future releases.

Likewise, feel free to fork this repository to parallelize or rewrite other subroutines.

## Documentation

The primary `poly3d` manual is Thomas, 1993. It covers the details of the boundary value problems behind `poly3d` and explains the input and output file formats.

Michele Cooke (University of Massachussetts, Amherst) maintains a `poly3d` tutorial page that previously included the thesis source code provided here. The original `poly3d` tutorial is available there, as are several examples and other software packages. Other examples are available on the Stanford Structural Geology and Geomechanics group website. 

### Compiling `poly3d`

Both versions of `poly3d` must be compiled with the `-O2` flag. Over-optimization (e.g. `-O3` or `-ffast-math`) can lead to unpredictable numerical errors and large singularities in the resulting displacement fields. Just say no!

### Using `poly3d`

`./poly3d -i <input> -o <output>`

### Contributors 

`poly3d` and a related commercial software package, ibem3d, were developed by researchers and students of David Pollard and Atilla Aydin associated with the Stanford Rock Fracture Project and the Structural Geology and Geomechanics research group. This repository is intended for academic use only, and is solely derived from the publicly available source code printed in Thomas, 1993 with modifications to fix the "shadow effect."

### References

Please cite this thesis and related publications if you use `poly3d` in published work.

Thomas, A. L., 1993, Poly3D: A three-dimensional, polygonal element, displacement discontinuity boundary element computer program with applications to fractures, faults, and cavities in the Earth's crust (Doctoral dissertation, Stanford University).

Questions? Comments? Complaints?
Robert Sare rmsare@stanford.edu
[STGL]()
