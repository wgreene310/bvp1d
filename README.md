# bvp1d
1D Boundary Value Problem solver for Octave and MATLAB

bvp1d solves systems of boundary value differential equations in a single
spatial variable. It is mostly compatible with the MATLAB function bvp4c,
the main difference being that it currently does not include the extra convenience
functions like bvpinit. Most bvp4c examples will work with bvp1d with only small
changes.

The easiest way to get started with bvp1d and Octave is to download one of the
pre-built versions from the Releases area. The bvp1d function is packaged as
a single mex file for Windows, Macintosh, and Linux platforms. Several examples
and basic documentation are included.

This page https://people.sc.fsu.edu/~jburkardt/m_src/bvp4c/bvp4c.html by John Burkhard 
has several examples and links to other resources for bvp4c.

bvp1d relies on the KINSOL library from Sundials 
http://computation.llnl.gov/projects/sundials
for solution of the nonlinear algebraic equations.
In turn, KINSOL relies on the KLU sparse solver from
SuiteSparse http://faculty.cse.tamu.edu/davis/suitesparse.html
for solving sparse, linear equations.
The Eigen C++ matrix class library
http://eigen.tuxfamily.org/index.php?title=Main_Page
is used throughout the code. There is also a small dependency on
the Boost C++ string package, http://www.boost.org/.
So all of these libraries are required to build the software
from scratch. The file, Makefile_octave, was used to build
the versions for the three platforms in the Releases area.
