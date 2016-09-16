# bvp1d
1D Boundary Value Problem solver for Octave and MATLAB

bvp1d solves systems of boundary value differential equations in a single
spatial variable. It is mostly compatible with the MATLAB function bvp4c,
one of the main differences being that it currently does not include the extra convenience functions like bvpinit. It also does not support multipoint
boundary value problems or equations with complex numbers.
Nevertheless, many bvp4c examples will work with bvp1d with only small
changes.

The theory for bvp4c (and bvp1d) is presented in the paper
by [Kierzenka and Shampine](http://dl.acm.org/citation.cfm?id=502801).

The easiest way to get started with bvp1d and Octave is to download one of the
pre-built versions from the Releases area. The bvp1d function is packaged as
a single mex file for Windows, Macintosh, and Linux platforms. Several examples
and basic documentation are included.

The page [bvp4c](https://people.sc.fsu.edu/~jburkardt/m_src/bvp4c/bvp4c.html) by John Burkhard 
has several examples and links to other resources for bvp4c.

bvp1d relies on the KINSOL library from [Sundials] 
(http://computation.llnl.gov/projects/sundials)
for solution of the nonlinear algebraic equations.
In turn, KINSOL relies on the KLU sparse solver from
[SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html)
for solving sparse, linear equations.
The [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) C++ matrix class library,
is used throughout the code. There is also a small dependency on
the [Boost](http://www.boost.org/) C++ string package.
So all of these libraries are required to build the software
from scratch. The file, Makefile_octave, was used to build
the versions for the three platforms in the Releases area.
