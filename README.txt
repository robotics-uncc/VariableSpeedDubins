Variable Speed Dubins C++ Library
Last Modified: 3-July-2019
Artur Wolek

Dependencies:
- IPOPT: https://github.com/coin-or/Ipopt
- cddlib: https://inf.ethz.ch/personal/fukudak/cdd_home/
- armadillo: http://arma.sourceforge.net/ (including LAPACK, BLAS, etc.)
- Ubuntu 16.04
- cmake
- (Optional) Octave/MATLAB
  for generating .m files for plotting

Build instructions:
- Run ./install_deps.sh to unpack/install tar files in ./external
- Run ./build.sh to create
  a) a library ./lib/libVarSpeedDubins.a
  b) executables ./bin

Usage instructions:
- Use the library in your own code
- See test programs for examples of how to incorporate into your own code
- clean.sh removes all compiled files, build files, and temporary files 

References:

[1] Wolek, A., Cliff, E. M., & Woolsey, C. A. (2016). Time-optimal path planning for a kinematic car with variable speed controls. Journal of Guidance, Control, and Dynamics, 39(10), 2374–2390.
https://doi.org/10.2514/1.G001317

