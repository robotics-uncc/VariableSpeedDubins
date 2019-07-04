# Variable Speed Dubins C++ Library

This code implements a minimum-time planar path planning algorithm for a kinematic car with variable (strictly positive) speed and symmetric turn rate limits. The optimal control problem was studied using both a geometric interpretation of the Minimum Principle (via the hodograph) and using an analytical approach (via the Karush-Kuhn-Tucker conditions). The extremal controls were found to consist of: maximum turn rate and maximum speed turns (denoted as B extremals), straight segments at maximum speed  (S extremals), and cornering turns (C extremals) with maximum turn rate and minimum speed. An extremal turn was found, in general, to be a sequence of three consecutive extremals of the form BCB. A finite sufficient set of candidate optimal controls was derived by analysis of the adjoint differential equations, and by identifying suboptimality conditions geometrically. Candidate paths were found to consist of (at most) a sequence of four turns, or a turn-straight-turn sequence.  A procedure was proposed to solve the path synthesis problem. In particular, a finite dimensional constrained optimization problem (of at most four parameters) was formulated for candidate path types in the finite sufficient set. Numerical solutions to this optimization problem give locally optimal solutions from which the lowest cost solution could be identified. It was found that, in some cases, the variable speed Dubins path is substantially faster paths than a (maximum speed) Dubins path.

<p align="center"> 
<img src="http://arturwolek.com/img/VariableSpeedDubins.png" width="300">
</p>

## Dependencies:
- IPOPT: https://github.com/coin-or/Ipopt
- cddlib: https://inf.ethz.ch/personal/fukudak/cdd_home/
- armadillo: http://arma.sourceforge.net/ (including LAPACK, BLAS, etc.)
- Ubuntu 16.04
- cmake
- (Optional) Octave/MATLAB
  for generating .m files for plotting

## Build instructions:
- Run `./install_deps.sh` to unpack/install tar files in `./external`
- Run `./build.sh` to create a library in `./lib` and executables in `./bin`

## Usage instructions:
- See test programs in `./programs` for examples of how to incorporate into your own code
- `clean.sh` removes all compiled files, build files, and temporary files 

## References:

1. Wolek, A., Cliff, E. M., & Woolsey, C. A. (2016). Time-optimal path planning for a kinematic car with variable speed controls. Journal of Guidance, Control, and Dynamics, 39(10), 2374–2390.
https://doi.org/10.2514/1.G001317

## Contact:
Artur Wolek, wolek@umd.edu

