# CGoGN_3

CGoGN is a C++ library that provides:
- an implementation of several mesh data structures (with a versatile cells attributes mechanism):
  - combinatorial maps (1D, 2D, 3D)
  - graphs (1D)
  - incidence graphs (1D, 2D)
- several geometric modeling and geometry processing algorithms (mostly generic w.r.t. the underlying data structure):
  - surface remeshing, subdivision, decimation
  - surface geometry filtering, curvature computation
  - surface skeleton extraction
  - surface deformation
  - convex hull
  - hex mesh generation
- an ImGui interface with modules and signals
- openGL rendering for points, lines, surface, volume meshes

# Building CGoGN

- create a build directory and go into it
- run `cmake <cgogn_path>`
- you can specify build type (and other options) by modifying CMAKE_BUILD_TYPE either by using cmake-gui or `ccmake`, or by specifying -DCMAKE_BUILD_TYPE="" when running cmake

### LINUX and MacOS

`make -jN` (N being the number of threads) in the build directory

### Windows

VS 2013 or better required

# Contribution HowTo

- Fork this GitHub repository
- Clone your GitHub fork on your working machine
- Checkout the "develop" branch
- From here, create your own "[user]" branch to commit your work in
- Push your branch with its commits in your GitHub for
- Create a Pull Request from your "[user]" branch to the "develop" branch of this repository
- If you want to update your repository with commits from other contributors, with your "[user]" branch as current branch, you can pull the "develop" branch of this repository. This will fetch new commits and merge them in your branch.
