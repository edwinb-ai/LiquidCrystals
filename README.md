# Liquid Crystals

This is a repository to explore the phase transitions of liquid crystals,
modeled with colloidal rods of length $L$ and spherical cap of diameter
$D$.

## Implementations

- The code `orientation.jl` is an implementation of the paper "The isotropic and nematic liquid crystal phase of colloidal rods" by van Roij (2005). The output of the current version is the equation of state for a particular aspect ratio $L / D$.
    - The code include a preliminary version of the corrections by the Parson-Lee theory to account for finite length rods.

- The code `orientation_gl.jl` is a modification of the paper "The isotropic and nematic liquid crystal phase of colloidal rods" by van Roij (2005). The modification is a complete re-write of the algorithm where we use Gauss-Legendre quadrature for all integration schemes. This improves the efficiency of the code, by reducing considerably the number of points to evaluate. The output of the current version is just the pressure for a given value of the concentration $c$.

## TODO

- Implement the root finding algorithm to obtain coexistence values for the istropic-nematic phase transition.
