# Liquid Crystals

This is a repository to explore the phase transitions of liquid crystals,
modeled with colloidal rods of length $L$ and spherical cap of diameter
$D$.

## Implementations

- The code `orientation.jl` is an implementation of the paper "The isotropic and nematic liquid crystal phase of colloidal rods" by van Roij (2005). The output of the current version is the equation of state for a particular aspect ratio $L / D$.
    - The code include a preliminary version of the corrections by the Parson-Lee theory to account for finite length rods.

## TODO

- Reimplement the current version of the code to use Gauss-Legendre quadratures instead of the trapezoidal rule. Since the orientation distribution is smooth, this method of integration should converge exponentially.
