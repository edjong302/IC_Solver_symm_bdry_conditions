/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYPHIFUNCTION_HPP
#define MYPHIFUNCTION_HPP

Real my_phi_function(RealVect loc, Real amplitude, Real wavelength, RealVect L)
{
    // Real x = loc[0];
    // Real y = loc[1];
    // Real z = loc[2];
    // Real rr = pow(pow(x, 2.) + pow(y, 2.) + pow(z, 2.), .5);
    // gaussian - use with Dirichlet BCs
    Real r2 = D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]);
    return amplitude * .5 * (1 + tanh((sqrt(r2) - 30) / 7.5));
}

#endif /* MYPHIFUNCTION_HPP */
