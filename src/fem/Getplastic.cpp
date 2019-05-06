//
//  Getplastic.cpp
//  hybrid_fem_bie
//
//  Created by Max on 3/20/18.
//
//

#include "Getplastic.hpp"

void Getplastic(double E, double nu , VectorXd &strain, VectorXd &stress , VectorXd &stress_n, VectorXd &strain_n,double sxx_initial,
                double syy_initial, double sxy_initial, double cohes, double blkfric)
{
//    double cohes = 0.0;
//    double blkfric = 0.6;
    double angfric = std::atan(blkfric);
    
    VectorXd strain_inc = strain-strain_n;
    double exx = strain_inc(0);
    double eyy = strain_inc(1);
    double exy = strain_inc(2);
    
    double sxx = stress_n(0);
    double syy = stress_n(1);
    double sxy = stress_n(2);
    
//    double density = 2670;
//    double vp = 6000.0 ;
//    double vs = 3464.0;
//    double mu = density*(vs*vs);
//    double lambda = density*(vp*vp-2.0*(vs*vs));
    double lambda = E*nu/((1+nu)*(1-2*nu));
    double mu = E/(2.0*(1+nu));
    
    double etrace = exx+eyy;
    sxx = sxx + lambda*etrace+2.0*mu*exx;
    syy = syy + lambda*etrace+2.0*mu*eyy;
    sxy = sxy + 2.0*mu*exy/2.0;

    
//    
//    double sxx_initial = -122.54e6;
//    double syy_initial = -50.0e6;
//    double sxy_initial = 19.285e6;
//
    sxx = sxx + sxx_initial;
    syy = syy + syy_initial;
    sxy = sxy + sxy_initial;
//    
    double szz = 0.5*(sxx+syy);
    double sm = (sxx+syy+szz)/3.0;
    // Find stress divatoric components
    double sdxx = sxx - sm;
    double sdyy = syy - sm;
    double sdzz = szz - sm;
    double sdxy = sxy;
    // Secondar invariant of stress
    double secinv = 0.5*(sdxx*sdxx+sdyy*sdyy+sdzz*sdzz)+sdxy*sdxy;
    // scalar measure of shear stress
    double tau = std::sqrt(secinv);
    double taulim = cohes*std::cos(angfric)-(sm)*std::sin(angfric);
    taulim = std::max(0.0, taulim);
//    
    if (tau>taulim)
    {
        double yldfac = taulim/tau;
        sxx = sdxx*yldfac + sm;
        syy = sdyy*yldfac + sm;
        sxy = sdxy*yldfac;
        
        sxx = sxx-sxx_initial;
        syy = syy-syy_initial;
        sxy = sxy-sxy_initial;
    }
    else
    {
        sxx = sxx-sxx_initial;
        syy = syy-syy_initial;
        sxy = sxy-sxy_initial;
    }
    
    stress(0) = sxx;
    stress(1) = syy;
    stress(2) = sxy;
}
