//
//  cal_strain_energy_density_hpp.hpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//

#ifndef cal_strain_energy_density_hpp
#define cal_strain_energy_density_hpp

#include <stdio.h>
#include <Eigen/Eigen>

using namespace Eigen;
void cal_strain_energy_density(std::vector<MatrixXd> &coord ,double E, double nu,std::vector<double>&se_out, std::vector<double> &exx_out,std::vector<double> &eyy_out,  int n_el, MatrixXi &index_store, double q, VectorXd &u_n);
#endif /* cal_strain_energy_density */
