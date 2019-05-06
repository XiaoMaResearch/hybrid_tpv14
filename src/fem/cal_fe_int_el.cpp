//
//  cal_fe_int_el.cpp
//  hybrid_fem_bie
//
//  Created by Max on 3/21/18.
//
//

#include "cal_fe_int_el.hpp"

void cal_fe_int_el(double E, double nu, std::vector<Eigen::MatrixXd> &B_mat, double &detJ,VectorXd &u_n_local, VectorXd &fe_int_el,
                   MatrixXd &strain_n, MatrixXd &stress_n,double sxx_initial,
                   double syy_initial, double sxy_initial, double cohes, double blkfric)

{
    for (int i=0;i<4;i++)
    {
        VectorXd strain = B_mat[i]*u_n_local;
       // VectorXd stress = D*strain;
        VectorXd stress = VectorXd::Zero(3,1);
        VectorXd strain_temp = strain_n.col(i);
        VectorXd stress_temp = stress_n.col(i);
        Getplastic(E, nu, strain, stress, stress_temp, strain_temp,sxx_initial,syy_initial,sxy_initial,cohes,blkfric);
        fe_int_el = fe_int_el+B_mat[i].transpose()*stress*detJ;
        //fe_int_el = fe_int_el+B_mat[i].transpose()*D*B_mat[i]*u_n_local*detJ;
        strain_n.col(i) = strain;
        stress_n.col(i) = stress;
    }
}
