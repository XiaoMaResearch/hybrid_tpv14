//
//  cal_ke.cpp
//  hybrid_fem_bie
//
//  Created by Max on 2/7/18.
//
//

#include "cal_ezz.hpp"
#include "maplocal.hpp"
#include <iostream>
void cal_ezz(MatrixXd coord ,double E, double nu,std::vector<double> &ezz_out, std::vector<double> &exx_out,std::vector<double> &eyy_out,  int n_el, MatrixXi &index_store, double q, VectorXd &u_n)
{
    MatrixXd D=MatrixXd::Zero(3,3);
    // Plane strain
    /*
    D << 1-nu, nu , 0,
    nu,   1-nu, 0,
    0,   0 ,  (1-2*nu)/2;
    D = E/((1+nu)*(1-2*nu))*D;
    */
    // Plane stress
    D << 1.0, nu , 0,
    nu,   1.0, 0,
    0,   0 ,  (1-nu)/2.0;
    D = E/(1-nu*nu)*D;
     VectorXd xi_gp=VectorXd::Zero(2,1);
    VectorXd eta_gp=VectorXd::Zero(2,1);
    VectorXd w_gp=VectorXd::Zero(2,1);
    xi_gp<<-1.0/sqrt(3),1.0/sqrt(3);
    eta_gp<<-1.0/sqrt(3),1.0/sqrt(3);
    w_gp<<1.0,1.0;
    
    for (int i=0;i<n_el;i++){
        VectorXd sigma_temp = VectorXd::Zero(3,1);
        VectorXd strain_temp = VectorXd::Zero(3,1);
    VectorXd u_n_local=VectorXd::Zero(8,1) ;
    ArrayXi index = index_store.col(i);
    maplocal(index,u_n,u_n_local);
    for (int i=0; i<2;i++)
    {
        for (int j=0;j<2;j++)
        {
            double xi = xi_gp(i);
            double eta = eta_gp(j);
            double wgp = w_gp(j);
            // Jacobian Matrix (2x2)
            MatrixXd dNdxi=MatrixXd::Zero(2,4);
            dNdxi << eta-1.0, 1.0-eta, 1.0+eta, -eta-1.0,
                     xi-1.0 , -xi-1.0, xi+1.0, 1.0-xi;
            dNdxi=dNdxi/4.0;
            MatrixXd Jac=MatrixXd::Zero(2,2);
            Jac = dNdxi*coord;
            double detJ = Jac.determinant();
            MatrixXd dNdx = MatrixXd::Zero(2,4);
            dNdx = Jac.inverse()*dNdxi;
            MatrixXd B= MatrixXd::Zero(3,8);
            B<< dNdx(0,0),         0, dNdx(0,1),         0, dNdx(0,2),         0,  dNdx(0,3),         0,
            0, dNdx(1,0),         0, dNdx(1,1),        0 , dNdx(1,2),          0,  dNdx(1,3),
            dNdx(1,0), dNdx(0,0), dNdx(1,1), dNdx(0,1), dNdx(1,2), dNdx(0,2),  dNdx(1,3),  dNdx(0,3);
            
            sigma_temp += 1.0/4.0*D* B*u_n_local;
            strain_temp +=1.0/4.0*B*u_n_local;
         //   ke = ke + B.transpose()*D*B*detJ*wgp;
        }
    }
        ezz_out[i] = -nu/E*(sigma_temp(0)+sigma_temp(1));
        exx_out[i] = strain_temp(0);
        eyy_out[i] = strain_temp(1);
      //  ezz_out[i] = sigma_temp(2);
    }
}
