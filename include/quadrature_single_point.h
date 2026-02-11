#include <Eigen/Dense>
#include <EigenTypes.h>
#include <iostream>

//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template<typename Ret, typename Integrand_Func>
inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q, 
                                               Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                               Integrand_Func integrand) {
    // integrand(integrated, q, element, Eigen::Ref<const Eigen::Vector3d> X);
    Eigen::Vector3d q0;
    q0 << q(3 * element(0)), q(3 * element(0) + 1), q(3 * element(0) + 2);
    // if (element(0) == 223) {
    //     std::cout << q0 << std::endl;
    // }
    // std::cout << "quadrature single point q0: " << q0 << std::endl;
    // integrated = volume * integrand(integrated, q, element, q0);
    integrand(integrated, q, element, q0);
    integrated = volume * integrated; //???????
    // std::cout << "quadrature single point result: " << integrated << std::endl;
}

