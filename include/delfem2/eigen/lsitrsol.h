

#include <Eigen/Core>


// {y_vec} = \beta * {y_vec} + \alpha * [mat] * {p_vec}
void AddMatVec(
    Eigen::VectorXd& y_vec,
    double beta,
    double alpha,
    const Eigen::MatrixXd& mat,
    const Eigen::VectorXd& x_vec)
{
  y_vec = alpha*mat*x_vec + beta*y_vec;
}

void AddScaledVec(
    Eigen::VectorXd& y,
    double alpha,
    const Eigen::VectorXd& x)
{
  y += alpha*x;
  //const std::size_t n = x.n;
  //for (unsigned int i = 0; i < n; i++) { y.p[i] += alpha * x.p[i]; }
}

void ScaleAndAddVec(
    Eigen::VectorXd& y,
    double beta,
    const Eigen::VectorXd& x)
{
  y = beta*y + x;
//  const std::size_t n = x.n;
//  for (unsigned int i = 0; i < n; i++) { y.p[i] = beta * y.p[i] + x.p[i]; }
}


#include "delfem2/lsitrsol.h"
