#ifndef LPDICTIONARY_H_
#define LPDICTIONARY_H_

#include <Eigen/Eigen>

class dictionary{
 public:
  dictionary();
  ~dictionary();

  Eigen::MatrixXd a_basic_;
  Eigen::MatrixXd a_nonbasic_;
  Eigen::VectorXd x_basic_;
  Eigen::VectorXd x_nonbasic_;
  Eigen::VectorXd c_basic_;
  Eigen::VectorXd c_nonbasic_;
  Eigen::VectorXd b_;
};

#endif  //LPDICTIONARY_H_
