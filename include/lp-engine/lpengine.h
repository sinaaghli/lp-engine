#ifndef LP_ENGINE_H_
#define LP_ENGINE_H_

#include <iostream>
#include <eigen3/Eigen/Eigen>
#include <lp-engine/lpproblem.h>
#include <cstddef>

class lpengine{
 public:
  lpengine(lpproblem& problem);
  ~lpengine();

  LPPRoblemStatus Solve();

 private:
  bool ReadDicSolution();
  bool IsSemiPositive(Eigen::VectorXd& vec);
  bool IsNegative(Eigen::VectorXd& vec);
  bool CheckDicfeasible(Eigen::VectorXd& x_b_hat);
  int ChooseEnteringVar(Eigen::VectorXd& vec);
  void ChooseLeavingVar(Eigen::VectorXd& delta_var, Eigen::VectorXd& var_hat, int* return_index, double* return_step_length);
  LPPRoblemStatus PrimalSimplexStep();
  LPPRoblemStatus DualSimplexStep();
  void SwapCols(Eigen::MatrixXd& matn, int index_matn, Eigen::MatrixXd& matb, int index_matb);
  void PrintStatusMsg(LPPRoblemStatus& status);
  lpproblem& p_;
};

#endif  // LP_ENGINE_H_
