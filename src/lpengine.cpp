#include <lp-engine/lpengine.h>

lpengine::lpengine(lpproblem& problem) : p_(problem) {}

lpengine::~lpengine() {}

bool lpengine::CheckDicfeasible(Eigen::VectorXd& x_b_hat) {
  for (int ii = 0; ii < x_b_hat.rows(); ii++) {
    if (!(x_b_hat[ii] >= 0)) {
      return false;
    }
  }
  return true;
}

int lpengine::ChooseEnteringVar(Eigen::VectorXd& vec) {
  // choose smallest index (Bland's rule)
  for (int ii = 0; ii < vec.rows(); ii++) {
    if (vec[ii] < 0) {
      return ii;
    }
  }
//  double var = vec[0];
//  int ind = 0;
//  for (int ii=1; ii<vec.rows(); ii++) {
//    if (vec[ii] < var) {
//      var = vec[ii];
//      ind = ii;
//    }
//  }
//  if (vec[ind] < 0) {
//    return ind;
//  }
  std::cout << "-- Error Didn't find entring variable. Halting ..."
            << std::endl;
  // TODO: this should never happen but figure out how to handle this situation.
  exit(EXIT_FAILURE);
}

void lpengine::ChooseLeavingVar(Eigen::VectorXd& delta_var,
                                Eigen::VectorXd& var_hat,
                                int* return_index,
                                double* return_step_length) {
  Eigen::VectorXd tmp(var_hat.rows());
  // choose smallest index (Bland's rule)
  for (int ii = 0; ii < var_hat.rows(); ii++) {
    if (var_hat[ii] == 0) {
      tmp[ii] = 0;
      *return_index = ii;
      return_step_length = 0;
      return;
    } else {
      tmp[ii] = delta_var[ii] / var_hat[ii];
    }
  }
  std::ptrdiff_t max_index;
  tmp.maxCoeff(&max_index);
  *return_index = max_index;

  // calculate primal step size
  if (tmp[max_index] == 0) {
    *return_step_length = 0;
  } else {
    *return_step_length = 1 / tmp[max_index];
  }
}

void lpengine::PrintStatusMsg(LPPRoblemStatus& status) {
  std::cout << "********************************************" << std::endl
            << "-- Problem Status : ";
  switch (status) {
    case STATUS_INFEASIBLE:
      std::cout << "INFEASIBLE" << std::endl;
      break;
    case STATUS_OPTIMAL:
      std::cout << "OPTIMAL" << std::endl
                << "-- number of itterations: " << p_.dic_number_ << std::endl
                // TODO: x_b_hat_ has different order of variables than ones Obj
                // func wants (it might be correct as well)
//                << "-- Optimal Objective Value: " << p_.obj_value_ << std::endl
                << "-- variable values:\n"<< std::endl
                << p_.x_b_hat_ << std::endl
                << p_.z_n_hat_ << std::endl;
      break;
    case STATUS_ITTERATING:
      std::cout << "SOLVING" << std::endl;
      break;
    case STATUS_UNBOUNDED:
      std::cout << "UNBOUNDED" << std::endl;
      break;
    case STATUS_NONE:
      std::cout << "NONE" << std::endl;
      break;
  }
}

bool lpengine::IsSemiPositive(Eigen::VectorXd& vec) {
  for (int ii = 0; ii < vec.rows(); ii++) {
    if (vec[ii] < 0) {
      return false;
    }
  }
  return true;
}

bool lpengine::IsNegative(Eigen::VectorXd& vec) {
  for (int ii = 0; ii < vec.rows(); ii++) {
    if (vec[ii] >= 0) {
      return false;
    }
  }
  return true;
}

void lpengine::SwapCols(Eigen::MatrixXd& matn, int index_matn,
                        Eigen::MatrixXd& matb, int index_matb) {
  Eigen::VectorXd tmp(matb.rows());
  tmp = matb.col(index_matb);
  matb.col(index_matb) = matn.col(index_matn);
  matn.col(index_matn) = tmp;
  int tmp_index;
  tmp_index = p_.basic_set_[index_matb];
  p_.basic_set_[index_matb] = p_.nonbasic_set_[index_matn];
  p_.nonbasic_set_[index_matn] = tmp_index;
}

LPPRoblemStatus lpengine::DualSimplexStep() {
  // Step1: Check if Dual is feasible then dictionary is optimal

  if (IsSemiPositive(p_.x_b_hat_)) {
    // current dictionary is optimal
    return LPPRoblemStatus::STATUS_OPTIMAL;
  }

  // Step2: Select Entring Variable
  p_.entering_var_ = ChooseEnteringVar(p_.x_b_hat_);

  // Step3: compute Dual step direction delta(z_n_)
  Eigen::MatrixXd tmp(p_.dim_m_, p_.dim_n_);
  tmp = p_.a_b_.inverse() * p_.a_n_;

  Eigen::VectorXd ei(p_.dim_m_);
  ei.setZero(ei.rows());
  ei[p_.entering_var_] = 1;
  p_.delta_z_n_ = -tmp.transpose() * ei;

  // Step4&5: Compute primal step length and leaving var
  ChooseLeavingVar(p_.delta_z_n_, p_.z_n_hat_, &p_.leaving_var_, &p_.s_);

  // if s_ is not positive then dual is unbounded (primal is infeasible)
  if (p_.s_ <= 0) {
    return LPPRoblemStatus::STATUS_UNBOUNDED;
  }

  // Step6: choose primal step direction. dual_enteringvar=dual_leavingvar
  Eigen::VectorXd ej(p_.dim_n_);
  ej.setZero(ej.rows());
  ej[p_.leaving_var_] = 1;
  p_.delta_x_b_ = tmp * ej;

  // Step7: Compute primal step length dual_leavingvar=primal_enteringvar
  p_.t_ = p_.x_b_hat_[p_.entering_var_] / p_.delta_x_b_[p_.entering_var_];

  // Step8: Update current primal and dual solutions
  p_.x_b_hat_ = p_.x_b_hat_ - p_.t_ * p_.delta_x_b_;
  p_.z_n_hat_ = p_.z_n_hat_ - p_.s_ * p_.delta_z_n_;

  // Step9: Update Basis
  SwapCols(p_.a_b_, p_.entering_var_,p_.a_n_ , p_.leaving_var_);

  p_.x_b_hat_[p_.entering_var_] = p_.t_;
  p_.z_n_hat_[p_.leaving_var_] = p_.s_;

  // Calc Objective Function Value
 // p_.obj_value_ = TODO: c_b_.transpose * a_b_.inverse() * b_;

  p_.dic_number_++;
  return LPPRoblemStatus::STATUS_ITTERATING;
}

LPPRoblemStatus lpengine::PrimalSimplexStep() {
  // Step1: Check if Dual is feasible then dictionary is optimal

  if (IsSemiPositive(p_.z_n_hat_)) {
    // current dictionary is optimal
    return LPPRoblemStatus::STATUS_OPTIMAL;
  }

  // Step2: Select Entring Variable
  p_.entering_var_ = ChooseEnteringVar(p_.z_n_hat_);

  // Step3: compute Primal step direction delta(x_b_)
  Eigen::MatrixXd tmp(p_.dim_m_, p_.dim_n_);
  tmp = p_.a_b_.inverse() * p_.a_n_;

  Eigen::VectorXd ej(p_.dim_n_);
  ej.setZero(ej.rows());
  ej[p_.entering_var_] = 1;
  p_.delta_x_b_ = tmp * ej;

  // Step4&5: Compute primal step length and leaving var
  ChooseLeavingVar(p_.delta_x_b_, p_.x_b_hat_, &p_.leaving_var_, &p_.t_);

  // if t_ is not positive then primal is unbounded
  if (p_.t_ <= 0) {
    return LPPRoblemStatus::STATUS_UNBOUNDED;
  }

  // Step6: choose dual Step Direction. dual_enteringvar=dual_leavingvar
  Eigen::VectorXd ei(p_.dim_m_);
  ei.setZero(ei.rows());
  ei[p_.leaving_var_] = 1;
  p_.delta_z_n_ = -tmp.transpose() * ei;

  // Step7: Compute dual step length dual_leavingvar=primal_enteringvar
  p_.s_ = p_.z_n_hat_[p_.entering_var_] / p_.delta_z_n_[p_.entering_var_];

  // Step8: Update current primal and dual solutions
  p_.x_b_hat_ = p_.x_b_hat_ - p_.t_ * p_.delta_x_b_;
  p_.z_n_hat_ = p_.z_n_hat_ - p_.s_ * p_.delta_z_n_;
  p_.x_b_hat_[p_.leaving_var_] = p_.t_;
  p_.z_n_hat_[p_.entering_var_] = p_.s_;

  // Step9: Update Basis
  SwapCols(p_.a_n_, p_.entering_var_, p_.a_b_, p_.leaving_var_);

  // Calc Objective Function Value
//  p_.obj_value_ = p_.x_b_hat_.head(p_.c_n_.rows()).dot(p_.c_n_);
  p_.dic_number_++;

  return LPPRoblemStatus::STATUS_ITTERATING;
}

LPPRoblemStatus lpengine::Solve() {
  LPPRoblemStatus current_status = LPPRoblemStatus::STATUS_NONE;

  if (IsSemiPositive(p_.b_) && IsNegative(p_.c_n_)) {
    std::cout << " -- Initial dictionary is Optimal" << std::endl;
    current_status = LPPRoblemStatus::STATUS_OPTIMAL;

  } else if (CheckDicfeasible(p_.x_b_hat_)) {
    std::cout << "-- Initial dictionary is primal feasible" << std::endl;
    p_.solver_type_ = LPSolverType::LP_SIMPLEX_PRIMAL;
    current_status = LPPRoblemStatus::STATUS_ITTERATING;
    while (current_status == LPPRoblemStatus::STATUS_ITTERATING) {
      current_status = PrimalSimplexStep();
    }

  } else if (CheckDicfeasible(p_.z_n_hat_)) {
    std::cout << "-- Initial dictionary is dual feasible" << std::endl;
    p_.solver_type_ = LPSolverType::LP_SIMPLEX_DUAL;
    current_status = LPPRoblemStatus::STATUS_ITTERATING;
    while (current_status == LPPRoblemStatus::STATUS_ITTERATING) {
      current_status = DualSimplexStep();
    }

  } else {
    std::cout << "-- Initial Primal and Dual dictionaries are not feasible"
              << std::endl
              << "-- Running PHASE1 algorithm" << std::endl;
    // modify the problem to feasible Dual problem and solve till optimal
    p_.ChangeToAuxiliary();
//    p_.PrintDictionary();
    std::cout << "--Aux dic is: ";
    if (CheckDicfeasible(p_.z_n_hat_)) {
      std::cout << "Dual feasible" << std::endl;

    } else {
      std::cout << "Error-> still not feasible" << std::endl;
      std::cout << "Exitting ..." << std::endl;
      current_status = LPPRoblemStatus::STATUS_INFEASIBLE;
      exit(EXIT_FAILURE);
    }
    p_.solver_type_ = LPSolverType::LP_SIMPLEX_DUAL;
    current_status = LPPRoblemStatus::STATUS_ITTERATING;
    while (current_status == LPPRoblemStatus::STATUS_ITTERATING) {
      current_status = DualSimplexStep();
    }

    p_.PrintDictionary();


  }

  // Problem is feasible and bounded at this point then do Simplex steps

  PrintStatusMsg(current_status);
  return current_status;
}
