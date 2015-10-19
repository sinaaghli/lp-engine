#include <lp-engine/lpengine.h>
#include <iostream>

int main() {
//////////////// Primal Feasible Problem
  lpproblem problem_Pfeasible;

  Eigen::MatrixXd A(3,2);
  A << 1.0f, -1.0f,
       2.0f, -1.0f,
       0.0f, 1.0f;

  Eigen::VectorXd b(3);
  b << 1.0f, 3.0f, 5.0f;

  Eigen::VectorXd c(2);
  c << 4.0f, 3.0f;

  problem_Pfeasible.Set_A(A);
  problem_Pfeasible.Set_b(b);
  problem_Pfeasible.Set_c(c);
//  problem_Pfeasible.UpdateProblem();

//////////////// Primal Feasible Problem2
  lpproblem problem_Pfeasible2;

  A.resize(3,2);
  A << -1.0f, 3.0f,
       1.0f, 1.0f,
       2.0f, -1.0f;

  b.resize(3);
  b << 12.0f, 8.0f, 10.0f;

  c.resize(2);
  c << 3.0f, 2.0f;

  problem_Pfeasible2.Set_A(A);
  problem_Pfeasible2.Set_b(b);
  problem_Pfeasible2.Set_c(c);
//  problem_Pfeasible2.UpdateProblem();

//////////////// Dual feasible Problem
  lpproblem problem_Dfeasible;

  A.resize(2,3);
  A << 1.0f, -1.0f, -1.0f,
       -3.0f, -1.0f, 1.0f;

  b.resize(2);
  b << -1.0f, -2.0f;

  c.resize(3);
  c << -10.0f, -6.0f, -2.0f;

  problem_Dfeasible.Set_A(A);
  problem_Dfeasible.Set_b(b);
  problem_Dfeasible.Set_c(c);
  problem_Dfeasible.UpdateProblem();

//////////////// Dual feasible Problem2
  lpproblem problem_Dfeasible2;

  A.resize(6,2);
  A <<  -2,7,
        -3,1,
        9,-4,
        1,-1,
        7,-3,
        -5,2;
  b.resize(6);
  b << 6.0f, -1.0f, 6.0f, 1, 6, -3;

  c.resize(2);
  c << -1.0f, -2.0f;

  problem_Dfeasible2.Set_A(A);
  problem_Dfeasible2.Set_b(b);
  problem_Dfeasible2.Set_c(c);
  problem_Dfeasible2.UpdateProblem();

//////////////// Dual feasible Problem3
  lpproblem problem_Dfeasible3;

  A.resize(3,2);
  A << -3.0f, -1.0f,
      -4.0f, -3.0f,
      -1.0f, -2.0f;

  b.resize(3);
  b << -3.0f, -6.0f, -3;

  c.resize(2);
  c << -29.0f, -10.0f;

  problem_Dfeasible3.Set_A(A);
  problem_Dfeasible3.Set_b(b);
  problem_Dfeasible3.Set_c(c);
  problem_Dfeasible3.UpdateProblem();

//////////////// Infeasible Problem
  lpproblem problem_infeasible;

  A.resize(3,2);
  A << -2.0f, -1.0f,
       -2.0f, 4.0f,
       -1.0f, 3.0f;

  b.resize(3);
  b << 4.0f, -8.0f, -7.0f;

  c.resize(2);
  c << -1.0f, 4.0f;

  problem_infeasible.Set_A(A);
  problem_infeasible.Set_b(b);
  problem_infeasible.Set_c(c);
  problem_infeasible.UpdateProblem();

//////////////// Cycling Problem
  lpproblem problem_cycling;

  A.resize(3,4);
  A <<   0.5f, -5.5f, -2.5f,  9.0f,
         0.5f, -1.5f, -0.5f,  1.0f,
         1.0f,  0.0f,  0.0f,  0.0f;

  b.resize(3);
  b <<  0.0f,  0.0f,  1.0f;

  c.resize(4);
  c << 10.0f, -57.0f, -9.0f, -24.0f;

  problem_cycling.Set_A(A);
  problem_cycling.Set_b(b);
  problem_cycling.Set_c(c);
//  problem_cycling.UpdateProblem();

  // Enable the requred problem to be solved in the following lines
  lpengine my_engine(problem_Pfeasible);
//  lpengine my_engine(problem_Pfeasible2);
//  lpengine my_engine(problem_Dfeasible);
//  lpengine my_engine(problem_Dfeasible2);
//  lpengine my_engine(problem_Dfeasible3);
//  lpengine my_engine(problem_infeasible);
//  lpengine my_engine(problem_cycling);

  my_engine.Solve();

  return 0;
}
