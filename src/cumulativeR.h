#ifndef CUMULATIVER_H_
#define CUMULATIVER_H_
#include "distribution.h"

class CumulativeR : virtual public distribution, Logistic, Normal, Cauchit, Student, Gumbel, Gompertz{
public:
  CumulativeR();

  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta1) const;
  virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_normal(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_normal(const Eigen::VectorXd& eta) const;

  virtual Eigen::VectorXd inverse_cauchit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_cauchit(const Eigen::VectorXd& eta) const;

  virtual Eigen::VectorXd inverse_gompertz(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_gompertz(const Eigen::VectorXd& eta) const ;


  // List GLMcum(std::string response,
  //             StringVector explanatory_complete,
  //             StringVector explanatory_proportional,
  //             std::string distribution,
  //             SEXP categories_order,
  //             DataFrame dataframe,
  //             StringVector beta_t,
  //             Eigen::VectorXd beta_init);


};

#endif


