#include "distribution.h"
#include "reference.h"
using namespace std;
using namespace Rcpp ;

// [[Rcpp::depends(RcppEigen)]]
ReferenceF::ReferenceF(void) {
  distribution dist;
}

Eigen::VectorXd ReferenceF::inverse_logistic(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(size_t j=0; j<eta.size(); ++j)
  {
    pi[j] = Logistic::cdf_logit( eta(j) ) / ( 1-Logistic::cdf_logit( eta(j) ) );
    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::VectorXd ReferenceF::inverse_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(size_t j=0; j<eta.size(); ++j)
  {
    pi[j] = cdf_normal( eta(j) ) / ( 1-cdf_normal( eta(j) ) );
    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_logistic(const Eigen::VectorXd& eta2 ) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_logistic(eta2);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  { D1(j,j) = pdf_logit( eta2(j) ) /
    (Logistic::cdf_logit(eta2(j)) * (1-Logistic::cdf_logit(eta2(j))));
  }
  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}

Eigen::MatrixXd ReferenceF::inverse_derivative_normal(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi = ReferenceF::inverse_normal(eta);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pi.rows(),pi.rows());
  for(size_t j=0; j<pi.rows(); ++j)
  { D(j,j) = pdf_normal( eta(j) ) /
    ( std::max(1e-10, std::min(1-1e-6,cdf_normal(eta(j)))) *
      std::max(1e-10, std::min(1-1e-6, 1-cdf_normal(eta(j)))) ); }
  return D * ( Eigen::MatrixXd(pi.asDiagonal()) - pi * pi.transpose().eval() );
}

Eigen::VectorXd ReferenceF::inverse_cauchit(const Eigen::VectorXd& eta) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(size_t j=0; j<eta.size(); ++j)
  {
    pi[j] = Cauchit::cdf_cauchit( eta(j) ) / ( 1-Cauchit::cdf_cauchit( eta(j) ) );
    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_cauchit(const Eigen::VectorXd& eta2) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_cauchit(eta2);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  { D1(j,j) = pdf_cauchit( eta2(j) ) /
    (Cauchit::cdf_cauchit(eta2(j)) * (1-Cauchit::cdf_cauchit(eta2(j))));
  }
  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}

Eigen::VectorXd ReferenceF::inverse_student(const Eigen::VectorXd& eta, double freedom_degrees) const
{
  Eigen::VectorXd pi( eta.size() );
  double norm1 = 1.;
  for(size_t j=0; j<eta.size(); ++j)
  {
    pi[j] = Student::cdf_student( eta(j) , freedom_degrees) / ( 1-Student::cdf_student( eta(j) , freedom_degrees ) );
    norm1 += pi[j];
  }
  return (pi/norm1);
}

Eigen::MatrixXd ReferenceF::inverse_derivative_student(const Eigen::VectorXd& eta2, double freedom_degrees) const
{
  Eigen::VectorXd pi1 = ReferenceF::inverse_student(eta2, freedom_degrees);
  Eigen::MatrixXd D1 = Eigen::MatrixXd::Zero(pi1.rows(),pi1.rows());
  for(int j=0; j<eta2.rows(); ++j)
  { D1(j,j) = pdf_student( eta2(j) , freedom_degrees) /
    (Student::cdf_student(eta2(j),freedom_degrees) * (1-Student::cdf_student(eta2(j),freedom_degrees)));
  }
  return D1 * ( Eigen::MatrixXd(pi1.asDiagonal()) - pi1 * pi1.transpose().eval() );
}

distribution dist1;

// [[Rcpp::export]]
List GLMref(std::string response,
            StringVector explanatory_complete,
            StringVector explanatory_proportional,
            std::string distribution,
            SEXP categories_order,
            DataFrame dataframe,
            double freedom_degrees = 2.0){

  int P_c = 0;
  if(explanatory_complete[0] != "NA"){P_c = explanatory_complete.size(); }
  int P_p = 0;
  if(explanatory_proportional[0] != "NA"){P_p = explanatory_proportional.size(); }
  int P =  P_c +  P_p ; // Number of explanatory variables without intercept

  const int N = dataframe.nrows() ; // Number of observations



  List Full_M = dist1.select_data(dataframe, response, explanatory_complete,
                                  explanatory_proportional, categories_order);

  Eigen::MatrixXd Y_init = Full_M["Y_ext"];
  Eigen::MatrixXd X_EXT = Full_M["X_EXT"];
  // CharacterVector levs1 = Full_M["levs1"];

  int Q = Y_init.cols();
  int K = Q + 1;
  // // // Beta initialization with zeros
  Eigen::MatrixXd BETA;
  BETA = Eigen::MatrixXd::Zero(X_EXT.cols(),1);
  //
  int iteration = 0;
  // double check_tutz = 1.0;
  double Stop_criteria = 1.0;
  Eigen::MatrixXd X_M_i ;
  Eigen::VectorXd Y_M_i ;
  Eigen::VectorXd eta ;
  Eigen::VectorXd pi ;
  Eigen::MatrixXd D ;
  Eigen::MatrixXd Cov_i ;
  Eigen::MatrixXd W_in ;
  Eigen::MatrixXd Score_i_2 ;
  Eigen::MatrixXd F_i_2 ;
  Eigen::VectorXd LogLikIter;
  LogLikIter = Eigen::MatrixXd::Zero(1,1) ;
  Eigen::MatrixXd var_beta;
  Eigen::VectorXd Std_Error;
  double LogLik;
  Eigen::MatrixXd pi_ma(N, K);
  Eigen::MatrixXd F_i_final = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());

  // for (int iteration=1; iteration < 18; iteration++){
  // while (check_tutz > 0.0001){
  double epsilon = 0.0001 ;
  while (Stop_criteria >( epsilon / N)){

    Eigen::MatrixXd Score_i = Eigen::MatrixXd::Zero(BETA.rows(),1);
    Eigen::MatrixXd F_i = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());
    LogLik = 0.;

    // Loop by subject
    for (int i=0; i < N; i++){
      // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_init.row(i);
      eta = X_M_i * BETA;

      ReferenceF ref;

      // Vector pi depends on selected distribution
      if(distribution == "logistic"){
        pi = ref.inverse_logistic(eta);
        D = ref.inverse_derivative_logistic(eta);
      }else if(distribution == "normal"){
        pi = ref.inverse_normal(eta);
        D = ref.inverse_derivative_normal(eta);
      }else if(distribution == "cauchit"){
        pi = ref.inverse_cauchit(eta);
        D = ref.inverse_derivative_cauchit(eta);
      }else if(distribution == "student"){
        pi = ref.inverse_student(eta, freedom_degrees);
        D = ref.inverse_derivative_student(eta, freedom_degrees);
      }

      Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
      W_in = D * Cov_i.inverse();
      Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
      F_i = F_i + F_i_2;
      LogLik = LogLik + (Y_M_i.transpose().eval()*Eigen::VectorXd(pi.array().log())) + ( (1 - Y_M_i.sum()) * std::log(1 - pi.sum()) );

      pi_ma.row(i) = pi.transpose();

    }

    Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(pi_ma.rows());
    pi_ma.col(Q) = Ones1 - pi_ma.rowwise().sum() ;

    LogLikIter.conservativeResize(iteration+2, 1);
    LogLikIter(iteration+1) = LogLik;
    Stop_criteria = (abs(LogLikIter(iteration+1) - LogLikIter(iteration))) / (epsilon + (abs(LogLikIter(iteration+1)))) ;
    Eigen::VectorXd beta_old = BETA;
    BETA = BETA + (F_i.inverse() * Score_i);
    // check_tutz = ((BETA - beta_old).norm())/(beta_old.norm()+check_tutz);
    iteration = iteration + 1;
    F_i_final = F_i;
  }

  // var_beta = (((X_EXT.transpose() * F_i_final) * X_EXT).inverse());
  var_beta = F_i_final.inverse();
  Std_Error = var_beta.diagonal();
  Std_Error = Std_Error.array().sqrt() ;

  std::vector<std::string> text=as<std::vector<std::string>>(explanatory_complete);
  std::vector<std::string> level_text=as<std::vector<std::string>>(categories_order);
  StringVector names(Q*P_c + P_p);
  if(P_c > 0){
    for(int var = 0 ; var < explanatory_complete.size() ; var++){
      for(int cat = 0 ; cat < Q ; cat++){
        names[(Q*var) + cat] = dist1.concatenate(text[var], level_text[cat]);
      }
    }
  }
  if(P_p > 0){
    for(int var_p = 0 ; var_p < explanatory_proportional.size() ; var_p++){
      names[(Q*P_c) + var_p] = explanatory_proportional[var_p];
    }
  }
  // TO NAMED THE RESULT BETAS
  NumericMatrix coef = wrap(BETA);
  rownames(coef) = names;

  // AIC
  double AIC = (-2*LogLik) + (2 *coef.length());

  // AIC
  double BIC = (-2*LogLik) + (coef.length() * log(N) );

  int df = (N*Q) - coef.length();

  Eigen::MatrixXd predicted = X_EXT * BETA;

  Eigen::VectorXd Ones2 = Eigen::VectorXd::Ones(Y_init.rows());
  Eigen::VectorXd vex1 = (Y_init.rowwise().sum()) ;
  Y_init.conservativeResize( Y_init.rows(), K);
  Y_init.col(Q) = (vex1 - Ones2).array().abs() ;

  Eigen::MatrixXd residuals = Y_init - pi_ma;

  Eigen::VectorXd pi_ma_vec(Eigen::Map<Eigen::VectorXd>(pi_ma.data(), pi_ma.cols()*pi_ma.rows()));
  Eigen::VectorXd Y_init_vec(Eigen::Map<Eigen::VectorXd>(Y_init.data(), Y_init.cols()*Y_init.rows()));

  Eigen::VectorXd div_arr = Y_init_vec.array() / pi_ma_vec.array();

  Eigen::VectorXd dev_r(Y_init.rows());

  int el_1 = 0;

  for (int element = 0 ; element < div_arr.size() ;  element++){
    if (div_arr[element] != 0){
      dev_r[el_1] = div_arr[element];
      el_1 = el_1 +1 ;
    }
  }

  Eigen::ArrayXd dev_log = dev_r.array().log();

  double deviance = dev_log.sum();
  deviance = -2*deviance;

  return List::create(
    // Named("Nb. iterations") = iteration-1 ,
    Named("coefficients") = coef,
    Named("AIC") = AIC,
    Named("BIC") = BIC,
    // Named("var_beta") = var_beta,
    Named("stderr") = Std_Error,
    Rcpp::Named("df") = df,
    Rcpp::Named("predicted") = predicted,
    Rcpp::Named("fitted") = pi_ma,
    Rcpp::Named("pi_ma_vec") = pi_ma_vec,
    Rcpp::Named("Y_init_vec") = Y_init_vec,
    Rcpp::Named("dev_log") = dev_log,
    Rcpp::Named("deviance") = deviance,
    Rcpp::Named("residuals") = residuals,
    Named("Log-likelihood") = LogLik
  );
}

List ReferenceF::GLMref_ec(std::string response, std::string actual_response,
                           std::string individuals,
                           StringVector explanatory_complete,
                           StringVector depend_y,
                           std::string distribution,
                           SEXP categories_order,
                           DataFrame dataframe,
                           std::string design,
                           double freedom_degrees){
  int P_c = 0;
  if(explanatory_complete[0] != "NA"){P_c = explanatory_complete.size(); }
  int P_p = 0;
  if(depend_y[0] != "NA"){P_p = depend_y.size(); }

  const int N = dataframe.nrows() ; // Number of observations

  List Full_M = distribution::select_data_nested(dataframe, response, actual_response,
                                                 individuals, explanatory_complete, depend_y , categories_order);

  Eigen::MatrixXd Y_init = Full_M["Y_ext"];
  Eigen::MatrixXd X_EXT;
  if(design == "tutz"){
    X_EXT = Full_M["X_M_dep_y"];
  }else if(design == "louviere"){
    X_EXT = Full_M["X_M_dep_y_alt"];
  }

  int Q = Y_init.cols();
  int K = Q + 1;
  // // // Beta initialization with zeros
  Eigen::MatrixXd BETA;
  BETA = Eigen::MatrixXd::Zero(X_EXT.cols(),1);

  int iteration = 0;
  // double check_tutz = 1.0;
  double Stop_criteria = 1.0;
  Eigen::MatrixXd X_M_i ;
  Eigen::VectorXd Y_M_i ;
  Eigen::VectorXd eta ;
  Eigen::VectorXd pi ;
  Eigen::MatrixXd D ;
  Eigen::MatrixXd Cov_i ;
  Eigen::MatrixXd W_in ;
  Eigen::MatrixXd Score_i_2 ;
  Eigen::MatrixXd F_i_2 ;
  Eigen::VectorXd LogLikIter;
  LogLikIter = Eigen::MatrixXd::Zero(1,1) ;
  double LogLik;
  // for (int iteration=1; iteration < 18; iteration++){
  // while (check_tutz > 0.0001){
  double epsilon = 0.0001 ;
  while (Stop_criteria >( epsilon / N)){

    Eigen::MatrixXd Score_i = Eigen::MatrixXd::Zero(BETA.rows(),1);
    Eigen::MatrixXd F_i = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());
    LogLik = 0.;

    // Loop by subject
    for (int i=0; i < N/K; i++){
      // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_init.row(i);
      eta = X_M_i * BETA;

      // Vector pi depends on selected distribution
      if(distribution == "logistic"){
        pi = ReferenceF::inverse_logistic(eta);
        D = ReferenceF::inverse_derivative_logistic(eta);
      }else if(distribution == "normal"){
        pi = ReferenceF::inverse_normal(eta);
        D = ReferenceF::inverse_derivative_normal(eta);
      }else if(distribution == "cauchit"){
        pi = ReferenceF::inverse_cauchit(eta);
        D = ReferenceF::inverse_derivative_cauchit(eta);
      }else if(distribution == "student"){
        pi = ReferenceF::inverse_student(eta, freedom_degrees);
        D = ReferenceF::inverse_derivative_student(eta, freedom_degrees);
      }

      Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
      W_in = D * Cov_i.inverse();
      Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
      F_i = F_i + F_i_2;
      LogLik = LogLik + (Y_M_i.transpose().eval()*Eigen::VectorXd(pi.array().log())) + ( (1 - Y_M_i.sum()) * std::log(1 - pi.sum()) );
    }

    LogLikIter.conservativeResize(iteration+2, 1);
    LogLikIter(iteration+1) = LogLik;
    Stop_criteria = (abs(LogLikIter(iteration+1) - LogLikIter(iteration))) / (epsilon + (abs(LogLikIter(iteration+1)))) ;
    Eigen::VectorXd beta_old = BETA;
    BETA = BETA + (F_i.inverse() * Score_i);
    // check_tutz = ((BETA - beta_old).norm())/(beta_old.norm()+check_tutz);
    iteration = iteration + 1;
  }
  CharacterVector levs1 = Full_M["levs1"];

  std::vector<std::string> text=as<std::vector<std::string>>(explanatory_complete);
  std::vector<std::string> level_text=as<std::vector<std::string>>(categories_order);
  StringVector names(Q*P_c + P_p);
  if(P_c > 0){
    for(int var = 0 ; var < explanatory_complete.size() ; var++){
      for(int cat = 0 ; cat < Q ; cat++){
        names[(Q*var) + cat] = distribution::concatenate(text[var], level_text[cat]);
      }
    }
  }
  if(P_p > 0){
    for(int var_p = 0 ; var_p < depend_y.size() ; var_p++){
      names[(Q*P_c) + var_p] = depend_y[var_p];
    }
  }

  Rcout << names << std::endl;

  // TO NAMED THE RESULT BETAS
  NumericMatrix BETA_2 = wrap(BETA);
  // rownames(BETA_2) = names;

  return List::create(Named("Nb. iterations") = iteration-1 , Named("Coefficients") = BETA_2,
                      Named("Log-likelihood") = LogLik);

}


// RCPP_MODULE(referencemodule){
//   Rcpp::function("GLMref", &GLMref,
//                  List::create(_["response"] = "a",
//                               _["explanatory_complete"] = CharacterVector::create( "A", NA_STRING),
//                               _["explanatory_proportional"] = CharacterVector::create( "A", NA_STRING),
//                               _["distribution"] = "a",
//                               _["categories_order"] = R_NaN,
//                               _["dataframe"] = NumericVector::create( 1, NA_REAL, R_NaN, R_PosInf, R_NegInf),
//                               _["freedom_degrees"] = 1.0),
//                               "Provides a simple vector norm");
//   Rcpp::class_<ReferenceF>("ReferenceF")
//   .constructor()
//   // .method( "GLMref", &ReferenceF::GLMref)
//      .method( "GLMref_ec", &ReferenceF::GLMref_ec )
//      .method( "inverse_logistic", &ReferenceF::inverse_logistic )
//   ;
// }
