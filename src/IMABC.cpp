// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "util_IMABC.h"
using namespace Rcpp;
bool is_trial = true;

class IMABC {
public:
    IMABC(size_t size_theta,
          std::vector<double> target,
          std::vector<double> target_sd,
          std::vector<double> alpha_start,
          std::vector<double> alpha_goal,
          int N0, int Nc, int Ngoal, int B, int LIM1, int LIM2, int LIM3) :
        size_theta(size_theta),
        target(target), target_sd(target_sd),
        alpha_vec(alpha_start), alpha_goal(alpha_goal),
        N0(N0), Nc(Nc), Ngoal(Ngoal), B(B), LIM1(LIM1), LIM2(LIM2), LIM3(LIM3) {
        ret = std::vector<Point>();
        theta_vec = std::vector<std::vector<double>>();
        sigma_vec = std::vector<NumericMatrix>();
        running_Nt = 0;
        n_selected_pts = 0;
        ESS = 0.0;
    };

    double d_prior(std::vector<double> x){
        double d_prior_val = R::punif(x[0], 70, 80, false, false)
                        + R::punif(x[1], 0, 1, false, false)
                        + R::punif(x[2], 0, 20, false, false)
                        + R::punif(x[3], 0, 4.5, false, false)
                        + R::punif(x[4], 0, 2, false, false)
                        + R::punif(x[5], 0, 30, false, false)
                        + R::punif(x[6], 0, 0.5, false, false)
                        + R::punif(x[7], 0, 3.5, false, false)
                        + R::punif(x[8], 0, 2.0, false, false)
                        + R::punif(x[9], 0, 3.5, false, false)
                        + R::punif(x[10], 0, 2.0, false, false);

        return d_prior_val;
    }

    // Microsimulation function
    std::vector<double> fsim(std::vector<double> x){
        return x;
    }

    double var_unif(double a, double b){
        return (b-a)*(b-a)/12;
    }

    void step1(){
        std::vector<Point>().swap(ret);
        std::vector<double> params(size_theta);
        double min_p, dist2_to_target;

        for(int it=0; it<N0; it++){
            /*
            for(int i=0; i<size_theta; ++i){
                params[i] = R::rnorm(3, 1);
            }
            */
            params[0] = R::runif(70, 80);
            params[1] = R::runif(0, 1);
            params[2] = R::runif(0, 20);
            params[3] = R::runif(0, 4.5);
            params[4] = R::runif(0, 2);
            params[5] = R::runif(0, 30);
            params[6] = R::runif(0, 0.5);
            params[7] = R::runif(0, 3.5);
            params[8] = R::runif(0, 2);
            params[9] = R::runif(0, 3.5);
            params[10] = R::runif(0, 2);

            std::vector<double> curr_sim_res = fsim(params);
            std::vector<double> pval_vec = cal_pval_vec(curr_sim_res, target, target_sd);
            min_p = *(std::min_element(pval_vec.begin(), pval_vec.end()));
            dist2_to_target = squared_dist_vec(curr_sim_res, target);
            ret.push_back(Point(params, min_p, dist2_to_target, pval_vec));
        }

        running_Nt += N0;
    }

    void step2(){
        std::sort(ret.begin(), ret.end(), better_fit_target());
    }

    void step3(){
        std::vector<Point> ret_copy(ret);
        std::vector<Point> new_points = std::vector<Point>();
        std::vector<double> curr_sim_res;
        NumericMatrix sig_mat;
        NumericMatrix new_vecs;
        int i;
        std::vector<Point>::iterator it, it_copy, it_begin, it_Nc;

        Environment mvrnorm_env("package:MASS");
        Function mvrnorm = mvrnorm_env["mvrnorm"];
        Function var("var");
        n_selected_pts = ret.size();

        if(!(n_selected_pts > LIM1*size_theta)){
            sig_mat = NumericMatrix::diag(size_theta, 1);
            sig_mat(0, 0) = var_unif(70, 80)/4;
            sig_mat(1, 1) = var_unif(0, 1)/4;
            sig_mat(2, 2) = var_unif(0, 20)/4;
            sig_mat(3, 3) = var_unif(0, 4.5)/4;
            sig_mat(4, 4) = var_unif(0, 2)/4;
            sig_mat(5, 5) = var_unif(0, 30)/4;
            sig_mat(6, 6) = var_unif(0, 0.5)/4;
            sig_mat(7, 7) = var_unif(0, 3.5)/4;
            sig_mat(8, 8) = var_unif(0, 2)/4;
            sig_mat(9, 9) = var_unif(0, 3.5)/4;
            sig_mat(10, 10) = var_unif(0, 2)/4;
        } else if((n_selected_pts > LIM1*size_theta) && (!(n_selected_pts > LIM2*size_theta))) {
            NumericMatrix param_mat(ret.size(), size_theta);
            for(i=0, it=ret.begin(); it!=ret.end(); it++, i++){
                NumericVector curr_params = wrap((*it).params);
                param_mat.row(i) = curr_params;
            }
            sig_mat = var(param_mat);
        } else {
        }

        std::vector<double> m, curr_new_vec, pval_vec;
        double min_pval, dist2_to_target;

        it_begin=ret.begin();
        it_Nc= ret.size()<Nc ? ret.end() : ret.begin()+Nc;

        for(std::vector<Point>::iterator it = it_begin; it < it_Nc; it++){
            m = (*it).params;
            curr_sim_res = fsim(m);
            pval_vec = cal_pval_vec(curr_sim_res, target, target_sd);
            min_pval = *(std::min_element(pval_vec.begin(), pval_vec.end()));
            dist2_to_target = squared_dist_vec(curr_sim_res, target);
            if(is_acceptable(pval_vec, alpha_vec)){
                new_points.push_back(Point(m, min_pval, dist2_to_target, pval_vec));
            }

            if(n_selected_pts > LIM2*size_theta){
                auto ref_point = (*it);
                int n_nearest = LIM2 * size_theta;

                std::sort(ret_copy.begin(), ret_copy.end(), less_dist_to_ref(ref_point));

                NumericMatrix param_mat(n_nearest, size_theta);

                for(i=0, it_copy=ret_copy.begin(); it_copy!=ret_copy.begin() + n_nearest; it_copy++, i++){
                    NumericVector curr_params = wrap((*it_copy).params);
                    param_mat.row(i) = curr_params;
                }

                sig_mat = var(param_mat);
            }

            if(B==1){
                new_vecs = NumericMatrix(1, size_theta);
                NumericVector temp = mvrnorm(B, m, sig_mat);
                new_vecs.row(0) = temp;
            } else {
                new_vecs = mvrnorm(B, m, sig_mat);
            }
            running_Nt += B;

            //Store the theta vec and sigma vec for the centre points
            theta_vec.push_back(m);
            sigma_vec.push_back(sig_mat);

            //Calculate the pval_vec and see if the point is acceptable
            for(i=0; i<B; i++){
                NumericVector temp_new_vec = new_vecs.row(i);
                curr_new_vec = as<std::vector<double>> (temp_new_vec);
                curr_sim_res = fsim(curr_new_vec);
                pval_vec = cal_pval_vec(curr_sim_res, target, target_sd);
                min_pval = *(std::min_element(pval_vec.begin(), pval_vec.end()));
                dist2_to_target = squared_dist_vec(curr_sim_res, target);

                if(is_acceptable(pval_vec, alpha_vec)){
                    new_points.push_back(Point(curr_new_vec, min_pval, dist2_to_target, pval_vec));
                }
            }
        }
        ret.erase(it_begin, it_Nc);
        ret.insert(ret.end(), new_points.begin(), new_points.end());
        n_selected_pts = ret.size();

        if(is_trial) Rcout << "running_Nt = " << running_Nt << "\n";
        if(is_trial) Rcout << "number of selected points: " << n_selected_pts << "\n";
    }

    // Calculate ESS
    void step4(){
        //double curr_d_prior, curr_dmixture;
        double curr_weight = 1;
        arma::mat params_mat(ret.size(), size_theta);
        arma::vec dprior(ret.size());
        std::vector<Point>::iterator it;
        int i;

        for(it = ret.begin(), i=0; it != ret.end(); it++, i++){
            std::vector<double> curr_params = (*it).params;
            dprior[i] = d_prior(curr_params);
            params_mat.row(i) = as<arma::vec>(wrap(curr_params)).t();
        }

        //Rf_PrintValue(wrap(dprior));

        arma::vec sum_dmixture = arma_sum_dnorm(params_mat, theta_vec, sigma_vec);
        arma::vec dmixture = (sum_dmixture*B + dprior * N0)/(B*theta_vec.size() + N0);
        arma::vec weight = dprior/dmixture;
        ESS = 1.0/arma::accu(arma::pow(weight, 2));
    }

    void IMABC_main(){
        step1();

        while(alpha_vec < alpha_goal){
            if(n_selected_pts > LIM3 * size_theta){
                std::sort(ret.begin(), ret.end(), better_fit_target());
                int median_idx = int(n_selected_pts/2);
                std::vector<Point>::iterator it = ret.begin() + median_idx;
                std::vector<double>pvals_median = (*it).pvals;
                for(int i=0; i<alpha_goal.size(); i++){
                    alpha_vec[i] = std::min(pvals_median[i], alpha_goal[i]);
                }
                if(is_trial) Rcout << n_selected_pts << "\n";
                int n_preserved = int(n_selected_pts/4);
                ret.erase(std::remove_if(ret.begin() + n_preserved, ret.end(), bad_point(alpha_vec)),
                          ret.end());
                n_selected_pts = ret.size();
                if(is_trial) Rcout << n_selected_pts << "\n";
            } else {
                step2();
                step3();
            }
        }

        step4();

        while(ESS < Ngoal){
            step2();
            step3();
            step4();
        }

    }

    // data
    size_t size_theta;
    std::vector<double> target;
    std::vector<double> target_sd;
    std::vector<double> alpha_vec;
    std::vector<double> alpha_goal;
    int N0, Nc, Ngoal, B, LIM1, LIM2, LIM3, running_Nt=0, n_selected_pts=0;
    double ESS = 0.0;

    std::vector<Point> ret;
    std::vector<std::vector<double>> theta_vec;
    std::vector<NumericMatrix> sigma_vec;
};

//' @export
//[[Rcpp::export]]
void create_IMABC(){
    size_t size_theta = 11;
    int N0=1e5, Nc=10, Ngoal=2000, B=1000, LIM1=5, LIM2=25, LIM3=50;
    std::vector<double>target{75, 0.5, 10, 2, 1, 15, 0.25, 1.75, 1, 1.75, 1};
    std::vector<double>target_sd(size_theta, 2);
    std::vector<double>alpha_start(size_theta, 0.01);
    std::vector<double>alpha_test(size_theta, 0.02);
    std::vector<double>alpha_goal(size_theta, 0.7);

    IMABC my_imabc = IMABC(size_theta,
                           target,
                           target_sd,
                           alpha_start,
                           alpha_goal,
                           N0, Nc, Ngoal, B, LIM1, LIM2, LIM3);
    my_imabc.IMABC_main();
    Rcout << "ESS = " << my_imabc.ESS << "\n";
    print_point_vec(my_imabc.ret);
    //write_csv_vpoints(my_imabc.ret, "./output/myparams.csv", size_theta, target.size());
}





