//[[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <math.h>
#include <fstream>
using namespace Rcpp;

//' @export
//[[Rcpp::export]]
void trial_read_csv_and_schedule_cancer(){
    std::string first_line, line;

    Environment env("package:base");
    Function sysfile = env["system.file"];
    CharacterVector vfname = sysfile("extdata", "pop_for_cancer.csv", Named("package", "microcancer"));
    std::string fname = as<std::string>(vfname[0]);
    Rcout << fname << "\n";

    std::ifstream myfile(fname);

    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            Rcout << line << "\n";
        }
        myfile.close();
    }

    else Rcout << "Unable to open file 1";
}
