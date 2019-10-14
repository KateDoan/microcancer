//[[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include "classA.h"
using namespace Rcpp;

class B{
public:
    B(int b1, int b2, int b3): b1(b1), b2(b2), b3(b3) {};
    int make_a_bigger_num(){
        A an_a = A(b1, b2);
        return an_a.bigger_than_a1();
    }
private:
    int b1, b2, b3;

};

//' @export
//[[Rcpp::export]]
void ab_work_together(){
    B a_b = B(1, 2, 3);
    Rcout << a_b.make_a_bigger_num() << "\n";
}

