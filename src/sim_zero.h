#include "sim.h"

#ifndef SIM_ZERO_H
#define SIM_ZERO_H

class Sim_zero : public Sim{
public:
    void schedule_cancer(double age_in, int year_in){
        double age_precursor, age_preclinical, age_clinical, age_death_other, age_death_cancer, age_death;
        int year_clinical = 9999, year_death = 9999;
        std::vector<double> res;

        // parameters
        double mean_death_other = (m1 + (year_in - age_in - 1988) * mt);

        age_death_other = gen_truncnorm(age_in, mean_death_other, sd_death_other);

        int i_cancer = 0;
        while(true){
            age_precursor = R::rlnorm(lm_precursor, ls_precursor);
            if(R::runif(0,1) < p_preclinical){
                age_preclinical = age_precursor + R::rexp(lam_preclinical);
            } else {
                age_preclinical = age_death_other + 1;
            }
            age_clinical = age_preclinical + R::rlnorm(lm_clinical, ls_clinical);
            age_death_cancer = age_clinical + R::rlnorm(lm_surv_cancer, ls_surv_cancer);
            if(age_death_cancer > age_in){
                break;
            }
            if(++i_cancer > 200){
                age_death_cancer = age_in;
                break;
            }
        }

        if(age_death_cancer < age_death_other){
            age_death = age_death_cancer;
            year_death = cal_year(age_in, year_in, age_death);
            map_cancer_death[find_lower_bound(years, year_death)] += 1;
        } else {
            age_death = age_death_other;
            year_death = cal_year(age_in, year_in, age_death);
            map_other_death[find_lower_bound(years, year_death)] += 1;
        }

        if(age_clinical < age_death){
            year_clinical = cal_year(age_in, year_in, age_clinical);
            map_cancer_incidence[find_lower_bound(years, year_clinical)] += 1;
        }
        /*
         if(year_death < year_in){
         Rcout << "age_in: " << age_in << " "
               << "age_clinical: " << age_clinical << " "
               << "age_death_cancer: " << age_death_cancer << " "
               << "age_death: " << age_death << "\n"
               << "year_clinical: " << year_clinical << " "
               << "year_death: " << year_death << "\n";
         }
         */
    }
};

#endif
