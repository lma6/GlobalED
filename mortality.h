#ifndef EDM_MORTALITY_H_
#define EDM_MORTALITY_H_

struct cohort;

void cohort_modifications_from_disturbance(int q, cohort** pcurrentc, 
                                           UserData* data);


#endif // EDM_MORTALITY_H_ 
