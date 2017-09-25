        #include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "edmodels.h"
#include "site.h"
#include "patch.h"
#include "read_site_data.h"

#ifdef ED
////////////////////////////////////////////////////////////////////////////////
//! Update_Water
//! When using the Dwdt formula for water, it was found that often the water value 
//! oscillated on a substep, sometimes with oscillations growing wlarger and larger 
//! until other errors occurred. To fix this, we now do a standard euler step if it
//! doesn't take us past the equilibrium, otherwise we find the equilibrium 
//! (min 0, max theta_max) and set the new water value to half way to equilibrium
//! (half way chosen to give trees a chance to change uptake in next substep)
//!   
//! Take standard time step, accept if haven't over shot equilibrium (derivative dwdt 
//! still same sign). Note that this function assumes that the function dwdt is smooth 
//! and monotonically increasing w.r.t water. For most reasonable functions for dwdt 
//! this is logical, but changes to dwdt function must be made bearing this in mind.
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void patch::Update_Water(double time, UserData* data, double deltat){
   double starting_dwdt = dwdt; 
   double starting_water = water;
   water += dwdt*deltat;
   if (Dwdt(time, data)*starting_dwdt>0) return;
   
   //Otherwise find approximation to equilibrium. Note this is not true equilibrium as trees have not had a chance to adjust.
   site* currents = siteptr;
   
   double theta_tol = 0.001; //Acceptable error in % saturation (difference off equilibrium)
   double water_tol = theta_tol*currents->sdata->soil_depth * currents->sdata->theta_max;
   water = starting_water; dwdt = starting_dwdt;
   
   double max_guess, min_guess;
   if (dwdt>0){
      min_guess = water; 
      max_guess = water + dwdt * deltat;
   } else {
      min_guess = water - dwdt * deltat; 
      max_guess = water;
   }
   //Binomial search
   while (max_guess-min_guess>water_tol){
      water = (max_guess+min_guess)/2;
      if (Dwdt(time, data)>0) min_guess = water;
      else max_guess = water;
   }
   //Adjust half way to equilibrium
   water = (water + starting_water) / 2.;
}

////////////////////////////////////////////////////////////////////////////////
//! Dwdt
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double patch::Dwdt (double time, UserData* data){  

   /* This function, and those that it calls, must use *w ! */
   /* this function calculates and updates dwdt a patch */
   //For Update_Water to work this function must be monotonically increasing w.r.t water
   site* currents = siteptr;
   
   /* calculate water evaporation from the soil, scale by plant cover */
   soil_evap = (currents->sdata->soil_evap_conductivity)
      * (radiative_flux(data)) * water;
  
   /* calculate water loss per unit area from patch */
   theta = water / (currents->sdata->soil_depth * currents->sdata->theta_max);
   if (water > 0.0) {
      perc = currents->sdata->k_sat * pow(theta, 2.0 * currents->sdata->tau + 2.0);
   } else {
      perc = 0.0;
   }
   dwdt = (currents->sdata->precip[(int) data->time_period] 
           - perc - total_water_uptake
           / area) - soil_evap;
   return dwdt;
}



////////////////////////////////////////////////////////////////////////////////
//! Dsdt
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void patch::Dsdt (unsigned int time_period, double time, UserData* data) {

   /*** This function must set time derivatives for each soil pool ***/
   /* Simplified 4-Box Century, based on Parton et al 1987 and 1993  */
   /* 4 carbon boxes + 2 nitrogen boxes                              */
   /* There is no: leaching, other n inputs, or effect of soil type  */
   double fast_C_loss=0.00,  fast_N_loss=0.00;
   double structural_C_loss=0.00;
   double slow_C_input=0.00,  slow_C_loss=0.00, slow_N_loss=0.00;
   double mineralized_N_input=0.00, mineralized_N_loss=0.00;
   double N_immobilization_demand=0.00;
   double passive_C_input=0.00, passive_C_loss=0.00, passive_N_loss=0.00;
   double structural_L_loss=0.00;

   double Ls; /* the fraction of structural material that is lignin */
   double Lc; /* decomp rate reduction due to lignin */

   /* CENTURY PARAM VALUES */
   /* values are simple averages from lumped Century Pools */
   /* Based on Parton et al 1993 GBC */
   /* 1-structural,2-fast,3=slow,4=passive */
   double K1=4.5, K2=11.0, K3=100.2, K4=0.0;  /* Max decay rate yr^-1;  *
                                               * from Century directly, *
                                               * or averages of lumped  *
                                               * Century pools          */
   /* K3 is high bc we wanted to added back the n immobilization story 
      without tracking the slow pool */

   /* std values= 1, .3, 1, 0 */
   double r_fsc=1.0, r_stsc=0.3, r_ssc=1.0, r_psc=0.0; 
   /* respiration rates of soil pools */
   
   /****************/
   /* OTHER FLUXES */
   /****************/
   /* calculate commonly used terms once */
   A = A_function(time_period, data);

   /* Slow down of decomp due to fraction of structural material that is lignin */
   /* From CENTURY Parton et al 1993 GBC 7(4):785-809 */
      Ls = structural_soil_L / structural_soil_C;

      Lc = exp(-3.0 * Ls);

   if(!data->n_competition) {
      total_plant_nitrogen_uptake = 0.0;
   }

   mineralized_N_loss = total_plant_nitrogen_uptake;

   if(data->n_decomp_limitation) {
      /* N immobilized by the decomp of structural material */
      N_immobilization_demand= A * Lc * K1 * structural_soil_C
         * ((1.0 - r_stsc) * (1.0 / data->c2n_slow) - (1.0 / data->c2n_structural));
      fstd = data->nitrogen2 * mineralized_soil_N
         / (N_immobilization_demand + data->nitrogen2 * mineralized_soil_N);
   } else {
      N_immobilization_demand = 0.0;
      fstd = 1.0;
   }

   mineralized_N_loss +=  N_immobilization_demand*fstd;

   /* STRUCTURAL POOL */
   /* compute associated decomposition of strucural C and lignin */ 
   structural_C_loss = A * Lc * K1 * structural_soil_C * fstd;  
   structural_L_loss = A * Lc * K1 * structural_soil_L * fstd; 

   /* fast pool */
   /* NOTE: one of the fast pool losses in Century includes a soil type 
      dependence thats not included here */ 
   fast_C_loss =  A * K2 * fast_soil_C;
   fast_N_loss =  A * K2 * fast_soil_N;
   mineralized_N_input += fast_N_loss;

   if(data->open_cycles) {
      /* deposition of N in rain */ 
      mineralized_N_input += siteptr->sdata->N_conc_in_rain 
         * siteptr->sdata->precip[data->time_period];

      /* leaching losses of C */
      fast_C_loss += data->NC_perc_coeff * perc * fast_soil_C;
      /* NOTE: CENTURY only has leaching losses from fast C&N and Min N! */ 
      /* leaching losses of N */
      fast_N_loss += data->NC_perc_coeff * perc*fast_soil_N;
      mineralized_N_loss += data->NC_perc_coeff * perc * mineralized_soil_N; 
   }


   /* slow pool */
   slow_C_input = (1.0 - r_stsc) * structural_C_loss;
   slow_C_loss =  A * K3 * slow_soil_C;
   slow_N_loss = (1.0 / data->c2n_slow) * slow_C_loss;
   mineralized_N_input += slow_N_loss;
  
   /* passive pool */
   passive_C_input = 0.0;
   passive_C_loss =  A * K4 * passive_soil_C;
   passive_N_loss = (1.0 / data->c2n_passive) * passive_C_loss;
   mineralized_N_input += passive_N_loss;

   /* C POOLS */
   dfsc = fsc_in - fast_C_loss;
   dstsc = ssc_in - structural_C_loss;
   dssc = slow_C_input - slow_C_loss;
   dpsc = passive_C_input - passive_C_loss;
   dstsl = ssl_in - structural_L_loss;

   /* loss of C to atmosphere and leaching */
   rh = r_fsc * fast_C_loss + r_stsc*structural_C_loss 
      + r_ssc*slow_C_loss + r_psc*passive_C_loss;

   /* NITROGEN POOLS */
   dfsn = fsn_in - fast_N_loss;
   dmsn = mineralized_N_input - mineralized_N_loss;
}
#endif /* ED */


#ifdef MIAMI_LU
////////////////////////////////////////////////////////////////////////////////
//! Dsdt
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void patch::Dsdt (UserData* data) {

   double fast_C_loss = 0.00;
   double structural_C_loss = 0.00;

   int time_period = data->time_period;

   /****************/
   /* OTHER FLUXES */
   /****************/

   /* calculate commonly used terms once */
   A = A_function(time_period, data);

   /*C POOLS*/
   fast_C_loss = fast_soil_C * (1.0 - exp(-data->K1 * A * data->deltat) );
   structural_C_loss = structural_soil_C * (1.0 - exp(-data->K2 * A * data->deltat) );
  
   fast_soil_C -= fast_C_loss;
   structural_soil_C -= structural_C_loss;
  
   /*loss of C to atmosphere and leaching */
   rh = fast_C_loss + structural_C_loss;
}
#endif /* MAIMI_LU */


////////////////////////////////////////////////////////////////////////////////
//! A_function
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double patch::A_function (unsigned int time_period, UserData* data) {
    
   /* The combined reduction in decomposition rate as a funciton of TEMP and MOIST */
   /* Based on CENTURY Parton et al 1993 GBC 7(4):785-809 */
   /* and Ben's copy of century code */

   /************/
   /* WARNINGS */
   /************/
   /* NOTE: This DOES include patch scale variability in soil moisture!, 
      a good thing */
  
   /* NOTE: this DOES NOT include patch scale variability in temperature 
      or soil temperature, a bad thing.*/
  
   site* currents = siteptr;

   double soil_temp;
#ifdef ED
   //double rain = currents->sdata->precip[time_period];
   //double pet = currents->sdata->pet[time_period];
   soil_temp = currents->sdata->soil_temp[time_period];
#endif
#ifdef MIAMI_LU
   soil_temp = currents->sdata->temp_average;
#endif /* MIAMI_LU */

   /* EFFECT OF TEMPERATURE   */
   /* from ben's century code */
   double Tmax = 45.0;
   double Topt = 35.0;
   double tshr = 0.2; 
   double tshl = 2.63;
   double t1 = (Tmax - soil_temp) / (Tmax - Topt);
   double t2 = exp( (tshr / tshl) * (1. - pow(t1, tshl)) );
   double Td = pow(t1, tshr) * t2; // rate multiplier due to temp 
   

   /*EFFECT OF MOISTURE*/
   double Wd;        /* rate reduction due to mositure */   
#if defined ED
   /*Linn and Doran, 1984, Soil Sci. Amer. J. 48:1267-1272*/
   /*This differs from the Century Wd*/
   /*gets rid of PET, and rainfall in Wd term*/
   if (theta <= 0.3) {
      Wd = 0.2;
   } else if (theta <= 0.6) {
      Wd = theta / 0.6;
   } else {
      Wd = 0.6 / (1.2 * theta);
   }
   return (Td * Wd);

#elif defined MIAMI_LU
   /* Crude Approximation to centruy function, ignoring stored water term. */
  
   double Rmax = 6.0; /* was 6.0, but caused problems at some sites */  
   double Ropt = 0.9;
   double Rshr = 90.; 
   double Rshl = 1.;
   double Rratio = currents->sdata->precip_average / currents->sdata->pet_average;
   if (Rratio > Rmax) Rratio = Rmax - 0.1; /* GCH Trap */
   double R1 = (Rmax - Rratio) / (Rmax - Ropt);
   double R2 = exp( (Rshr / Rshl) * (1. - pow(R1, Rshl)) );
   Wd = pow(R1, Rshr) * R2;

   double term = Td * Wd;
   if (term < 0.005) term = 0.005; /* GCH Trap */

   return term; /* the combined (multiplicative) effect of temp and water on decom rates */

#endif /* ED vs. MIAMI_LU */
}

