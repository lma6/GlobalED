#include <stdio.h>
#include "edmodels.h"
#include "site.h"
#include "patch.h"
#include "fire.h"
#include "read_site_data.h"
#include "cohort.h"
#include <cstring>

#define THETA_CRIT 0.2      /* theta threshold for leaf drop (mm)*/
#define L_FRACT 0.5 /* 0.5 fraction of leaves retained after leaf fall */

using namespace std;

bool temperature_leaf_drop(int t, site**);
bool leaf_on_date(int, site**);

////////////////////////////////////////////////////////////////////////////////
//! phenology
//! 
//!

//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void phenology (int t, patch** patchptr, UserData* data) {   
   patch* currentp = *patchptr;
   site* currents = currentp->siteptr;
   
   /*update_site_water(&currentp,data);*/
   
   update_fuel(t,&currentp,data); /* calc. fuel levels and set site fire flag */
   
   while (currentp != NULL) {
      /* go to every cohort */
      cohort* currentc = currentp->tallest;
      while (currentc != NULL) {
         size_t spp = currentc->species;
         currentc->status = 0; /* set status to indicate no leaf drop */  
         /* if drought-deciduous or cold-deciduous species    */
         /* temp=10 degrees gives good growing season pattern */
     
          //test_larch
          int if_trigger = 0;
          if((data->phenology[spp] == 1) && (currentp->theta < THETA_CRIT))
              if_trigger = 1;
          if((data->phenology[spp] == 2) && (currents->sdata->temp[data->time_period] < 10.0))
              if_trigger = 1;
          if((data->phenology[spp] == 4) && ((currentp->theta < THETA_CRIT) || (currents->sdata->temp[data->time_period] < 10.0)))
              if_trigger = 1;
          if((data->phenology[spp] == 3) && ((currentp->theta < THETA_CRIT) || (currents->sdata->temp[data->time_period] < 10.0)))
              if_trigger = 1;
          
          //test_larch
          if(!strcmp(data->title[spp],"early_succ") and currents->is_frozen_early_succ)
          {
              if_trigger = 2;
          }
          if(!strcmp(data->title[spp],"mid_succ") and currents->is_frozen_mid_succ)
          {
              if_trigger = 2;
          }
          if(!strcmp(data->title[spp],"late_succ") and currents->is_frozen_late_succ)
          {
              if_trigger = 2;
          }
          if(!strcmp(data->title[spp],"cold_decid") and currents->is_frozen_cold_decid)
          {
              if_trigger = 2;
          }
          else if (!strcmp(data->title[spp],"evergreen_short") and currents->is_frozen_evergreen_short)
          {
              if_trigger = 2;
          }
          
          
          //test_larch
//         if( ((data->phenology[spp] == 1) && (currentp->theta < THETA_CRIT))
//            || ((data->phenology[spp] == 2)
//                && (currents->sdata->temp[data->time_period] < 10.0)))  // 10.0
        if(if_trigger)
         {
            currentc->status = 5; /* set status to indicate leaf drop */
             
             //test_larch
             double leaf_litter = 0.0;
             
//            double leaf_litter = (1.0 - L_FRACT) * currentc->bl * (currentc->nindivs / currentp->area);
//            /* add to patch litter flux terms */
//            currentp->litter += leaf_litter;
//
//            /* add litter to soil pools */
//            currentp->fast_soil_N += data->fraction_balive_2_fast * leaf_litter
//                                   / data->c2n_leaf[currentc->species];
//            currentp->fast_soil_C += data->fraction_balive_2_fast * leaf_litter;
//            currentp->structural_soil_C += (1.0 - data->fraction_balive_2_fast) * leaf_litter;
//            currentp->structural_soil_L += (1.0 - data->fraction_balive_2_fast)
//                                         * (data->l2n_stem / data->c2n_stem) * leaf_litter;
//
//            /* decrement balive for leaf litterfall */
//            currentc->balive -= (1.0 - L_FRACT) * currentc->bl;
//            /* move retained leaf matter to virtual leaf pool */
//            currentc->blv = L_FRACT * currentc->bl;
//            /* reset leaf pool to zero */
//            currentc->bl = 0.0;
             
             //test_larch
             //further reduce balive from blv when bl is zero. This is particular for cold-deciduout in leaf-droped month,if frost happen, balive should continue decreasing even if has no leaves.
             //Here I assume blv acts same function as bud, so winter frose/freezing will cause bud injury even no leaves exist
             if((if_trigger==2) and (currentc->bl<1e-10))
             {
                 leaf_litter = (1.0 - L_FRACT) * currentc->blv * (currentc->nindivs / currentp->area);
                 currentp->litter += leaf_litter;
                 currentc->balive -= (1.0 - L_FRACT) * currentc->blv;
                 currentc->blv = L_FRACT * currentc->blv;
                 currentc->bl = 0.0;

                 /* add litter to soil pools */
                 currentp->fast_soil_N += data->fraction_balive_2_fast * leaf_litter
                 / data->c2n_leaf[currentc->species];
                 currentp->fast_soil_C += data->fraction_balive_2_fast * leaf_litter;
                 currentp->structural_soil_C += (1.0 - data->fraction_balive_2_fast) * leaf_litter;
                 currentp->structural_soil_L += (1.0 - data->fraction_balive_2_fast)
                 * (data->l2n_stem / data->c2n_stem) * leaf_litter;
             }
             else
             {
                leaf_litter = (1.0 - L_FRACT) * currentc->bl * (currentc->nindivs / currentp->area);
                 /* add to patch litter flux terms */
                 currentp->litter += leaf_litter;

                 /* add litter to soil pools */
                 currentp->fast_soil_N += data->fraction_balive_2_fast * leaf_litter
                 / data->c2n_leaf[currentc->species];
                 currentp->fast_soil_C += data->fraction_balive_2_fast * leaf_litter;
                 currentp->structural_soil_C += (1.0 - data->fraction_balive_2_fast) * leaf_litter;
                 currentp->structural_soil_L += (1.0 - data->fraction_balive_2_fast)
                 * (data->l2n_stem / data->c2n_stem) * leaf_litter;

                 /* decrement balive for leaf litterfall */
                 currentc->balive -= (1.0 - L_FRACT) * currentc->bl;
                 /* move retained leaf matter to virtual leaf pool */
                 currentc->blv = L_FRACT * currentc->bl;
                 /* reset leaf pool to zero */
                 currentc->bl = 0.0;
             }
             
             
             
         } /* end if */
         currentc = currentc->shorter;
      } /* end while cohort */
      currentp = currentp->older;
   } /* end while patch */
}

////////////////////////////////////////////////////////////////////////////////
//! new_phenology
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void new_phenology (int t, patch** patchptr, UserData* data) {
   
   patch* currentp = *patchptr;
   site* currents = currentp->siteptr;
   printf("HERE \n");
   
   /*update_site_water(&currentp,data);*/
   
   update_fuel(t,&currentp,data); /* calc. fuel levels and set site fire flag */
   
   int leaf_on;
   if (((t%N_CLIMATE)>0.38356*N_CLIMATE)&((t%N_CLIMATE)<0.7124*N_CLIMATE)) leaf_on = 1;
   //if (((t%N_CLIMATE)>0.3698*N_CLIMATE)&((t%N_CLIMATE)<0.7261*N_CLIMATE)) leaf_on = 1;
   //if (((t%N_CLIMATE)>0.378*N_CLIMATE)&((t%N_CLIMATE)<0.721*N_CLIMATE)) leaf_on = 1;
   else leaf_on =0;
   //leaf_on = leaf_on_date(t, &currents);
   //printf("%d, %d, %d, %f\n", t, leaf_on_date(t, &currents), leaf_on, t*1./N_CLIMATE);
   /*if (abs((t%N_CLIMATE)*1./N_CLIMATE -0.5)< 0.24) leaf_on = 1;
   else leaf_on = 0;*/
   
   while (currentp != NULL) {
      /* go to every cohort */
      cohort* currentc = currentp->tallest;
      while (currentc != NULL) {
         size_t spp = currentc->species;
         currentc->status = 0; /* set status to indicate no leaf drop */ 
         temperature_leaf_drop(t, &currents);  
         
         if( ((data->phenology[spp] == 1) && (currentp->theta < THETA_CRIT))
            || ((data->phenology[spp] == 2)
                && !leaf_on))
         {
            currentc->status = 5;
            double leaf_litter = (1.0 - L_FRACT) * currentc->bl * (currentc->nindivs / currentp->area);
            currentp->litter += leaf_litter;
            
            currentp->fast_soil_N += data->fraction_balive_2_fast * leaf_litter
                                   / data->c2n_leaf[currentc->species];
            currentp->fast_soil_C += data->fraction_balive_2_fast * leaf_litter;
            currentp->structural_soil_C += (1.0 - data->fraction_balive_2_fast) * leaf_litter;
            currentp->structural_soil_L += (1.0 - data->fraction_balive_2_fast)
                                         * (data->l2n_stem / data->c2n_stem) * leaf_litter;
               
            currentc->balive -= (1.0 - L_FRACT) * currentc->bl;
            currentc->blv = L_FRACT * currentc->bl;
            currentc->bl = 0.0;
         } 
         currentc = currentc->shorter;
      } /* end while cohort */
      currentp = currentp->older;
   } /* end while patch */
}

////////////////////////////////////////////////////////////////////////////////
//! temperature_leaf_drop
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
bool temperature_leaf_drop(int t, site** siteptr){
   
   site* currents = *siteptr;
   double average = 0.0;
   for (size_t i=0;i<24;i++) {
      average += currents->ave_temp[i]/24.;
   }
   return (average<10.);
}
#if FTS
// Not being used right now
////////////////////////////////////////////////////////////////////////////////
//! leaf_on_date
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
bool leaf_on_date(int t, site** siteptr){
   /*Decides what day of year for leaf on/leaf off using 30 day running average*/
   /*ONLY FOR NORTHERN HEMISPHERE AT THE MOMENT!!!!!*/
   site* currents = *siteptr;
   double period_average[N_CLIMATE_INPUT], running_average;
   double leaf_on=-1, leaf_off=-1;
   
   
   for (size_t period=0; period<N_CLIMATE_INPUT; period++){
      period_average[period] = 0.;
      for (size_t interval=0; interval<CLIMATE_INPUT_INTERVALS; interval++) {
         period_average[period] += currents->sdata->Input_Temperature[interval+period*CLIMATE_INPUT_INTERVALS]/CLIMATE_INPUT_INTERVALS;
      }
   }
   for (size_t period=0; period<N_CLIMATE_INPUT; period++) {
      running_average = 0;
      for (size_t delta_period=-15; delta_period<15; delta_period++) {
         running_average+=period_average[(period+delta_period+N_CLIMATE_INPUT)%N_CLIMATE_INPUT];
      }
      running_average = running_average/30;
      if (leaf_on==-1&&running_average>10) leaf_on = period;
      if (running_average>10) leaf_off = period;
   }
   if (t%N_CLIMATE<leaf_on*1./N_CLIMATE_INPUT*N_CLIMATE-0.5) return false;
   else if (t%N_CLIMATE_INPUT<leaf_off*1./N_CLIMATE_INPUT*N_CLIMATE-0.5) return true;
   else return false;
}
#endif

/******************************************************************************/
