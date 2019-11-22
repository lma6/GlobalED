#include <cmath>
#include <cstring>
#include "edmodels.h"
#include "site.h"
#include "read_site_data.h"
#include "patch.h"
#include "cohort.h"

////////////////////////////////////////////////////////////////////////////////
//! survivorship_from_disturbance
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double cohort::survivorship_from_disturbance ( int q, UserData* data ) {
   /*q is the disturbance track*/

   /*printf("sfd: entering...\n");*/
   double csurv = 0.0;

   /*treefalls*/    
   if ( q == 0 ) { 
      if ( hite < data->treefall_hite_threshold ) {
         csurv = data->treefall_s_ltht[species]; 
         /*csurv = data->treefall_s_gtht[spp] + data->treefall_s_ltht[spp] * 
            ( 1.0 - cc->hite/data->treefall_hite_threshold );*/
      } else  {
         csurv = data->treefall_s_gtht[species];      
      }
   }  

   /*fires*/    
   else if ( q == 1 ) { 

      if ( data->is_grass[species] ) {
          csurv = 1.00;
      } else {
          csurv = 0.0;
          //test_larch
          //topkill rate from Hoffmann et al ., 2009
//          double base_csurv = 0.0;
          csurv = 0.3;
//          if (species==6)
//              csurv = 0.7;
//          else
//          {
//              base_csurv = 0.5;
//              if(dbh>20.0)
//                  csurv = base_csurv;
//              else
//                  csurv = base_csurv * dbh/20.0;
//          }
//
          
      }

   } /* end if on fires */
  
   /*printf("sfd: exiting...\n");*/
   return csurv;
}

////////////////////////////////////////////////////////////////////////////////
//! cohort_modifications_from_disturbance
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void cohort_modifications_from_disturbance ( int q, cohort** pcurrentc, patch** newp,
                                             UserData* data ) {   
   cohort* cc = *pcurrentc;
    ///CarbonConserve
    patch* np = *newp;
   /* q is disturbance track */

   /* grass species lose all leaf biomass and are returned to HMIN in Fires */

   if ( q == 1 ) {
      if ( data->is_grass[cc->species] ) {
          ///CarbonConserve
          double old_balive = cc->balive, old_bdead = cc->bdead;
         cc->hite = data->hgt_min[cc->species];
         cc->dbh = cc->Dbh(data);     
         cc->bdead = cc->Bdead(data);
         cc->balive = cc->balive - cc->bl;
         cc->bl = 0.0;
         cc->b = cc->balive + cc->bdead;
          
          ///CarbonConserve
          np->fire_emission += (old_balive+old_bdead-cc->balive-cc->bdead)*cc->nindivs*LANDUSE_FREQ;
      }
   }

   /*NOTE, CARBON LEAK IN CLOSED MODE BECAUSE BIOMASS CHANGE NOT YET LOADED INTO
     LITTER, GOES TO ATM MBY IMPLICATION*/
    /* NORE, THERE WILL NO CARBON LEAK IF CHAGE OF BIOMASS IS PUT INTO fire_emission OF PATCH IN CarbonConserve Modification by Lei */
}

////////////////////////////////////////////////////////////////////////////////
//! den_dep_death
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double cohort::den_dep_death (UserData* data ) {
   double dmort = 0.0;

//   // Mortality is computed differently for temperate vs tropical species
//   if (data->is_tropical[species] or data->is_grass[species]) {
//      dmort = data->m1 * ( data->rho_max1 - data->rho[species] );
//   } else if (!data->is_tropical[species]){
//      dmort = data->m1 * ( data->rho_max2 - data->rho[species] );
//   }
    
    //test_larch
//    if(siteptr->is_tropical_site==1)
//    {
//        dmort = data->m1 * ( data->rho_max1 - data->rho[species] );
//    }
//    else
//    {
//        if (!strcmp(data->title[species],"early_succ"))
//            dmort += 0.067;
//        else if (!strcmp(data->title[species],"mid_succ"))
//            dmort += 0.059;
//        else if (!strcmp(data->title[species],"late_succ"))
//            dmort += 0.011;
//
//    }
    
    //test_larch
    data->m1 = 0.15;
    dmort = data->m1 * ( data->rho_max1 - data->rho[species] );
    dmort += 0.011;
//    if(siteptr->is_tropical_site)
//        dmort += 0.02;
//    if(siteptr->climate_zone==1)
//        dmort += 0.02;
    
    //test_mor
//    dmort = data->den_ind_mort[species];
//    dmort = 1.0;
//    if (species==0)
//        dmort = 0.0;  // 0.006144 used in ED2
//    if (species==1)
//        dmort = 0.0;  // 0.00428 used in ED2
//    if (species==2)
//        dmort = 0.03;  // 0.006144 used in ED2
//    if (species==3)
//        dmort = 0.02;  // 0.00428 used in ED2
//    if (species==4)
//        dmort = 0.002; // 0.0033928 used in ED2
//    if (species==5)
//        dmort = 0.002;  // 0.002
//    if (species==6)
//        dmort = 0.002; // 0.002
    

   double lat = siteptr->sdata->lat_;
   if ( hite < data->treefall_hite_threshold ) {  
      if ( (lat < data->tropic_n_limit) && (lat > data->tropic_s_limit) ) {
         dmort += data->treefall_max_disturbance_rate_trop; 
      } else {
         dmort += data->treefall_max_disturbance_rate_temp; 
      }
   } 
   
   dmort += data->m2 * 1.0 / ( 1.0 + exp( data->m3 * cbr_bar ));
   /*dmort += data->m2 / ( 1.0 + exp( data->m3 * cbr_bar[data->time_period] ) );*/

   /* do mean field treefall disturbance when fire *
    * is the dominant disturbance mode             */
   if ( patchptr->treefall_as_dndt_flag == 1 ) {
      if ( (lat < data->tropic_n_limit) && (lat > data->tropic_s_limit) ) {
         dmort += data->treefall_max_disturbance_rate_trop; 
      } else {
         dmort += data->treefall_max_disturbance_rate_temp; 
      }
   }
   /* do mean field fire when treefall is dominant mode */
   if ( patchptr->fire_as_dndt_flag == 1 ) {
      dmort += patchptr->fire_dndt_factor 
             * ( 1.0 - survivorship_from_disturbance(1, data) );
   }

   if(data->additional_mort and ( !data->mort_s_hemi or 
           (lat > data->tropic_s_limit) )) {
      dmort += data->den_ind_mort[species];
   } else if (data->additional_mort and data->mort_s_hemi 
           and (lat < data->tropic_s_limit) and data->is_cold_deciduous[species]) {
      // Only applicable to Southern Hemisphere, Cold Deciduous PFT
      dmort += data->den_ind_mort_s_hemi[species];
   }

   return dmort;
}


/******************************************************************************/
