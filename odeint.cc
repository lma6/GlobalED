/*5/1/00 modified to never crash by PRM*/
/*4/May/99 modified to handle meterology*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <cmath>

#include "edmodels.h"
#include "site.h"
#include "patch.h"
#include "read_site_data.h"
#include "cohort.h"

static void f(double t, void *f_data);

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//! cm_sodeint
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
int cm_sodeint (patch** patchptr, int time_step, double t1, double t2, 
                UserData* data) {

   /* Now using RK2 method. */
   /* Split step in half and re-integrate if one value is negative */
   
   //Return 0 if integrated successfully, 1 if encountered problem
  
   int iout, flag;
   double deltat;
   int split_factor, split_count;
    

   patch* currentp = *patchptr;
    //checkstep
    currentp->rh_avg = 0.0;
    currentp->gpp_avg = 0.0;
    currentp->npp_avg = 0.0;
    currentp->fire_emission = 0.0;
    
#if LANDUSE
    currentp->forest_harvested_c = 0.0;
    currentp->past_harvested_c = 0.0;
    currentp->crop_harvested_c = 0.0;
#endif
    cohort* currentc = currentp->shortest;
    while (currentc != NULL) {
        currentc->p[0] = 0.0;
        currentc->p[1]  = 0.0;
        currentc->p_avg[0]  = 0.0;
        currentc->p_avg[1]  = 0.0;
        currentc->gpp_avg = 0.0;
        currentc->npp_avg = 0.0;
        currentc->md_avg = 0.0;
        currentc = currentc->taller;
    }
   iout = 0;
   while (iout<data->substeps){ 
      split_factor=1;
      split_count = 0;
      int iout2 = 0;
      while (iout2 <split_factor){
         currentp->save_old();
         f(t1, currentp);
         currentp->copy_derivatives();
          
#if CHECK_C_CONSERVE
          ///CarbonConserve
          cohort* mlcc= NULL;
          double all_tb_before=0.0, all_sc_before=0.0, all_tc_before=0.0, all_fire_emission_before = 0.0;
          double all_tb_after=0.0, all_sc_after = 0.0, all_tc_after=0.0, all_repro=0.0, all_npp_after = 0.0, all_rh_after = 0.0, all_fire_emission_after = 0.0;
          double actual_dt_tc = 0.0, esti_dt_tc = 0.0;
          mlcc = currentp->shortest;
          while (mlcc!=NULL) {
              all_tb_before += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
              mlcc = mlcc->taller;
          }
          all_fire_emission_before = currentp->fire_emission;
          all_sc_before = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
          all_tc_before = all_tb_before + all_sc_before;
#endif
          
         do {//Repeat until compare_derivatives returns true, halving dt each time
            deltat = data->deltat * 1. / (data->substeps*split_factor);
            currentp->Update_Water(t1, currentp->siteptr->data, deltat/2.);
            currentp->fast_soil_C +=currentp->dfsc * deltat / 2.;
            currentp->structural_soil_C += currentp->dstsc * deltat / 2.;
            currentp->slow_soil_C +=  currentp->dssc * deltat / 2.;
            currentp->mineralized_soil_N += currentp->dmsn * deltat / 2.;
            currentp->fast_soil_N += currentp->dfsn * deltat / 2.;
            currentp->passive_soil_C += currentp->dpsc * deltat / 2.;
            currentp->structural_soil_L += currentp->dstsl * deltat / 2.;
            cohort* currentc = currentp->shortest;
            while (currentc != NULL) {
               currentc->nindivs += currentc->dndt * deltat / 2.;
               currentc->dbh += currentc->ddbhdt * deltat / 2.;
               currentc->balive += currentc->dbalivedt * deltat / 2.;
               currentc->bdead += currentc->dbdeaddt * deltat / 2.;
               currentc = currentc->taller;
            }
            f(t1 + deltat / 2., currentp);
            flag = currentp->check_for_negatives(deltat);
            if (flag==0) {
               break;
            }
            else if ((flag==2)|(split_count>15)) {
                //if (flag==2) printf("Failed for flag=2\n");
                //if (split_count>15) printf("Failed for split_count>15\n");
               return 1;
            }
            //Otherwise split step
            split_factor*=2;
            split_count++;
            iout2*=2;
            currentp->load_old();
            currentp->load_derivatives();
         } while (true);
          
          ///CarbonConserve -- Lei Ma
          /// Here, the oder of updating soil carbon pool and cohort attribures are reversed because Litter() should be called after update of all cohorts.
          /// Because the above do_while for loop just check whether deriative in current substep are valid. If yes, then update cohorts using the checked deriatives. If not, make more substep.
          /// However, calculation of litter should use final cohort density after deltat rather than ones in deltat/2.0. Therefore, we need to update cohorts first, then call Listter() and Dsdt() again.
          
          /// This block is original
//         currentp->Update_Water(t1, currentp->siteptr->data, deltat/2.);
//         currentp->fast_soil_C = currentp->old_fast_soil_C + currentp->dfsc * deltat;
//         currentp->structural_soil_C = currentp->old_structural_soil_C + currentp->dstsc * deltat;
//         currentp->slow_soil_C = currentp->old_slow_soil_C + currentp->dssc * deltat;
//         currentp->mineralized_soil_N = currentp->old_mineralized_soil_N + currentp->dmsn * deltat;
//         currentp->fast_soil_N = currentp->old_fast_soil_N + currentp->dfsn * deltat;
//         currentp->passive_soil_C = currentp->old_passive_soil_C + currentp->dpsc * deltat;
//         currentp->structural_soil_L = currentp->old_structural_soil_L + currentp->dstsl * deltat;
          /// The above is original
          
         cohort* currentc = currentp->shortest;
          double tmp_gpp_avg = 0.0, tmp_npp_avg = 0.0, tmp_repro_avg = 0.0;
         while (currentc != NULL) {
            currentc->nindivs = currentc->old_nindivs + currentc->dndt * deltat;
            currentc->dbh = currentc->old_dbh + currentc->ddbhdt * deltat;
            currentc->balive = currentc->old_balive + currentc->dbalivedt * deltat; 
            currentc->bdead = currentc->old_bdead + currentc->dbdeaddt * deltat;
             ///CarbonConserve
             currentc->Allocate_Biomass(data);
             //checkstep
             //currentc->gpp_avg += currentc->gpp*1./(data->substeps*split_factor);
             //currentc->npp_avg += currentc->npp*1./(data->substeps*split_factor);
             //currentc->md_avg += currentc->md*1./(data->substeps*split_factor);
             /// Here, cohorts die before photosynthesis, there updated cohort density is used than old_nindivs
             tmp_gpp_avg += currentc->gpp * currentc->nindivs/currentp->area;
             tmp_npp_avg += currentc->npp * currentc->nindivs/currentp->area;
             tmp_repro_avg += currentc->p[0] * currentc->nindivs/currentp->area;
             currentc->p_avg[0] += currentc->p[0]*currentc->nindivs*1./(data->substeps*split_factor);
             currentc->p_avg[1] += currentc->p[1]*currentc->nindivs*1./(data->substeps*split_factor);
            currentc = currentc->taller;
         }
          /// This block is added by Lei, should delelte if it does not work
          ///CarbonConserve
          currentp->Litter(t1+deltat, data);
          currentp->Dsdt(data->time_period, t1+deltat, data);
          currentp->Update_Water(t1, currentp->siteptr->data, deltat/2.);
          currentp->fast_soil_C = currentp->old_fast_soil_C + currentp->dfsc * deltat;
          currentp->structural_soil_C = currentp->old_structural_soil_C + currentp->dstsc * deltat;
          currentp->slow_soil_C = currentp->old_slow_soil_C + currentp->dssc * deltat;
          currentp->mineralized_soil_N = currentp->old_mineralized_soil_N + currentp->dmsn * deltat;
          currentp->fast_soil_N = currentp->old_fast_soil_N + currentp->dfsn * deltat;
          currentp->passive_soil_C = currentp->old_passive_soil_C + currentp->dpsc * deltat;
          currentp->structural_soil_L = currentp->old_structural_soil_L + currentp->dstsl * deltat;
          /// The above is added by Lei.
          
          //CarbonConserve
          currentp->gpp_avg +=tmp_gpp_avg*1./(data->substeps*split_factor);
          currentp->npp_avg +=tmp_npp_avg*1./(data->substeps*split_factor);
          currentp->rh_avg +=currentp->rh*1./(data->substeps*split_factor);
          currentp->fire_emission += currentp->fire_c_loss*1./(data->substeps*split_factor);
//          total_litter += currentp->litter *1./(data->substeps*split_factor);
//          total_repro += tmp_repro_avg*1./(data->substeps*split_factor);
         iout2++;

#if CHECK_C_CONSERVE
          ///CarbonConserve
          mlcc = currentp->shortest;
          while (mlcc != NULL) {
              all_tb_after += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
              all_npp_after += mlcc->npp*mlcc->nindivs/currentp->area;
              all_repro += (mlcc->p[0]+mlcc->p[1])*mlcc->nindivs/currentp->area;
              mlcc = mlcc->taller;
          }
          all_fire_emission_after = currentp->fire_emission;
          all_rh_after = currentp->rh;
          all_sc_after = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
          all_tc_after = all_tb_after + all_sc_after;
          actual_dt_tc = all_tc_after - all_tc_before;
          esti_dt_tc = (all_npp_after-all_rh_after-all_repro*(1-data->sd_mort)-currentp->fire_c_loss)*deltat;
          
          if (abs(actual_dt_tc - esti_dt_tc)>1e-9)
          {
              printf("Carbon leakage in substep_integ%d : imbalance    %.15f actual_dt_tc %.15f esti_dt_tc  %.15f\n",iout,actual_dt_tc-esti_dt_tc,actual_dt_tc,esti_dt_tc);
              printf("                                 : patch_tc_bf  %.15f patch_sc_bf  %.15f patch_tb_bf %.15f\n",all_tc_before,all_sc_before,all_tb_before);
              printf("                                 : patch_tc_af  %.15f patch_sc_af  %.15f patch_tb_af %.15f patch_fire_emi_bf %.15f\n",all_tc_after,all_sc_after,all_tb_after,all_fire_emission_before);
              printf("                                 : patch_npp_af %.15f patch_rh_af  %.15f patch_repro %.15f patch_fire_emi_af %.15f\n",all_npp_after,all_rh_after,all_repro,all_fire_emission_after);
              printf("                                    : site_lat %.3f site_lon %.3f area %.6f\n",currentp->siteptr->sdata->lat_,currentp->siteptr->sdata->lon_,currentp->area);
              printf(" --------------------------------------------------------------------------------------\n");
          }
#endif
      }
      iout++;
   }
    currentc =  currentp->shortest;
    double test_repro = 0.0;
    while (currentc!=NULL) {
        currentc->p[0] = currentc->p_avg[0]/currentc->nindivs;
        currentc->p[1] = currentc->p_avg[1]/currentc->nindivs;
        test_repro +=(currentc->p_avg[0]+currentc->p_avg[1])/currentp->area;
        currentc = currentc->taller;
    }
   return 0;   /* return to community dynamics */
}

////////////////////////////////////////////////////////////////////////////////
//! f
//! Functions Called by the CVODE Solver. f routine. Compute f(t,u).
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
static void f(double t, void *f_data){

   /* Extract constants from data structure */
   patch* currentp = (patch*)(f_data);    
   site* currents = currentp->siteptr;
   UserData* data = currents->data;
   
   currents->function_calls++;
   
   cohort* currentc = currentp->shortest;
   // The order of functions is critical. Important to call Allocate_Biomass first 
   // so that functions following it are using accurate biomass/allometry estimates
   while(currentc!=NULL){
      currentc->Allocate_Biomass(data);
      currentc=currentc->taller;
   }   

   if(data->stiff_light) {
      sort_cohorts(&currentp,data);
      light_levels(&currentp,data);
   }

   currentp->Water_and_Nitrogen_Uptake(data->time_period,t,data);
   currentp->Dwdt(t,data); 

#if DOES_COMPILE
   if(!data->patch_dynamics) /* if no pd update distrates to be applied as mort */
      calculate_patch_disturbance_rates(t,&currentp,data); /* update dist rates */ 
#endif

   /* call growth function */
   /* printf("f: calculating growth function \n"); */
   currentc = currentp->shortest;
   while(currentc!=NULL){
      /* printf("f: currentc %p \n",currentc); */
      currentc->Growth_Derivatives(t,data);  
      currentc=currentc->taller;
   }

   currentp->Litter(t,data);
   currentp->Dsdt(data->time_period,t,data);
   return;
}


////////////////////////////////////////////////////////////////////////////////
//! Water_and_Nitrogen_Uptake
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void patch::Water_and_Nitrogen_Uptake (unsigned int time_period, double time, UserData* data) {

   site* currents= siteptr;
   
   /*TRANSPIRATION****************************/
   cohort* currentc=shortest;
   while(currentc != NULL){
      size_t spp = currentc->species;
      size_t pt = currentc->pt;
      currentc->get_cohort_vm0(data);
      currentc->Vm0_bin= currentc->get_cohort_vm0_bin(currentc->Vm0, data);
      currentc->leaf_area = currentc->bl*data->specific_leaf_area[spp];
      
      size_t time_index = data->time_period;
      size_t light_index = 0; /*calc light index*/
       while(currents->sdata->light_levels[spp][light_index] > currentc->lite){
           light_index++;
       }

       if (currents->sdata->E[spp][time_index][light_index]<-1000)
       {
           currents->sdata->compute_mech(pt,spp,currentc->Vm0,currentc->Vm0_bin,time_index,light_index,data);
       }
       if (currents->sdata->Eb[spp][time_index][N_LIGHT-1]<-1000)
       {
           currents->sdata->compute_mech(pt,spp,currentc->Vm0,currentc->Vm0_bin,time_index,N_LIGHT-1,data);
       }
       currentc->E_pot = currents->sdata->E[spp][time_index][light_index];
       currentc->E_pot *= currentc->leaf_area*N_CLIMATE/1000.0;
       currentc->E_shut = currents->sdata->Eb[spp][time_index][N_LIGHT-1];
       currentc->E_shut *= currentc->leaf_area*N_CLIMATE/1000.0;
 
      /*NITROGEN UPTAKE**************************/
      /* set fs_open to 1 and */
      /* calc implied nitrogen uptake (kgN/yr) per plant*/
      currentc->fs_open=1.0;
      currentc->N_uptake_pot= currentc->nitrogen_demand_function(time,data);
      /*convert to kgN/m2/yr per plant*/
      currentc->N_uptake_pot /= currentc->patchptr->area;

      /*nitrogen uptake with stomates shut*/
      currentc->N_uptake_shut = 0.0;

      if(data->water_competition) {
        double water_supply;
        double wilt_factor;
        water_supply=data->water1*currentc->br*water*data->mass_of_water;        
        currentc->fsw = (water_supply-currentc->E_shut)/(currentc->E_pot + water_supply-currentc->E_shut);
        if (currentc->fsw<0) {
           //If demand>supply, stomates shut and leaves wilt to the extent that supply := demand
           //printf("1: %f, %f\n", water_supply, currentc->E_shut);
           currentc->fsw = 0.; 
           wilt_factor = water_supply*1./currentc->E_shut;
           currentc->blv = currentc->blv + (1-wilt_factor)*currentc->bl;
           currentc->bl *= wilt_factor;
           currentc->leaf_area = currentc->bl*data->specific_leaf_area[spp];
        
            if (currents->sdata->E[spp][time_index][light_index]<-1000)
            {
                currents->sdata->compute_mech(pt,spp,currentc->Vm0,currentc->Vm0_bin,time_index,light_index,data);
            }
            if (currents->sdata->Eb[spp][time_index][N_LIGHT-1]<-1000)
            {
                currents->sdata->compute_mech(pt,spp,currentc->Vm0,currentc->Vm0_bin,time_index,N_LIGHT-1,data);
            }
            currentc->E_pot = currents->sdata->E[spp][time_index][light_index];
            currentc->E_shut = currents->sdata->Eb[spp][time_index][N_LIGHT-1];

           currentc->E_pot *= currentc->leaf_area*12/1000.0;
           currentc->E_shut *= currentc->leaf_area*12/1000.0;
        } //Water_supply<E_shut

      }
      else
        currentc->fsw = 1.0;

      if(data->n_competition) {
        double nitrogen_supply;
        nitrogen_supply=data->nitrogen1*currentc->br*mineralized_soil_N;
        if(currentc->N_uptake_pot >0.00)
           currentc->fsn = nitrogen_supply/(currentc->N_uptake_pot+nitrogen_supply);
        else
           currentc->fsn = 1.00;
           /*so that fsn !> 1.0 when plant n demand is < 0*/
        }
      else
         currentc->fsn = 1.0;

      currentc->fs_open = currentc->fsw*currentc->fsn;
     
      /*printf("wnu1: fs_open %f\n",currentc->fs_open);*/

      currentc->water_uptake = currentc->fs_open*currentc->E_pot + (1.0-currentc->fs_open)*currentc->E_shut; 
      /*recompute n uptake with reduced npp, not simple function like evap
        because allocation in npp dependent*/
      currentc->nitrogen_uptake = currentc->nitrogen_demand_function(time,data);
      currentc->nitrogen_uptake/=currentc->patchptr->area;      
     
      /*printf("wnu2: fs_open %f\n",currentc->fs_open);*/
     
      /*printf("%f\n",currentc->nitrogen_uptake);*/
     
      currentc=currentc->taller;
   } /* end loop over cohorts */
  
   /* printf("wnu: calculating total uptake rates \n"); */
   total_water_uptake=0.0;
   total_water_demand=0.0;
   total_plant_nitrogen_uptake=0.0;  
   currentc=shortest;
   while(currentc != NULL){
      if(data->water_competition) {
        total_water_demand += currentc->E_pot*currentc->nindivs;
        total_water_uptake += currentc->water_uptake*currentc->nindivs;
      }
      if(data->n_competition)
         total_plant_nitrogen_uptake += currentc->nitrogen_uptake*currentc->nindivs;

      currentc=currentc->taller;
   }
   return;
  
}
////////////////////////////////////////////////////////////////////////////////
//! check_for_negatives
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
int patch::check_for_negatives(double dt){
   //Return: 0 if no problems detected
   //        1 if negative/impending negative (fixed by smaller dt)
   //        2 if Nan - bug somewhere, stop modelling site
   if ((dwdt+dfsc+dstsc+dssc+dmsn+dfsn+dpsc+dstsl)*0!=0)
   {
       return 2;
   }

   if ((water<0)|(fast_soil_C<0)|(structural_soil_C<0)|(slow_soil_C<0)|(mineralized_soil_N<0)|(fast_soil_N<0)|(structural_soil_L<0))
   {
#if 0
       if (water<0) printf("Neg-wt ");
       if (fast_soil_C<0) printf("Neg-fsc ");
       if (structural_soil_C<0) printf("Neg-stsc ");
       if (slow_soil_C<0) printf("Neg-slsc ");
       if (mineralized_soil_N<0) printf("Neg-msn ");
       if (fast_soil_N<0) printf("Neg-fsn ");
       if (structural_soil_L<0) printf("Neg-stsl ");
#endif
       return 1;
   }
   
   if ((old_water+dwdt*dt<0)|(old_fast_soil_C+dfsc*dt<0)|(old_structural_soil_C+dstsc*dt<0)|(old_slow_soil_C+dssc*dt<0)|(old_mineralized_soil_N+dmsn*dt<0)|(old_fast_soil_N+dfsn*dt<0)|(old_passive_soil_C+dpsc*dt<0)|(old_structural_soil_L+dstsl*dt<0))
   {
#if 0
       if (old_water+dwdt*dt<0) printf("Neg-wtder %f %f\n",old_water,dwdt);
       if (old_fast_soil_C+dfsc*dt<0) printf("Neg-fscder ");
       if (old_structural_soil_C+dstsc*dt<0) printf("Neg-stscder ");
       if (old_slow_soil_C+dssc*dt<0) printf("Neg-slscder ");
       if (old_mineralized_soil_N+dmsn*dt<0) printf("Neg-msnder ");
       if (old_fast_soil_N+dfsn*dt<0) printf("Neg-fsnder ");
       if (old_passive_soil_C+dpsc*dt<0) printf("Neg-pscder ");
       if (old_structural_soil_L+dstsl*dt<0) printf("Neg-stslder ");
#endif
       return 1;
   }
   
   if (water+dwdt*dt<0)
   {
#if 0
       printf("Neg-wtgerC ");
#endif
       return 1; //update_water uses old_water, not water, so check for negatives there as well.
   }
   
   cohort* cc = shortest;
   while (cc != NULL) {
      if ((cc->dndt+cc->ddbhdt+cc->dbalivedt+cc->dbdeaddt)*0!=0) return 2;
      if ((cc->nindivs<0)|(cc->dbh<0)|(cc->balive<0)|(cc->bdead<0))
      {
#if 0
          if (cc->nindivs<0) printf("Neg-nind ");
          if (cc->dbh<0) printf("Neg-dbh ");
          if (cc->balive<0) printf("Neg-balive ");
          if (cc->bdead<0) printf("Neg-bdead ");
#endif
          return 1;
      }
      if ((cc->nindivs+cc->dndt*dt<0)|(cc->dbh+cc->ddbhdt*dt<0)|(cc->balive+cc->dbalivedt*dt<0)|(cc->bdead+cc->dbdeaddt*dt<0))
      {
#if 0
          if (cc->nindivs+cc->dndt*dt<0) printf("Neg-nindder ");
          if (cc->dbh+cc->ddbhdt*dt<0) printf("Neg-dbhder ");
          if (cc->balive+cc->dbalivedt*dt<0) printf("Neg-baliveder ");
          if (cc->bdead+cc->dbdeaddt*dt<0) printf("Neg-bdeadder ");
#endif
          return 1;
      }
      cc = cc->taller;
   }
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! compare_derivatives
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
bool patch::compare_derivatives(double dt){
   //printf("%f\n", dt);
   // printf("%d, %d, %d, %d, %d, %d, %d\n", -dwdt*dt/old_water>0.2, -dfsc*dt/old_fast_soil_C>0.2, -dstsc*dt/old_structural_soil_C>0.2, -dssc*dt/old_slow_soil_C>0.2, -dmsn*dt/old_mineralized_soil_N>0.2, -dfsn*dt/old_fast_soil_N>0.2, -dstsl*dt/old_structural_soil_L>0.2);

   //Check if any nans
   if ((dwdt+dfsc+dstsc+dssc+dmsn+dfsn+dpsc+dstsl)*0!=0) {return false;}
   
   /*If upcoming time step is going to reduce quantity by more than 20%, halve timestep*/
   if ((-dwdt*dt/old_water>0.2)|
       (-dfsc*dt/old_fast_soil_C>0.2)|
       (-dstsc*dt/old_structural_soil_C>0.2)|
       (-dssc*dt/old_slow_soil_C>0.2)|
       (-dmsn*dt/old_mineralized_soil_N>0.2)|
       (-dfsn*dt/old_fast_soil_N>0.2)|
       //(-dpsc*dt/old_passive_soil_C>0.2)|
       (-dstsl*dt/old_structural_soil_L>0.2)) {return false;}
   /*Else if difference in derivatives is more than 5% of quantity or more than 50% of derivatives, halve timestep*/
   if ((abs(dwdt1 - dwdt)>max(0.5*abs(dwdt), 0.05*water/dt))|
       (abs(dfsc1 - dfsc)>max(0.5*abs(dfsc), 0.05*fast_soil_C/dt))|
       (abs(dstsc1-dstsc)>max(0.5*abs(dstsc),0.05*structural_soil_C/dt))|
       (abs(dssc1 - dssc)>max(0.5*abs(dssc), 0.05*slow_soil_C/dt))|
       (abs(dmsn1 - dmsn)>max(0.5*abs(dmsn), 0.05*mineralized_soil_N/dt))|
       (abs(dfsn1 - dfsn)>max(0.5*abs(dfsn), 0.05*fast_soil_N/dt))|
       //(abs(dpsc1 - dpsc)>max(0.1*abs(dpsc), 0.001*passive_soil_C/dt))|
       (abs(dstsl1-dstsl)>max(0.5*abs(dstsl),0.005*structural_soil_L/dt))) {return false;}
   
   cohort* cc = shortest;
   while (cc != NULL) {
      if ((cc->dndt+cc->ddbhdt+cc->dbalivedt+cc->dbdeaddt)*0!=0) {return false;}
      if ((-cc->dndt*dt/cc->nindivs>0.2)|
          (-cc->ddbhdt*dt/cc->dbh>0.2)|
          (-cc->dbalivedt*dt/cc->balive>0.2)|
          (-cc->dbdeaddt*dt/cc->bdead>0.2)) {return false;}
      if ((abs(cc->dndt1-cc->dndt)>max(0.5*abs(cc->dndt), 0.05*cc->nindivs/dt))|
          (abs(cc->ddbhdt1-cc->ddbhdt)>max(0.5*abs(cc->ddbhdt), 0.05*cc->dbh/dt))|
          (abs(cc->dbalivedt1-cc->dbalivedt)>max(0.5*abs(cc->dbalivedt), 0.05*cc->balive/dt))|
          (abs(cc->dbdeaddt1-cc->dbdeaddt)>max(0.5*abs(cc->dbdeaddt), 0.05*cc->bdead/dt))) {
         //printf("%f, %f, %f, %f\n", cc->dbalivedt, cc->dbalivedt1, cc->balive, dt);
         //printf("%d, %d, %d, %d\n", abs(cc->dndt1-cc->dndt)>max(0.1*abs(cc->dndt), 0.001*cc->nindivs/dt), abs(cc->ddbhdt1-cc->ddbhdt)>max(0.1*abs(cc->ddbhdt), 0.001*cc->dbh/dt), abs(cc->dbalivedt1-cc->dbalivedt)>max(0.1*abs(cc->dbalivedt), 0.001*cc->balive/dt), abs(cc->dbdeaddt1-cc->dbdeaddt)>max(0.1*abs(cc->dbdeaddt), 0.001*cc->bdead/dt));
         return false;  
      }
      
      cc = cc->taller;
   }
   return true;
}

////////////////////////////////////////////////////////////////////////////////
//! save_old
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void patch::save_old(){
   
   old_water = water; 
   old_fast_soil_C = fast_soil_C; 
   old_structural_soil_C = structural_soil_C;
   old_slow_soil_C = slow_soil_C; 
   old_mineralized_soil_N = mineralized_soil_N; 
   old_fast_soil_N = fast_soil_N;
   old_passive_soil_C = passive_soil_C; 
   old_structural_soil_L = structural_soil_L; 
   
   cohort* cc = shortest;
   while (cc != NULL) {
      cc->old_nindivs = cc->nindivs; 
      cc->old_dbh = cc->dbh; 
      cc->old_balive = cc->balive;
      cc->old_bdead = cc->bdead;
      cc = cc->taller;
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
//! load_old
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void patch::load_old(){
   
   water = old_water; 
   fast_soil_C = old_fast_soil_C; 
   structural_soil_C = old_structural_soil_C;
   slow_soil_C = old_slow_soil_C; 
   mineralized_soil_N = old_mineralized_soil_N; 
   fast_soil_N = old_fast_soil_N;
   passive_soil_C = old_passive_soil_C; 
   structural_soil_L = old_structural_soil_L; 
   
   cohort* cc = shortest;
   while (cc != NULL) {
      cc->nindivs = cc->old_nindivs; 
      cc->dbh = cc->old_dbh; 
      cc->balive = cc->old_balive;
      cc->bdead = cc->old_bdead;
      cc = cc->taller;
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
//! copy_derivatives
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void patch::copy_derivatives(){
   
   dwdt1 = dwdt;
   dfsc1 = dfsc;
   dstsc1 = dstsc;
   dssc1 = dssc;
   dmsn1 = dmsn;
   dfsn1 = dfsn;
   dpsc1 = dpsc;
   dstsl1 = dstsl;
   
   cohort* cc = shortest;
   while (cc!=NULL){
      cc->dndt1 = cc->dndt;
      cc->ddbhdt1 = cc->ddbhdt;
      cc->dbalivedt1 = cc->dbalivedt;
      cc->dbdeaddt1 = cc->dbdeaddt;
      cc = cc->taller;
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
//! load_derivatives
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void patch::load_derivatives(){
   
   dwdt = dwdt1;
   dfsc = dfsc1;
   dstsc = dstsc1;
   dssc = dssc1;
   dmsn = dmsn1;
   dfsn = dfsn1;
   dpsc = dpsc1;
   dstsl = dstsl1;
   
   cohort* cc = shortest;
   while (cc!=NULL){
      cc->dndt = cc->dndt1;
      cc->ddbhdt = cc->ddbhdt1;
      cc->dbalivedt = cc->dbalivedt1;
      cc->dbdeaddt = cc->dbdeaddt1;
      cc = cc->taller;
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
//! check_quantities
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
int patch::check_quantities(){
   //printf("%f, %f, %f, %f, %f, %f, %f\n", water, fast_soil_C, structural_soil_C, slow_soil_C, mineralized_soil_N, fast_soil_N ,structural_soil_L<0);
   if (!((water>=0)&(fast_soil_C>=0)&(structural_soil_C>=0)&(slow_soil_C>=0)&(mineralized_soil_N>=0)&(fast_soil_N>=0)&(passive_soil_C>=0)&(structural_soil_L>=0))) return 1;
   cohort* cc = shortest;
   while (cc!=NULL){
      //printf("%f, %f, %f, %f\n", cc->nindivs, cc->dbh, cc->balive, cc->bdead);
      if (!((cc->nindivs>=0)&(cc->dbh>=0)&(cc->balive>=0)&(cc->bdead>=0))) return 2;
      cc = cc->taller;
   }
   return 0;
}
