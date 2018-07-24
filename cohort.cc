#include <cmath>
#include <cstdio>

#include "edmodels.h"
#include "site.h"
#include "read_site_data.h"
#include "patch.h"
#include "disturbance.h"
#include "phenology.h"

#include "cohort.h"

////////////////////////////////////////////////////////////////////////////////
//! get_cohort_vm0
//! Downregulate Vm0 based on (1) light profile, and (2) day length
//!
//! @param  data Userdata structure
//! @return Nothing, updates cohort Vm0 (carboxylation rate internally)
////////////////////////////////////////////////////////////////////////////////
void cohort::get_cohort_vm0(UserData *data) {
   double cum_lai = 0.0;
   double Kn      = 0.0;
   
   Vm0 = data->Vm0_max[species];

   // If number of Vm0 bins is 1, Vm0 stays at maximum value for PFT
   if(data->num_Vm0 > 1 and data->do_downreg) {    
      // For expression relating Kn to Vm0, see Lloyd et al. 2010 (Figure 10)
      // http://www.biogeosciences.net/7/1833/2010/bg-7-1833-2010.pdf
      // Essentially, higher the Vm0 value, greater the extinction coefficient
      cohort *cc = this; 
      Kn = exp(0.00963*Vm0 - 2.43);
      if(Kn>0.0) {
         Kn *= -1.0;
      }

      // Half of the current cohort's LAI contributes to shading as well, this is
      // based on assumption that the chloroplasts on the leaves are placed mid-way
      cum_lai += lai/2.0;
      while (cc->taller != NULL) {
         cum_lai += cc->taller->lai;
         cc = cc->taller;
      }

      // Scale Vm0 of cohort based on cumulative LAI of cohorts above it
      Vm0 *= exp(Kn*cum_lai);

      // Scale Vm0 based on day length, i.e. downregulation in fall, higher latitudes
      Vm0 *= pow(siteptr->dyl_factor[data->time_period],2.0);
   }   
}

////////////////////////////////////////////////////////////////////////////////
//! get_cohort_vm0_bin
//! Get the mechanism bin corresponding to the cohort's Vm0
//!
//! @param  Vm0 Carboxylation rate for current cohort
//! @param  data Userdata structure
//! @return Integer specifying the mechanism bin
////////////////////////////////////////////////////////////////////////////////
int cohort::get_cohort_vm0_bin(double Vm0, UserData* data) {
   unsigned int index = 0;
   double high, low;
    

   // Only execute if number of Vm0 bins is > 1
   if (data->num_Vm0 > 1) {
      // Find the index corresponding to the bin above 
      for(index = 0; index < data->num_Vm0; index++) {
        high = data->Vm0_bins.at(index);
        if(high > Vm0) {
           break;
        }
      }
           
      if(index > 0) {
         low =  data->Vm0_bins.at(index-1);
         // Find distance of Vm0 from lower and higher bin
         if( abs(Vm0 - low) <= abs(high - Vm0)) {
            // If closer to lower bin then adjust index
            index--;
         }
      }
   }

   // Result is the closest bin to current Vm0
   return index;
}

////////////////////////////////////////////////////////////////////////////////
//! Cohort dynamics
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void cohort_dynamics(unsigned int t, double t1, double t2, 
                     patch** patchptr, FILE* outfile, UserData* data) {

   patch* currentp = *patchptr;
   site* currents = currentp->siteptr;
#if FTS 
   /* new_phenology is an incomplete project.
   new_phenology(t, &currentp, data);
   phenology(t, &currentp, data);
   for(currentp=*patchptr; currentp!=NULL; currentp=currentp->older)
      light_levels(&currentp, data);*/
   if(t%(COHORT_FREQ) == 0){    
      /*************/
      /* leaf fall */
      /*************/
      phenology(t, &currentp, data);  
      while (currentp != NULL) {
         light_levels(&currentp, data);
         currentp = currentp->older;
      } 
   } 
#else
   if(t%(COHORT_FREQ) == 0){    
      /*************/
      /* leaf fall */
      /*************/
      phenology(t, &currentp, data);  /* compute plant phenology */ 
      while (currentp != NULL) {
         light_levels(&currentp, data);
         currentp = currentp->older;
      } 
   } /* end t % COHORT_FREQ */
#endif


   if(!data->patch_dynamics) { /* if no pd update dist rates as will be applied in mort */
      currentp = *patchptr;
      while (currentp != NULL) {  
         calculate_disturbance_rates(t, &currentp,data);  
         currentp = currentp->older;
      } /* end loop over patches */ 
   }

   /*************************************************************/
   /*  call stiff ode integrator  */
   /*************************************************************/
   currentp = *patchptr;
   while (currentp != NULL) {    
      if(data->cd_file) {
         fprintf(outfile, "cd: integrating site %s %p patch %p \n",
                 currents->sdata->name_, currentp->siteptr, currentp); 
      }

      if ( (currentp->tallest != NULL) && (currentp->shortest != NULL) ) {
          ///CarbonConserve
          cohort* mlcc = NULL;
          double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0, all_repro_before = 0.0;
          double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0, all_repro_after = 0.0;
          double actual_dt_tc = 0.0, esti_dt_tc = 0.0, all_npp_avg_after = 0.0, all_rh_avg_after = 0.0;
          mlcc = currentp->shortest;
          while (mlcc!=NULL) {
              all_tb_before += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
              all_repro_before += (mlcc->p[0]+mlcc->p[1])*mlcc->nindivs/currentp->area;
              mlcc = mlcc->taller;
          }
          all_sc_before = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
          all_tc_before = all_tb_before + all_sc_before;
          
         currentp->okint = cm_sodeint(&currentp,t, t1, t2, data); /* call ode intgrtr */
          
          ///CarbonConserve
          mlcc = currentp->shortest;
          while (mlcc!=NULL) {
              all_tb_after += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
              all_repro_after += (mlcc->p[0]+mlcc->p[1])*mlcc->nindivs/currentp->area;
              mlcc = mlcc->taller;
          }
          all_sc_after = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
          all_tc_after = all_tb_after + all_sc_after;
          
          actual_dt_tc = all_tc_after-all_tc_before;
          esti_dt_tc = (currentp->npp_avg-currentp->rh_avg-(1.0-data->sd_mort)*all_repro_after)*data->deltat;
          
          if (abs(actual_dt_tc - esti_dt_tc)>1e-9)
          {
              printf("Carbon leakage in cm_sodeint : imbalance    %.15f actual_dt_tc %.15f esti_dt_tc     %.15f\n",actual_dt_tc-esti_dt_tc,actual_dt_tc,esti_dt_tc);
              printf("                             : patch_tc_bf  %.15f patch_sc_bf  %.15f patch_tb_bf    %.15f\n",all_tc_before,all_sc_before,all_tb_before);
              printf("                             : patch_tc_af  %.15f patch_sc_af  %.15f patch_tb_af    %.15f\n",all_tc_after,all_sc_after,all_tb_after);
              printf("                             : patch_npp_af %.15f patch_rh_af  %.15f patch_repro_af %.15f\n",currentp->npp_avg,currentp->rh_avg,all_repro_after);
              printf(" --------------------------------------------------------------------------------------\n");
          }
         if (currentp->okint) {
            printf("Fatal error integrating site %s. Site skipped okint No is %d\n", currents->sdata->name_,currentp->okint);
#if 0
             exit(0);
#endif
            currents->skip_site = 1;
            return;
         }
         if(data->cd_file) {
            if (currentp->okint != 0) 
               fprintf(outfile, "cd: exited from sode patch= %p ok=%d\n", currentp, currentp->okint);
            fprintf(outfile, "cd: exited from sode patch= %p ok=%d\n", currentp, currentp->okint);
         }
      } else {
         currentp->okint = -99;
          ///CarbonConserve
          currentp->rh = 0.0;
          currentp->rh_avg = 0.0;
          currentp->gpp = 0.0;
          currentp->gpp_avg = 0.0;
          currentp->npp = 0.0;
          currentp->npp_avg = 0.0;
         if(data->cd_file)
            fprintf(outfile, "cd: skipped integrating empty patch= %p ok=%d\n", currentp, currentp->okint); 
      }     
      currentp = currentp->older;
   } /* end loop over patches */    
  
   if((t)%(COHORT_FREQ) == 0){

      /****************************************/
      /* Stochastic Mortality and Disturbance */
      /****************************************/
       if(data->cohort_termination) {
         /***  drop cohorts  ***/      
         currentp = *patchptr;
         while (currentp != NULL) {
             ///CarbonConserve
             cohort* mlcc = NULL;
             double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0, all_repro_before = 0.0;
             double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0, all_repro_after = 0.0;
             double actual_dt_tc = 0.0, esti_dt_tc = 0.0;
             mlcc = currentp->shortest;
             while (mlcc!=NULL) {
                 all_tb_before += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
                 all_repro_before += (mlcc->p[0]+mlcc->p[1])*(1-data->sd_mort)*mlcc->nindivs/currentp->area;
                 //printf("ck chort before spp %d ba %.15f bd %.15f nin %.15f\n",mlcc->species,mlcc->balive,mlcc->bdead,mlcc->nindivs);
                 mlcc = mlcc->taller;
             }
             all_sc_before = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
             all_tc_before = all_tb_before + all_sc_before + all_repro_before;
             terminate_cohorts(&currentp->tallest,&currentp->shortest,data);
            /* following termination check for empty linked list */
             
             ///CarbonConserve
             mlcc = currentp->shortest;
             while (mlcc!=NULL) {
                 all_tb_after += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
                 all_repro_after += (mlcc->p[0]+mlcc->p[1])*(1-data->sd_mort)*mlcc->nindivs/currentp->area;
                 //printf("ck chort after spp %d ba %.15f bd %.15f nin %.15f\n",mlcc->species,mlcc->balive,mlcc->bdead,mlcc->nindivs);
                 mlcc = mlcc->taller;
             }
             all_sc_after = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
             all_tc_after = all_tb_after + all_sc_after + all_repro_after;
             
             if (abs(all_tc_after-all_tc_before)>1e-9)
             {
                 printf("Carbon leakage in terminate_cohorts : imbalance   %.15f patch_tc_af %.15f patch_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
                 printf("                                    : patch_sc_bf %.15f patch_tb_bf %.15f patch_repro_bf %.15f\n",all_sc_before,all_tb_before,all_repro_before);
                 printf("                                    : patch_sc_af %.15f patch_tb_af %.15f patch_repro_af %.15f\n",all_sc_after,all_tb_after,all_repro_after);
                 printf(" --------------------------------------------------------------------------------------\n");
             }
      
            if(data->cd_file)
               if((currentp->tallest == NULL)&&(currentp->shortest == NULL)) 
                  fprintf(outfile,"**** warning: no cohorts left in patch %p age %f area %f !\n",currentp,currentp->age,currentp->area); 
            currentp = currentp->older;
         }
      }
    
      /* sort remaining cohorts */      
      if(data->cd_file)
         fprintf(outfile,"sort...  \n");

      currentp=*patchptr;
      while (currentp != NULL){
          ///CarbonConserve
          cohort* mlcc = NULL;
          double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0;
          double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0, all_repro = 0.0;
          double actual_dt_tc = 0.0, esti_dt_tc = 0.0;
          mlcc = currentp->shortest;
          while (mlcc!=NULL) {
              all_tb_before += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
              all_repro += (mlcc->p[0]+mlcc->p[1])*mlcc->nindivs/currentp->area;
              mlcc = mlcc->taller;
          }
          all_sc_before = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
          all_tc_before = all_tb_before + all_sc_before;
         sort_cohorts(&currentp, data);
          
          ///CarbonConserve
          mlcc = currentp->shortest;
          while (mlcc!=NULL) {
              all_tb_after += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
              mlcc = mlcc->taller;
          }
          all_sc_after = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
          all_tc_after = all_tb_after + all_sc_after;

          actual_dt_tc = all_tc_after - all_tc_before;
          esti_dt_tc = (currentp->npp_avg-currentp->rh_avg-(1.0-data->sd_mort)*all_repro)/12.0;

          if (abs(all_tc_after-all_tc_before)>1e-9)
          {
              printf("Carbon leakage in sort_cohorts : imbalance   %.15f patch_tc_af %.15f patch_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
              printf("                               : patch_sc_bf %.15f patch_tb_bf %.15f\n",all_sc_before,all_tb_before);
              printf("                               : patch_sc_af %.15f patch_tb_af %.15f\n",all_sc_after,all_tb_after);
              printf(" --------------------------------------------------------------------------------------\n");
          }
         currentp = currentp->older;
      }    

   }/*end t%cohort_freq*/
 
   /****************/
   /* Reproduction */
   /****************/
   if((t)%(COHORT_FREQ) == 0){ 
   if(data->cd_file)
      fprintf(outfile,"repro...  \n");
    
      currentp=*patchptr;
       ///CarbonConserve
       patch* mlcp = NULL;
       double all_tb_before = 0.0, all_sc_before = 0.0, all_tc_before = 0.0;
       double all_tb_after = 0.0, all_sc_after = 0.0, all_tc_after = 0.0;
       mlcp = *patchptr;
       while (mlcp!=NULL) {
           cohort* mlcc = mlcp->shortest;
           double tmp_tb = 0.0;
           while (mlcc!=NULL) {
               tmp_tb += (mlcc->balive+mlcc->bdead+(mlcc->p[0]+mlcc->p[1])*(1.0-data->sd_mort)*data->deltat)*mlcc->nindivs/mlcp->area;
               mlcc=mlcc->taller;
           }
           all_tb_before +=tmp_tb*mlcp->area/data->area;
           all_sc_before += (mlcp->fast_soil_C+mlcp->slow_soil_C+mlcp->structural_soil_C+mlcp->passive_soil_C)*mlcp->area/data->area;
           mlcp = mlcp->older;
       }
       all_tc_before = all_tb_before + all_sc_before;
      reproduction(t,&currentp,data);
       ///CarbonConserve
       mlcp = *patchptr;
       while (mlcp!=NULL) {
           cohort* mlcc = mlcp->shortest;
           double tmp_tb = 0.0;
           while (mlcc!=NULL) {
               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/mlcp->area;
               mlcc=mlcc->taller;
           }
           all_tb_after += tmp_tb*mlcp->area/data->area;
           for (int spp=0;spp<NSPECIES;spp++)
           {
                all_tb_after +=mlcp->repro[spp]*mlcp->area/data->area;;
           }
           all_sc_after += (mlcp->fast_soil_C+mlcp->slow_soil_C+mlcp->structural_soil_C+mlcp->passive_soil_C)*mlcp->area/data->area;
           mlcp = mlcp->older;
       }
       all_tc_after = all_tb_after+all_sc_after;

       if (abs(all_tc_after-all_tc_before)>1e-9)
       {
           printf("Carbon leakage in reproduction : imbalance  %.15f site_tc_af %.15f site_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
           printf("                               : site_sc_bf %.15f site_tb_bf %.15f\n",all_sc_before,all_tb_before);
           printf("                               : site_sc_af %.15f site_tb_af %.15f\n",all_sc_after,all_tb_after);
           printf(" --------------------------------------------------------------------------------------\n");
       }
       
   }
  

   if((t)%(COHORT_FREQ) == 0){ 
      /* spawn new cohorts */
   if(data->cd_file)
      fprintf(outfile,"spawn...\n");

      currentp = *patchptr;
       patch* mlcp = NULL;
       double all_tb_before = 0.0, all_sc_before = 0.0, all_tc_before = 0.0;
       double all_tb_after = 0.0, all_sc_after = 0.0, all_tc_after = 0.0;
       
       mlcp = *patchptr;
       while (mlcp!=NULL) {
           cohort* mlcc = mlcp->shortest;
           double tmp_tb = 0.0;
           while (mlcc!=NULL) {
               tmp_tb += (mlcc->balive+mlcc->bdead+(mlcc->p[0]+mlcc->p[1])*(1-data->sd_mort)*data->deltat)*mlcc->nindivs/mlcp->area;
               mlcc=mlcc->taller;
           }
           all_tb_before +=tmp_tb*mlcp->area/data->area;
           all_sc_before +=(mlcp->fast_soil_C+mlcp->slow_soil_C+mlcp->structural_soil_C+mlcp->passive_soil_C)*mlcp->area/data->area;
           mlcp = mlcp->older;
       }
       all_tc_before = all_tb_before + all_sc_before;
      while (currentp != NULL){
          spawn_cohorts(t,&currentp,data);
         currentp = currentp->older;
      }
       
       ///CarbonConserve
       mlcp = *patchptr;
       while (mlcp!=NULL) {
           cohort* mlcc = mlcp->shortest;
           double tmp_tb = 0.0;
           while (mlcc!=NULL) {
               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/mlcp->area;
               mlcc=mlcc->taller;
           }
           all_tb_after += tmp_tb*mlcp->area/data->area;
           all_sc_after += (mlcp->fast_soil_C+mlcp->slow_soil_C+mlcp->structural_soil_C+mlcp->passive_soil_C)*mlcp->area/data->area;
           mlcp = mlcp->older;
       }
       all_tc_after = all_tb_after + all_sc_after;

       if (abs(all_tc_after-all_tc_before)>1e-9)
       {
           printf("Carbon leakage in spawn_cohorts : imbalance  %.15f site_tc_af %.15f site_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
           printf("                                : site_sc_bf %.15f site_tb_bf %.15f\n",all_sc_before,all_tb_before);
           printf("                                : site_sc_af %.15f site_tb_af %.15f\n",all_sc_after,all_tb_after);
           printf(" --------------------------------------------------------------------------------------\n");
       }
    
      if(data->cohort_fusion) {
         /*fuse cohorts*/
         if(data->cd_file)
            fprintf(outfile,"fusing cohorts... \n");

         currentp=*patchptr;
         fuse_cohorts(&currentp, data);
      }
    
      if(data->cohort_fission) {
         /*split cohorts*/
         if(data->cd_file)
            fprintf(outfile,"splitting cohorts... \n");   

         currentp=*patchptr;
         split_cohorts(&currentp,data);
      }
   } /*end t%COHORT_FREQ*/
  
   return;
  
}

////////////////////////////////////////////////////////////////////////////////
//! init_cohorts
//! initialize new cohort
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void init_cohorts(patch** patchptr, UserData* data){
   /*this function is called only during model initialization from init_patches.*/

   patch* cp = *patchptr;
   site* cs = cp->siteptr;

   /* add new cohorts to linked list of cohorts */
   for (int spp=0; spp<NSPECIES; spp++) {
      /*initial number of seedlings of each spp in patch */    
      double nindivs = data->initial_density[spp] * (*patchptr)->area; 

      double lat = cs->sdata->lat_;
      /* rule out tropical spp out of tropics */
      if ( ((lat > data->tropic_n_limit) || (lat < data->tropic_s_limit))   && (data->is_tropical[spp])) {
         nindivs = 0.000001;
      }
      /* rule out temp spp in tropics */
      else if ( ((lat <= data->tropic_n_limit) && (lat >= data->tropic_s_limit)) && (!data->is_tropical[spp] and !data->is_grass[spp]) ) {
         nindivs = 0.000001;
      }
    
#if LANDUSE
      /* plant gasses only cropland, and at much higher density */
      if (cp->landuse == LU_CROP) { 
         if (spp < 2) {
            nindivs *= 100.0;
         }
         if (spp > 1) {
            nindivs = 0.000001;
         }
      }
#endif

      /* creat a dummy cohort to figure out size of biomass compartments */  
      cohort* dc = (cohort*) malloc(sizeof(cohort));
      if (dc == NULL) {
         printf("out of memory\n");
         while(getchar() != '\n');
      }

      dc->species = spp; 
      dc->nindivs = nindivs;
      if(data->allometry_type == 0 or data->is_tropical[spp] or data->is_grass[spp]) {
         dc->hite = (data->hgt_min[spp]); 
      } else {
         dc->hite = (data->min_hgt[spp]);
      }             
      dc->dbh = dc->Dbh(data);    
      /* calculate size of biomass compartments */
      dc->bdead = dc->Bdead(data);
      dc->balive= dc->Bleaf(data)*(1.0 + data->q[spp] + data->qsw[dc->species]*dc->hite);
      dc->b = dc->balive + dc->bdead;

      /* create appropritaely sized cohort then free dummy cohort */
      create_cohort(dc->species,dc->nindivs,dc->hite,dc->dbh,dc->balive,dc->bdead,patchptr,data);  
      free(dc);
   }  /* end loop over species */
}

////////////////////////////////////////////////////////////////////////////////
//! create_cohort
//! create new cohort
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void create_cohort (unsigned int spp, double nindivs, double hite, double dbh, 
                    double balive, double bdead, patch** patchptr, 
                    UserData* data) {   
   /* this function is called in model initialization from init_cohorts and *
    * read cohort distribution and during the model run from spawn_cohorts  */
  
   cohort* newcohort = (cohort*) malloc(sizeof(cohort));
   if (newcohort == NULL) {
      printf("out of memory\n");
      while( getchar() != '\n');
   }
   /* assign cohort attributes */
   newcohort->siteptr  = (*patchptr)->siteptr;
   newcohort->patchptr = *patchptr;
   newcohort->species  = spp;
   newcohort->pt       = data->pt[spp];
   newcohort->gpp      = 0.0;
   newcohort->npp      = 0.0;
    //checkstep
    newcohort->gpp_avg = 0.0;
    newcohort->npp_avg = 0.0;
    newcohort->md_avg = 0.0;
   newcohort->npp2     = 0.0;
   newcohort->md       = 0.0;
   newcohort->fs_open  = 1.0;
   newcohort->fsn      = 1.0;
   newcohort->fsw      = 1.0;
   for(size_t i=0; i<N_CLIMATE; i++) {
      newcohort->cbr[i]    = 1.0;
      newcohort->cb[i]     = 1.0;
      newcohort->cb_toc[i] = 1.0;
   }
   newcohort->cbr_bar= 1.0;
   newcohort->hite               = hite;
   newcohort->dbh                = dbh;
   newcohort->p[0]               = 0.0;
   newcohort->p[1]               = 0.0;
   newcohort->status             = 0;

   /* initialize derivatives to zero */
   newcohort->dndt               = 0.0;
   newcohort->ddbhdt             = 0.0;
   newcohort->dhdt               = 0.0;
   newcohort->dbalivedt          = 0.0;
   newcohort->dbdeaddt           = 0.0;
   newcohort->payment_to_Nfixers = 0.0;
  
   /* calculate size of biomass compartments */
   newcohort->bdead = bdead;
   newcohort->balive = balive;
   
   // if(data->restart) {
   //   if(data->allometry_type == 0) {
   //      newcohort->hite = (data->hgt_min[spp]); 
   //   } else {
   //      newcohort->hite = (data->min_hgt[spp]);
   //   }             
   //   newcohort->dbh = newcohort->Dbh(data);
   //   newcohort->bdead = newcohort->Bdead(data);
   //   newcohort->balive= newcohort->Bleaf(data)*(1.0 + data->q[spp] + data->qsw[spp]*newcohort->hite);
   // }

   newcohort->b = newcohort->balive + newcohort->bdead;
  
   newcohort->bl = (1.0 / (1.0 + data->q[spp] + data->qsw[newcohort->species]
                           * newcohort->hite) ) * newcohort->balive;
   newcohort->blv = 0.0;
   newcohort->br = (data->q[newcohort->species] 
                    / (1.0 + data->q[spp] + data->qsw[newcohort->species] * newcohort->hite) )
      * newcohort->balive;
   newcohort->bsw = (data->qsw[newcohort->species] * newcohort->hite
                     / (1.0 + data->q[spp] + data->qsw[newcohort->species] * newcohort->hite) ) 
      * newcohort->balive;

   newcohort->bs  = newcohort->bsw + newcohort->bdead;   
   newcohort->bstem = data->agf_bs*newcohort->bs;

   newcohort->babove = newcohort->bl + data->agf_bs * newcohort->bs; 
   newcohort->bbelow = newcohort->br + (1.0 - data->agf_bs) * newcohort->bs; 

   /***********************/
   /* Nitrogen in recruit */
   /***********************/ 
   /*101398-GCH: this calc really only needs to be done once for each 
     species*/

   data->c2n_recruit[spp] = newcohort->b 
      / (newcohort->balive / data->c2n_leaf[spp] + newcohort->bdead/data->c2n_stem);
  
   /*printf("spp %d c2n_recruit %f\n",spp,data->c2n_recruit[spp]);*/
   /****************************/
   newcohort->resp         = 0.0;
   newcohort->gr_resp      = 0.0;
   newcohort->An_pot       = 0.0;
   newcohort->An_shut      = 0.0;
   newcohort->E_pot        = 0.0;
   newcohort->E_shut       = 0.0;
   newcohort->water_uptake = 0.0;

   newcohort->nindivs = nindivs; 
  
   newcohort->lai = (newcohort->nindivs) 
      * (newcohort->bl * data->specific_leaf_area[newcohort->species])
      * (1.0 / ((*patchptr)->area));
   insert_cohort(&newcohort, &(*patchptr)->tallest, &(*patchptr)->shortest,data);
   if(data->num_Vm0  > 1) {
      // Multiple mechanism file case
      newcohort->get_cohort_vm0(data);      
      newcohort->Vm0_bin = newcohort->get_cohort_vm0_bin(newcohort->Vm0, data);
   } else {
      // Single mechanism file case
      newcohort->Vm0_bin = 0;
   }
}

////////////////////////////////////////////////////////////////////////////////
//! next_taller
//! finds next taller cohort with patch
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
cohort* next_taller(cohort* current, double* stp){

   /***  starting with shortest tree on the grid, find tree just  ***/
   /***  taller than tree being considered and return its pointer ***/
   while (current != NULL && current->hite < *stp)
      current = current->taller;
   
   return (current);
}

////////////////////////////////////////////////////////////////////////////////
//! terminate_cohorts
//! terminates cohorts
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void terminate_cohorts(cohort** ptallest, cohort** pshortest, UserData* data){
   cohort* currentc = *ptallest;
   while (currentc != NULL){
      cohort* nextc = currentc->shorter;
      /***  cohort size is below threshold   ***/
      /***  remove & adjust size relations among remaining cohorts  ***/
      if(((currentc->nindivs/currentc->patchptr->area)*currentc->b) < data->btol){
         if (currentc->taller == NULL) *ptallest = currentc->shorter;
         else (currentc->taller)->shorter = currentc->shorter;
         if (currentc->shorter == NULL) *pshortest = currentc->taller;
         else (currentc->shorter)->taller = currentc->taller;
         free (currentc);
      }
      currentc = nextc;
   }
   return;
}

////////////////////////////////////////////////////////////////////////////////
//! spawn_cohorts
//! spawn new cohorts of jueveniles
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void spawn_cohorts(unsigned int t, patch** patchptr, UserData* data){  
   patch* cp= *patchptr;  
   /* add new cohorts to linked list of cohorts */
   /* dummy cohort is used to figure out biomass */
   for(int spp=0;spp<NSPECIES;spp++){
      cohort* dc = (cohort*) malloc(sizeof(cohort));
      if (dc == NULL) {
         printf("out of memory\n");
         while(getchar() != '\n');
      }
      dc->species = spp; 
      if(data->allometry_type == 0 or data->is_tropical[spp] or data->is_grass[spp]) {
         dc->hite = (data->hgt_min[spp]); 
      } else {
         dc->hite = (data->min_hgt[spp]);
      }
      dc->dbh = dc->Dbh(data);     
      /* calculate size of biomass compartments */
      dc->bdead = dc->Bdead(data);
      dc->balive= dc->Bleaf(data)*(1.0 + data->q[spp] + data->qsw[dc->species]*dc->hite);
      dc->b = dc->balive + dc->bdead;
      if(data->internal_recruitment)
         dc->nindivs = cp->area*cp->repro[spp]/dc->b;        

      if(data->external_recruitment) {
        /* add in external recruitment */
        dc->nindivs += data->seed_rain[dc->species];
      }
      if(dc->nindivs  > 0.0) {
         create_cohort(dc->species,dc->nindivs,dc->hite,dc->dbh,dc->balive,dc->bdead,&cp,data); 
      }

      cp->repro[spp] = 0.0;   /* reset reprodictive array */    
      free(dc);        
   }  /* end loop over species */
}

////////////////////////////////////////////////////////////////////////////////
//! split_cohorts
//! split cohorts that get too big
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void split_cohorts (patch** patchptr, UserData* data) {
   // check to see if any cohorts are too big, then split them in two
   // if they are and assign new cohort a height = height + epsilon

   patch* currentp = *patchptr;
   while(currentp != NULL){
       ///CarbonConserve
       cohort* mlcc = NULL;
       double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0;
       double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0, all_repro = 0.0;
       double actual_dt_tc = 0.0, esti_dt_tc = 0.0;
       mlcc = currentp->shortest;
       while (mlcc!=NULL) {
           all_tb_before += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
           all_repro += (mlcc->p[0]+mlcc->p[1])*mlcc->nindivs/currentp->area;
           mlcc = mlcc->taller;
       }
       all_sc_before = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
       all_tc_before = all_tb_before + all_sc_before;
       
      cohort* currentc=currentp->tallest;
      while(currentc != NULL){
         if((currentc->lai > data->lai_tol)){           
            cohort* copyc = (cohort*) malloc(sizeof(cohort));
            if (copyc == NULL) {
               printf("sc: out of memory\n");
               while(getchar() != '\n');
            }

            /*copy cohort*/
            copy_cohort(&currentc,&copyc);
            // Give each cohort half of the individuals and half lai
            // Only those variables that have units per area (m2) should be 
            // included below for halving. Variables with units per plant should
            // not be included since we are halving area not plants.
            copyc->nindivs *= 0.5;
            copyc->lai *= 0.5;
            currentc->nindivs *= 0.5;
            currentc->lai *= 0.5;            
            copyc->dndt *=0.5;
            currentc->dndt *=0.5;
            currentc->dbh -= 0.001; 
            copyc->dbh += 0.001;
            
            /*insert*/
            insert_cohort(&copyc,&(currentp->tallest),&(currentp->shortest),data);  
            // Update the Vm0 value for the cohort since changing the height
            // changes the light profile
            copyc->get_cohort_vm0(data);
            currentc->get_cohort_vm0(data);
         } /* end if on lai */
         currentc=currentc->shorter;
      }  /*  end loop over cohorts */
       ///CarbonConserve
       mlcc = currentp->shortest;
       while (mlcc!=NULL) {
           all_tb_after += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
           mlcc = mlcc->taller;
       }
       all_sc_after = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
       all_tc_after = all_tb_after + all_sc_after;

       actual_dt_tc = all_tc_after-all_tc_before;
       esti_dt_tc = (currentp->npp_avg-currentp->rh_avg-(1.0-data->sd_mort)*all_repro)*data->deltat;
       
       if (abs(all_tc_after-all_tc_before)>1e-9)
       {
           printf("Carbon leakage in split_cohorts : imbalance    %.15f patch_tc_af %.15f patch_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
           printf("                                : patch_sc_bf  %.15f patch_tb_bf %.15f\n",all_sc_before,all_tb_before);
           printf("                                : patch_sc_af  %.15f patch_tb_af %.15f\n",all_sc_after,all_tb_after);
           printf("                                : actual_dt_tc %.15f esti_dt_tc  %.15f\n",actual_dt_tc,esti_dt_tc);
           printf(" --------------------------------------------------------------------------------------\n");
       }
       
      currentp=currentp->older;
   }   /* end loop over patches */
}


////////////////////////////////////////////////////////////////////////////////
//! copy_cohort
//! copy cohorts
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void copy_cohort(cohort** currentc, cohort** copyc){

   /* The easier way to do this would be to bitmap it */

   cohort* o = *currentc;
   cohort* n = *copyc;

   for (size_t i=0; i<N_CLIMATE; i++) {
      n->cbr[i]    = o->cbr[i];
      n->cb[i]     = o->cb[i];
      n->cb_toc[i] = o->cb_toc[i];
   }

   n->cbr_bar = o->cbr_bar;
   n->species = o->species;
   n->pt = o->pt;
   n->nindivs = o->nindivs;     /* number of individuals in cohort      */
   n->dndt = o->dndt;           /* rate of change of cohort size        */
   n->dbh = o->dbh;              /* dbh in cm                            */
   n->hite= o->hite;            /* height in meters                     */
   n->b = o->b;                    /* total biomass per indiv              */
   n->babove = o->babove;       /* total above ground biomass           */
   n->bbelow = o->bbelow;       /* total below ground biomass           */
   n->balive = o->balive;       /* total living biomass per indiv       */
   n->bdead = o->bdead;         /* dead biomass per indiv               */
   n->bsw = o->bsw;              /* sapwood in stem and roots            */
   n->bl = o->bl;                  /* leaf biomass per indiv               */
   n->blv = o->blv;             /* leaf biomass per indiv               */
   n->bs = o->bs;               /* structural biomass per indiv stem + structual roots */
   n->bstem = o->bstem;         /* stem biomass per indiv               */
   n->br = o->br;               /* fine root biomass per indiv          */
   n->leaf_area = o->leaf_area;  /* leaf area of plant                   */
   n->lai = o->lai;              /* leaf area index of plant             */
   n->lite = o->lite;            /* light level for cohort               */
   n->gpp = o->gpp;             /* net primary PER PLANT!(kg/plant/yr)  */
   n->npp = o->npp;             /* net primary PER PLANT!(kg/plant/yr)  */
   n->npp2 = o->npp2;              /* net primary PER PLANT!(kg/plant/yr)  */
    //checkstep
    n->gpp_avg = o->gpp_avg;
    n->npp_avg = o->npp_avg;
    n->md_avg = o->md_avg;
   n->resp = o->resp;              /* plant respiration                    */
   n->gr_resp = o->gr_resp;
   n->md = o->md;
   n->An_pot = o->An_pot;
   n->An_shut = o->An_shut;
   n->fs_open = o->fs_open; /* fraction of month with stomates open (dimensionless) */
   n->fsn = o->fsn;         /* degree of n limitation (dimensionless) */
   n->fsw = o->fsw;         /* degree of water limitation (dimensionless) */



   n->carbon_balance = o->carbon_balance; /* relative carbon balance: ratio of bl to bl* */
   n->p[0]           = o->p[0];           /* seed reproduction */
   n->p[1]           = o->p[1];           /* clonal reproduction */
   n->status         = o->status;         /* growth status of plant */

   /* water fields */
   n->water_uptake = o->water_uptake;
   n->E_pot        = o->E_pot;
   n->E_shut       = o->E_shut;

   /* nitrogen fields */
   n->nitrogen_uptake = o->nitrogen_uptake;
   n->N_uptake_pot    = o->N_uptake_pot;
   n->N_uptake_shut   = o->N_uptake_shut;

   /* variables needed for integration */
   n->dndt      = o->dndt;      /* time derivative of cohort size             */
   n->dhdt      = o->dhdt;      /* time derivative of height                  */
   n->ddbhdt    = o->ddbhdt;    /* time derivative of dbh                     */
   n->dbalivedt = o->dbalivedt; /* time derivative of total living biomass    */
   n->dbdeaddt  = o->dbdeaddt;  /* time derivative of dead biomass            */
  
  
   /*linked list fields*/
   n->taller     = NULL;        /* pointer to next tallest cohort             */
   n->shorter    = NULL;        /* pointer to next shorter cohort             */
   n->patchptr   = o->patchptr; /* pointer to patch that cohort is in         */
   n->siteptr    = o->siteptr;  /* pointer to site that cohort is in          */
}

////////////////////////////////////////////////////////////////////////////////
//! fuse_cohorts
//! join similar cohorts
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void fuse_cohorts (patch** patchptr, UserData* data) {
   /* fuse cohorts of same spp and similar heights     */
   /* then resort lists by height if a fusion occurred */

   unsigned int fusion_took_place = 0;

   patch* cp = *patchptr;
   /* loop over all cohorts */
   while (cp != NULL) {
       ///CarbonConserve
       cohort* mlcc = NULL;
       double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0;
       double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0, all_repro = 0.0;
       double actual_dt_tc = 0.0, esti_dt_tc = 0.0;
       mlcc = cp->shortest;
       while (mlcc!=NULL) {
           all_tb_before += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/cp->area;
           all_repro += (mlcc->p[0]+mlcc->p[1])*mlcc->nindivs/cp->area;
           mlcc = mlcc->taller;
       }
       all_sc_before = cp->fast_soil_C+cp->slow_soil_C+cp->structural_soil_C+cp->passive_soil_C;
       all_tc_before = all_tb_before + all_sc_before;
       
      cohort* cc = cp->tallest;
      while (cc != cp->shortest) {
         cohort* nextc = cc->shorter;
         cohort* nextnextc = nextc->shorter;
         /* calculate total density within adjacent cohorts */ 
         /* totaln= (cc->nindivs + nextc->nindivs)/cp->area; */
         double total_lai;
         // If status == 5 i.e. leaves have fallen, compute lai using blv instead of bl
         if(cc->status == 0) {
            total_lai = cc->lai + nextc->lai;
         } else {
            total_lai = (cc->nindivs*cc->blv + nextc->nindivs*nextc->blv)*
                    data->specific_leaf_area[cc->species]/cp->area;
         }
         
         while ( (ABS(cc->dbh - nextc->dbh) / (0.5 * (cc->dbh + nextc->dbh)) < data->fusetol) 
                 && (nextc != cc)) {
            
            if ((cc->species == nextc->species) && (total_lai < data->lai_tol)) {
#if 0
               printf("fusion: FUSION TAKING PLACE %s\n", cp->siteptr->name);
               printf("fusion: fusing c %p with nextc %p\n",cc,nextc);
               printf("fusion: c %p: spp= %d h= %f n=%f\n",cc,cc->species,cc->hite,cc->nindivs);
               printf("fusion: nextc %p: spp= %d h= %f n=%f\n",
                      nextc,nextc->species,nextc->hite,nextc->nindivs);
#endif
               fusion_took_place = 1;
                
               /* update current: add individuals */
               double newn = cc->nindivs + nextc->nindivs;
                
               /* update current: weighted average of biomasses */
               double newbalive = (cc->nindivs * cc->balive + nextc->nindivs * nextc->balive) / newn;
               double newbdead  = (cc->nindivs * cc->bdead + nextc->nindivs * nextc->bdead) / newn;
               double newblv    = (cc->nindivs * cc->blv + nextc->nindivs * nextc->blv) / newn;
               double newb      = (cc->nindivs * cc->b + nextc->nindivs * nextc->b) / newn;
               double newbl     = (cc->nindivs * cc->bl + nextc->nindivs * nextc->bl) / newn;
               double newbs     = (cc->nindivs * cc->bs + nextc->nindivs * nextc->bs) / newn;
               double newbr     = (cc->nindivs * cc->br + nextc->nindivs * nextc->br) / newn;
               double newbstem  = (cc->nindivs * cc->bstem + nextc->nindivs * nextc->bstem) / newn;
               double newbsw    = (cc->nindivs * cc->bsw + nextc->nindivs * nextc->bsw) / newn;
               double newh      = (cc->nindivs * cc->hite + nextc->nindivs * nextc->hite) / newn;
               double newdbh    = (cc->nindivs * cc->dbh + nextc->nindivs * nextc->dbh) / newn;
               
               /* update current: weighted average of derivatives  */
               cc->dndt      = (cc->nindivs*cc->dndt + nextc->nindivs*nextc->dndt)/newn;
               cc->dhdt      = (cc->nindivs*cc->dhdt + nextc->nindivs*nextc->dhdt)/newn;          
               cc->ddbhdt    = (cc->nindivs*cc->ddbhdt + nextc->nindivs*nextc->ddbhdt)/newn;    
               cc->dbalivedt = (cc->nindivs*cc->dbalivedt + nextc->nindivs*nextc->dbalivedt)/newn;
               cc->dbdeaddt  = (cc->nindivs*cc->dbdeaddt + nextc->nindivs*nextc->dbdeaddt)/newn;

               /* set cbr to average of the two cohorts */
               cc->cbr_bar = (cc->cbr_bar*cc->nindivs + nextc->nindivs*nextc->cbr_bar)/newn;
               for(int i=0;i<N_CLIMATE;i++) { 
                  cc->cb[i] = (cc->cb[i]*cc->nindivs + nextc->nindivs*nextc->cb[i])/newn;
               }
               for(int i=0;i<N_CLIMATE;i++) {
                  cc->cb_toc[i] =(cc->cb_toc[i]*cc->nindivs + nextc->nindivs*nextc->cb_toc[i])/newn;
               }
               
               /* update current: implied characteristics */
               cc->balive  = newbalive;
               cc->bdead   = newbdead;
               cc->blv     = newblv;
               cc->bstem   = newbstem;
               cc->bsw     = newbsw;
               cc->bs      = newbs;
               cc->bl      = newbl;
               cc->br      = newbr;
               cc->b       = newb;
               cc->babove  = cc->bl + data->agf_bs*cc->bs;
               cc->nindivs = newn;    
               cc->hite    = newh;    
               cc->dbh     = newdbh;
               
               /*delete nextc*/  
               (nextc->taller)->shorter = nextnextc;  
               if (nextc->shorter == NULL){
                  cp->shortest = nextc->taller;
               }
               else{
                  nextnextc->taller = nextc->taller;
               }
        
               /* free memory of nextc */
               free (nextc);
            } /* end if */

            if(nextnextc != NULL){
               nextc = nextnextc;
               nextnextc = nextc->shorter;
            } else {
               nextc = cc;
            }
         } /* end checking other cohorts loop */
         if (cc != cp->shortest) {
            cc = cc->shorter;
         }

      } /* end cohort loop */
       
       ///CarbonConserve
       mlcc = cp->shortest;
       while (mlcc!=NULL) {
           all_tb_after += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/cp->area;
           mlcc = mlcc->taller;
       }
       all_sc_after = cp->fast_soil_C+cp->slow_soil_C+cp->structural_soil_C+cp->passive_soil_C;
       all_tc_after = all_tb_after + all_sc_after;
       
       actual_dt_tc = all_tc_after-all_tc_before;
       esti_dt_tc = (cp->npp_avg-cp->rh_avg-(1.0-data->sd_mort)*all_repro)/12.0;
       
       if (abs(all_tc_after-all_tc_before)>1e-9)
       {
           printf("Carbon leakage in fuse_cohorts  : imbalance    %.15f patch_tc_af %.15f patch_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
           printf("                                : patch_sc_bf  %.15f patch_tb_bf %.15f\n",all_sc_before,all_tb_before);
           printf("                                : patch_sc_af  %.15f patch_tb_af %.15f\n",all_sc_after,all_tb_after);
           printf("                                : actual_dt_tc %.15f esti_dt_tc  %.15f\n",actual_dt_tc,esti_dt_tc);
           printf(" --------------------------------------------------------------------------------------\n");
       }

      cp = cp->older;
   } /* ends patch loop */

   if (fusion_took_place == 1) {  /* if fusion(s) occured sort cohorts */
      cp = *patchptr;
      while (cp != NULL) {  /* loop over patches */
         sort_cohorts(&cp, data);
         cp = cp->older;
      } /* end loop over patches */
   } /* end if */
    
}

////////////////////////////////////////////////////////////////////////////////
//! sort_cohorts
//! sort cohorts
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void sort_cohorts (patch** patchptr, UserData* data) {  
   /* resort cohort lists within each patch by straight insertion */
  
   patch* current_patch= *patchptr;
   /* loop over cohorts */
   cohort* tallestc  = NULL;
   cohort* shortestc = NULL;
  
   cohort* current_c = current_patch->tallest; 
  
   while (current_c != NULL){  /* loop over current list */     
      cohort* next_c = current_c->shorter;
      insert_cohort(&current_c,&tallestc,&shortestc,data);
      current_c=next_c;
   }
   current_patch->tallest = tallestc;
   current_patch->shortest = shortestc;
}

////////////////////////////////////////////////////////////////////////////////
//! insert_cohort
//! insert cohort into linked list
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void insert_cohort (cohort** pcurrentc, cohort** ptallest, 
                    cohort** pshortest, UserData* data) {

   cohort* icohort = *pcurrentc; /* assign address to icohort local name */

   /* place in the correct place in the linked list of heights */
   /* begin by finding cohort that is just taller than the new cohort */
   double tsp = icohort->hite;
   cohort* tallptr = next_taller(*pshortest, &tsp);
   cohort* shortptr = NULL;

   /* new cohort is tallest */
   if(tallptr==NULL){
      /* new shorter cohort to the new cohort is the old tallest cohort */ 
      shortptr = *ptallest;
      /* new cohort is tallest cohort and next taller remains null */
      *ptallest = icohort;
   }
   /* new cohort is not tallest */
   else{
      /* next shorter cohort to new cohort is the next shorter cohort */
      /* to the cohort just taller than the new cohort */
      shortptr = tallptr->shorter;
    
      /* new cohort becomes the next shorter cohort to the cohort */
      /* just taller than the new cohort */
      tallptr->shorter = icohort;
   }

   /* new cohort is shortest */ 
   if(shortptr==NULL){
      /*  next shorter reamins null */
      /* cohort is placed at the bottom of the list */
      *pshortest = icohort; 
   } 
   else{
      /* new cohort is not shortest and becomes next taller cohort */
      /* to the cohort just below it as defined in the previous block */
      shortptr->taller = icohort;
   }

   /* assign taller and shorter links for the new cohort */ 
   icohort->taller = tallptr;
   icohort->shorter = shortptr;

#if 0
   printf("** ic: address of icohort taller = %p \n",icohort->taller);
   printf("** ic: address of icohort shorter = %p \n",icohort->shorter);
   printf("** ic: address of shorter = %p \n",shortptr);
   printf("** ic: address of taller = %p \n",tallptr);
   printf("** ic: address shortest cohort  = %p \n",*pshortest);
   printf("** ic: address tallest cohort = %p \n",*ptallest);
#endif 
}

////////////////////////////////////////////////////////////////////////////////
//! reproduction
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void reproduction (unsigned int t, patch** patchptr, UserData* data) {

   patch* cp = *patchptr;
#if LANDUSE
   site* cs = cp->siteptr;
#endif

   /* for each target patch */
   patch* target = cp;
   while (target != NULL) {
      /* goto every other patch */
      patch* source = cp;
      while (source !=NULL) {
         /* add in repro from all cohorts weighted by the fraction of total area that is the focal patch */
         cohort* cc = source->tallest;
         while (cc != NULL) {  
            /* global dispersal */
            size_t spp = cc->species;  

#if LANDUSE
            /* this says dispersal stays within landuse type      */
            /* first term is abs amount produced                  */
            /* second term is fraction that lands on target site  */
            /* last term converts xamount to amount per unit area */
            target->repro[spp] += ((data->m[spp] * cc->nindivs
                                    * (1.0 - data->sd_mort) * cc->p[0])
                                   * (target->area / (cs->area_fraction[target->landuse] * data->area))
                                   / (target->area)) * (data->deltat*COHORT_FREQ);
#else
            /* this assumes only one landuse type */
            target->repro[spp] += ((data->m[spp] * cc->nindivs 
                                    * (1.0 - data->sd_mort) * cc->p[0])
                                   * (target->area / (data->area)) / (target->area)) * (data->deltat*COHORT_FREQ);
#endif
            if(source==target) /* if source = target add in local dispersal */
               target->repro[spp] += (cc->nindivs * ((1.0 - data->m[spp])
                                                     * (1.0 - data->sd_mort)
                                                     * cc->p[0] + cc->p[1])
                                      / (target->area)) * (data->deltat*COHORT_FREQ);

            cc = cc->shorter;  
         }
         source = source->older;
      }
      target = target->older;
   }  /* end loop over targets */
}


/******************************************************************************/
