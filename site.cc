 #include <cmath>
#include <cstring>

#include "edmodels.h"
#include "site.h"
#include "patch.h"
#ifdef ED
#include "cohort.h"
#endif
#include "disturbance.h"
#include "restart.h"
#include "read_site_data.h"
#if LANDUSE
#include "landuse.h"
#endif
#ifdef COUPLED // TODO: this makes things messy
#include "../iGLM/glm_coupler.h"
#endif


int model_site(size_t lat, size_t lon, size_t counter, UserData* data);
void update_site_landuse(site** siteptr, size_t lu, UserData* data);
#ifdef ED
void species_site_size_profile(site** pcurrents, unsigned int nbins, UserData* data);
void update_height_profiles(site** current_site, UserData* data);
void modeled_height_distribution(site** current_site, UserData* data); 
void observed_height_distribution(site** current_site, UserData* data);
#endif


using namespace std;

////////////////////////////////////////////////////////////////////////////////
//! radians
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double radians(double degree) {
   return degree * M_PI/180.0;
}

////////////////////////////////////////////////////////////////////////////////
//! degrees
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double degrees(double radians) {
   return radians * 180.0/M_PI;   
}

////////////////////////////////////////////////////////////////////////////////
//! get_day_length
//! Formula for calculating day length is from
//! ocean.stanford.edu/courses/EESS151/activity2/daylength.xls
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double get_day_length(double lat, double month, int get_max) {
   double day_length = 0.0;
   double data_angle_rad = 0.0;
   double decl_deg = 0.0;
   double decl_rad = 0.0;
   double lat_rad = 0.0;
   double cond = 0.0;
   
   // Convert date into angle-rad
   data_angle_rad = radians(360.0*month/365.0);
   // Compute declination from day of month
   if (!get_max) {
   decl_deg = 0.39637 - 22.9133*cos(data_angle_rad) + 4.02543*sin(data_angle_rad)
                     - 0.3872*cos(2.0*data_angle_rad) + 0.052*sin(2.0*data_angle_rad);
   } else {
      if(lat > 0.0) {
         decl_deg = MAX_DECL_ANGLE;
      } else {
         decl_deg = -MAX_DECL_ANGLE;
      }
   }
   // Convert decl_deg to radians
   decl_rad = -radians(decl_deg);
   // Convert latitude to radians
   lat_rad  = -radians(lat);
   
   cond = 0.133*degrees(acos(-tan(lat_rad)*tan(decl_rad)));
   if(cond != cond) {
      if(fabs(-tan(lat_rad)*tan(decl_rad)) == -tan(lat_rad)*tan(decl_rad)) {
         day_length = 0.0;
      } else {
         day_length = 24.0;
      }
   } else {
      day_length = 0.133*degrees(acos(-tan(lat_rad)*tan(decl_rad)));
   }
   
   return day_length;
}

////////////////////////////////////////////////////////////////////////////////
//! compute_dyl_factor
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double compute_dyl_factor (double lat, double day) {
   double dyl = 0.0;
   double dyl_max = 0.0;
   double factor = 0.0;
   
   // Set 3rd parameter to true to compute maximum day length for a given latitude
   // Formula for computing factor is from
   // http://www.cesm.ucar.edu/models/ccsm4.0/clm/CLM4_Tech_Note.pdf (Page 169)
   dyl = get_day_length(lat, day*30.0 + 15.0, false);
   dyl_max = get_day_length(lat, day*30.0 + 15.0, true);
   
   if(dyl_max == 0.0 or dyl == 0.0) {
      factor = 0.01;
   } else if(dyl > dyl_max){
      factor = 1.0;
   } else {
      factor = dyl/dyl_max;
   }
   
   return factor;
}

////////////////////////////////////////////////////////////////////////////////
//! init_sites
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void init_sites (site** firsts, UserData* data) {
   FILE *outfile;

   int first_site_flag = 1; /* flag for first site */
   int model_site_flag = 0; /* flag to model site  */
   site *new_site, *last_site;

   if(data->cd_file) {
      char filename[STR_LEN];

      strcpy(filename, data->base_filename);
      strcat(filename, ".cd");
      if ( !(outfile = fopen(filename, "a")) ) {
         fprintf(stderr, "init_sites: Can't open cd file: %s \n", filename);
         exit(1);
      }
   }

   printf("intializing sites \n");

   /* loop over region */
   size_t counter = 0;
   for (size_t y=0; y<data->n_lat; y++) {
      for (size_t x=0; x<data->n_lon; x++) {
         data->map[y][x] = NULL;
         counter ++;
         model_site_flag = model_site(y, x, counter, data);
         if (model_site_flag == 1) {
            /* allocate memory for new site */
            new_site = (site *) malloc (sizeof(site));
            if (new_site == NULL) {
               fprintf(stderr,"init_sites: malloc site: out of memory\n");
               exit(1);
            }
            new_site->next_site = NULL;
      
            // assign site attributes 
            new_site->data = data;

            // allocate memory for site data 
            new_site->sdata = new SiteData(y, x, *data);
            if(data->do_downreg) {
               for(size_t z=0;z<N_CLIMATE;z++) {
                  new_site->dyl_factor[z] = compute_dyl_factor(new_site->sdata->lat_, z);
               }
            }
            if(data->cd_file) {
               fprintf(outfile, "new site name: %s\n", new_site->sdata->name_);
            }

            // check to see if site of interest
            new_site->finished = 0;
            new_site->skip_site = 0;
            
            if ( ! new_site->sdata->readSiteData(*data) ) {
               // missing inputs... free memory and move on 
               delete new_site->sdata;
               free(new_site);
               continue;
            }

            new_site->area_burned                   = 0.0;
            new_site->last_site_total_c             = 0.0;

            for (size_t lu=0; lu<N_LANDUSE_TYPES; lu++) {
               new_site->youngest_patch[lu]         = NULL;
               new_site->oldest_patch[lu]           = NULL;
               new_site->new_patch[lu]              = NULL;

               new_site->months_below_rain_crit[lu] = 0.0;
               new_site->fire_flag[lu]              = 0;    /* initialize site fire flag  */
               new_site->fuel[lu]                   = 0.0;  /* initialize site fuel level */
               new_site->last_total_c[lu]           = 0.0;
            }

            new_site->area_fraction[LU_NTRL]              = 1.0;

            if (data->restart) {
               if (data->old_restart_read) {
                  // read in inital patch distribution for the site 
                  read_patch_distribution(&new_site,data);
               } else if (data->new_restart_read) {
                  data->restartReader->readPatchDistribution(new_site, *data);
               } else {
                  fprintf (stderr, "No restart read-type specified\n");
                  exit(1);
               }
            } else {
               data->start_time = 0;
               // create inital patches for the site 
               init_patches(&new_site,data);
            }

#if LANDUSE
            for (size_t i=0; i<2; i++) {
               for (size_t j=0; j<N_SBH_TYPES; j++) {
                  new_site->area_harvested[i][j]          = 0.0;
                  new_site->biomass_harvested[i][j]       = 0.0;
                  new_site->biomass_harvested_unmet[i][j] = 0.0;
               }
            }
            /* if data->start_year is > 0, then we are restarting from *
             * previous landuse... no need to init_landuse_patches     */
            if (data->start_time == 0) {
               init_landuse_patches(&new_site, data);
            }
#ifndef COUPLED
            read_transition_rates(&new_site, data);
#endif
            update_landuse(new_site, *data);
#endif
             
            /* calculate disturbance rates */
             //ml-modified
            //calculate_disturbance_rates(0, &(new_site->oldest_patch[LU_NTRL]), data);
             
            update_site(&new_site, data);

            /* link sites */
            data->map[y][x] = new_site;
            if (first_site_flag == 1) { /* if first site do this, otherwise link it */
               *firsts = new_site;
               new_site->next_site = NULL;
               last_site = new_site;
               first_site_flag = 0; /* reset flag */
            } else {
               new_site->next_site = NULL;
               last_site->next_site = new_site;
               last_site = new_site;
            }

            data->number_of_sites++; /* increment counter */
            new_site->function_calls = 0;
#ifndef USEMPI
            if (data->number_of_sites % 50 == 0) {
               printf("N sites = %d\n",data->number_of_sites);
            }
#endif

         } /* end model_site_flag loop */
      } /* end loops over grid */
   } /* end loops over grid */

   printf("Total N sites = %d\n",data->number_of_sites);
   fprintf(stdout,"Init sites complete\n");
   if(data->cd_file) {
      fclose(outfile);
   }
}


////////////////////////////////////////////////////////////////////////////////
//! model_site
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
int model_site (size_t y, size_t x, size_t counter, UserData* data) {
   /* This function returns 1 if site should be modeled, zero otherwise */

   int model_site_flag = 0;

   if(!data->is_site) {
#if 0
      if (data->mask[y][x] && data->grid_cell_area[y][x] > 0) 
#else
      if (data->grid_cell_area[y][x] > 0) 
#endif
         model_site_flag = 1;  /*land mask*/

      if ((counter % data->rarify_factor)) {
         model_site_flag = 0; /*rarify*/
      }
      if (is_soi(data->lats[y],data->lons[x], *data)) {
         model_site_flag = 1; /*retain sois*/
      }
#if 0
      if(!strcmp(data->region,"US")) {
         if (data->gcode[y][x] != 221)  { /* not us          */
            model_site_flag = 0;       /* dont model site */
         }
      }
#endif
   } else {
      if (data->grid_cell_area[y][x] > 0) {
         if (is_soi(data->lats[y], data->lons[x], *data)) {
            model_site_flag = 1;
         }
      }
   }
  
   return model_site_flag;
}


////////////////////////////////////////////////////////////////////////////////
//! community_dynamics
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void community_dynamics (unsigned int t, double t1, double t2, 
                         site** siteptr, UserData* data) {
  
   FILE* outfile = NULL;
   if(data->cd_file) {
      char filename[STR_LEN];
   
      strcpy(filename, data->base_filename);
      strcat(filename, ".cd");
      if(data->long_cd_file) {
         if( !(outfile = fopen(filename, "a")) ) {
            fprintf(stderr, "cd: Can't open file: %s \n", filename);
            while(getchar() != '\n');
         }
      }
      else {
         if( !(outfile = fopen(filename, "w")) ) {
            fprintf(stderr, "cd: Can't open file: %s \n", filename);
            while(getchar() != '\n');
         }
      }
   } /* if CD_FILE */

   site* currents = *siteptr;
   
#if FTS
   currents->Update_FTS(data->time_period); 
#endif
   if (currents->skip_site) return;

   /* assign current site to site strucutre in UserData structure data *
    * (address of site containing patch to be integrated)              */  
   if(data->cd_file) {
      fprintf(outfile, "Community dynamics for %s %p t=%d t1=%f t2=%f \n",
              currents->sdata->name_, currents, t, t1, t2); 
   }

#ifdef ED
   /*****************************/
   /****** COHORT DYNAMICS ******/
   /*****************************/
   for (size_t lu=0; lu<N_LANDUSE_TYPES; lu++) {
      patch* currentp = currents->youngest_patch[lu];
      if (currentp != NULL) 
         cohort_dynamics(t, t1, t2, &currentp, outfile, data);
      if (currents->skip_site)
         return;
   } 
#endif

   if(data->patch_dynamics) {
       /************************************/
       /******* PATCH DYNAMICS   ***********/
       /************************************/

       /*reinitialized every time step, value calculated in pd*/
       if ( (t % PATCH_FREQ == 0) && (t > 0) ) {
          currents->area_burned = 0.0;
          if(data->do_hurricane) {
             currents->hurricane_litter = 0.0;
          }
       }

       for (size_t lu=0; lu<N_LANDUSE_TYPES; lu++) {
          if ( (lu != LU_CROP) && (currents->youngest_patch[lu] != NULL) ) {
             patch_dynamics(t, &(currents->youngest_patch[lu]), outfile, data);
          }
       }
   } /*PATCH_DYNAMICS*/
   
#if LANDUSE    
   /************************************/
   /******* LANDUSE DYNAMICS ***********/
   /************************************/
   landuse_dynamics(t, &currents, data);
#endif

   if(data->cd_file) {
      fprintf(outfile, "\n");
      fprintf(outfile, "end of community dynamics\n");
      fclose(outfile);
   }
}

#if FTS
////////////////////////////////////////////////////////////////////////////////
//! Update_FTS
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void site::Update_FTS(unsigned int t){
   int shade; //Percent
   int pt, hr;
   int temp_index, hum_index, rad_index;
   double temp, hum, rad;
   int t1, t2;
   double p1, p2;
   int first_interval, last_interval;
   int current_interval, hrs_per_interval;
   double factor;
   first_interval = t*N_CLIMATE_INPUT/N_CLIMATE;
   last_interval = ((t+1)*N_CLIMATE_INPUT)/N_CLIMATE;
   if (t==11){last_interval -=1;}
   hrs_per_interval = 24/CLIMATE_INPUT_INTERVALS;
   double *Anptr, *Anbptr, *Eptr, *Ebptr, *ptr1, *ptr2, *ptr3, *ptr4;
   unsigned int i;
   
   //Init everything to 0
   for (pt=0;pt<2;pt++){
      tf[pt]=0;
      for (shade=0;shade<=120;shade++){
         An[pt][shade] = 0;
         An_shut[pt][shade] = 0;
         E[pt][shade] = 0;
         E_shut[pt][shade] = 0;
      }
   }
   
   //Loop through days and get the hourly value for each input using correct factor 
   //if the day is only partially weighted
   for (current_interval=first_interval;current_interval<=last_interval;current_interval++){
      if (current_interval==first_interval){
         factor = 1.+current_interval-t*N_CLIMATE_INPUT*1./float(N_CLIMATE);
      }
      else if (current_interval==last_interval) {
         if (t <11){factor = (t+1)*N_CLIMATE_INPUT*1./float(N_CLIMATE)-last_interval;}
         else{factor=1;}
      }
      else {factor = 1;}
      factor*=N_CLIMATE*1./(float)N_CLIMATE_INPUT;
      // convert to hourly
      for (hr=0;hr<24;hr++){
         t1 = hr/hrs_per_interval;
         t2 = (t1+1)%CLIMATE_INPUT_INTERVALS;
         p2 = (hr%hrs_per_interval)*1./hrs_per_interval;
         p1 = 1-p2;
         
         temp = p1*sdata->Input_Temperature[current_interval*CLIMATE_INPUT_INTERVALS+t1]+p2*sdata->Input_Temperature[current_interval*CLIMATE_INPUT_INTERVALS+t2];
         hum = p1*sdata->Input_Specific_Humidity[current_interval*CLIMATE_INPUT_INTERVALS+t1]+p2*sdata->Input_Specific_Humidity[current_interval*CLIMATE_INPUT_INTERVALS+t2];
         rad = p1*sdata->Input_Par[current_interval*CLIMATE_INPUT_INTERVALS+t1]+p2*sdata->Input_Par[current_interval*CLIMATE_INPUT_INTERVALS+t2];
         
         // look up indices
         temp_index = (int)(TEMP_TO_INDEX(temp)+.5);
         hum_index = (int)(HUM_TO_INDEX(hum)+.5);
         
         //perform calculations
         for (pt=0;pt<2;pt++){
            if (pt) tf[pt]+=exp(3000.0 * (1.0 / 288.2 - 1.0 / (temp + 273.2)))/(1.0 + exp(0.4 * (10.0 - temp))) * (1.0 + exp(0.4 * (temp - 50.0)))*factor/24.;
            else tf[pt]+=exp(3000.0 * (1.0 / 288.2 - 1.0 / (temp + 273.2)))/(1.0 + exp(0.4 * (5.0 - temp))) * (1.0 + exp(0.4 * (temp - 45.0)))*factor/24.;
            //The following pointers make the code run faster than directly accessing relevant memory
            /*
             Anptr = data->An[pt][temp_index][hum_index];
             Anbptr = data->Anb[pt][temp_index][hum_index];
             Eptr = data->E[pt][temp_index][hum_index];
             Ebptr = data->Eb[pt][temp_index][hum_index];
             ptr1 = An[pt]; ptr2 = An_shut[pt]; ptr3 = E[pt]; ptr4 = E_shut[pt];*/
            
            //Weight indices for light levels
            for (shade=0;shade<121;shade++){ 
               double logval = exp(-1.0*shade/20.0);
               if (shade==120){logval=0;}
               rad_index = (int)(RAD_TO_INDEX(2*rad*logval)+.5);
               An[pt][shade] += data->An[pt][temp_index][hum_index][rad_index]*factor;
               An_shut[pt][shade] += data->Anb[pt][temp_index][hum_index][rad_index]*factor;
               E[pt][shade] += data->E[pt][temp_index][hum_index][rad_index]*factor;
               E_shut[pt][shade] += data->Eb[pt][temp_index][hum_index][rad_index]*factor;
               
            }               
         }
      }
   }
   
   
}
#endif


////////////////////////////////////////////////////////////////////////////////
//! update_site
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void update_site (site** current_site, UserData* data) {  
   /* updating all info before graphics and printing */

   site* cs = *current_site;
  
#if LANDUSE
   update_landuse(cs, *data);
   update_mean_age_secondary_lands(&cs, data);
#endif
  
   /* disturbance rates */
   for (size_t i=0; i<NUM_TRACKS; i++) {
      cs->site_disturbance_rate[i]  = 0.0;
   }
   cs->site_total_disturbance_rate  = 0.0;
   /* biomass */
   cs->site_total_ag_biomass        = 0.0;
   cs->site_total_biomass           = 0.0;
#ifdef ED
   for (size_t i=0; i<NSPECIES; i++) {
      cs->site_total_spp_biomass[i] = 0.0;
      cs->site_total_spp_babove[i]  = 0.0;
      cs->site_basal_area_spp[i]    = 0.0;
   }
   cs->site_basal_area              = 0.0;
   cs->site_lai                     = 0.0;
#endif
   /* soil pools */
   cs->site_total_soil_c            = 0.0;
   cs->site_litter                  = 0.0;
   cs->site_fast_soil_C             = 0.0;
   cs->site_structural_soil_C       = 0.0;
#ifdef ED
   cs->site_slow_soil_C             = 0.0;
   cs->site_passive_soil_C          = 0.0;
   cs->site_mineralized_soil_N      = 0.0;
   cs->site_fast_soil_N             = 0.0;
   cs->site_structural_soil_L       = 0.0;
   cs->site_total_soil_N            = 0.0;
#endif
   /* total c */
   cs->site_total_c                 = 0.0;
   /* c fluxes */
   cs->site_npp                     = 0.0;
   cs->site_nep                     = 0.0;
   cs->site_rh                      = 0.0;
   cs->site_dndt                    = 0.0;
   cs->site_aa_lai                  = 0.0;
   cs->site_aa_npp                  = 0.0;
   cs->site_aa_nep                  = 0.0;
   cs->site_aa_rh                   = 0.0;
#ifdef ED
   cs->site_gpp                     = 0.0;
   cs->site_npp2                    = 0.0;
   cs->site_aa_gpp                  = 0.0;
   cs->site_aa_npp2                 = 0.0;
   /* water */
   cs->site_total_water             = 0.0;
   cs->site_total_theta             = 0.0;
   cs->site_total_water_uptake      = 0.0;
   cs->site_total_water_demand      = 0.0;
   cs->site_total_perc              = 0.0;
   cs->site_total_soil_evap         = 0.0;
   /* height */
   cs->site_avg_height              = 0.0;
#endif

   for (size_t lu=0; lu<N_LANDUSE_TYPES; lu++) {
      update_site_landuse(&cs, lu, data);

      /* disturbance rates */
      for (size_t i=0; i<NUM_TRACKS; i++) {
         cs->site_disturbance_rate[i] += cs->disturbance_rate[i][lu];
         cs->site_total_disturbance_rate += cs->disturbance_rate[i][lu];
      }
      /* biomass */
      cs->site_total_ag_biomass   += cs->total_ag_biomass[lu] * cs->area_fraction[lu];         
      cs->site_total_biomass      += cs->total_biomass[lu] * cs->area_fraction[lu];
#ifdef ED
      cs->site_basal_area         += cs->basal_area[lu] * cs->area_fraction[lu];         
      cs->site_lai                += cs->lai[lu] * cs->area_fraction[lu];         
      for(size_t i=0;i<NSPECIES;i++){
         cs->site_total_spp_biomass[i] += cs->total_spp_biomass[i][lu] * cs->area_fraction[lu];
         cs->site_total_spp_babove[i]  += cs->total_spp_babove[i][lu] * cs->area_fraction[lu];
         cs->site_basal_area_spp[i]    += cs->basal_area_spp[i][lu] * cs->area_fraction[lu];
      }
#endif
      cs->site_total_soil_c       += cs->total_soil_c[lu] * cs->area_fraction[lu];
      cs->site_litter             += cs->litter[lu] * cs->area_fraction[lu];
      cs->site_fast_soil_C        += cs->fast_soil_C[lu] * cs->area_fraction[lu];;
      cs->site_structural_soil_C  += cs->structural_soil_C[lu] * cs->area_fraction[lu];
#ifdef ED
      cs->site_slow_soil_C        += cs->slow_soil_C[lu] * cs->area_fraction[lu];
      cs->site_passive_soil_C     += cs->passive_soil_C[lu] * cs->area_fraction[lu];
      cs->site_mineralized_soil_N += cs->mineralized_soil_N[lu] * cs->area_fraction[lu];
      cs->site_fast_soil_N        += cs->fast_soil_N[lu] * cs->area_fraction[lu];
      cs->site_structural_soil_L  += cs->structural_soil_L[lu] * cs->area_fraction[lu];
#endif
      cs->site_aa_lai             += cs->aa_lai[lu] * cs->area_fraction[lu];
      /* total c */
      cs->site_total_c            += cs->total_c[lu] * cs->area_fraction[lu];         
      /* c fluxes */
      cs->site_npp                += cs->npp[lu] * cs->area_fraction[lu];
      cs->site_nep                += cs->nep[lu] * cs->area_fraction[lu];
      cs->site_rh                 += cs->rh[lu] * cs->area_fraction[lu];      
      cs->site_aa_npp             += cs->aa_npp[lu] * cs->area_fraction[lu];
      cs->site_aa_nep             += cs->aa_nep[lu] * cs->area_fraction[lu];
      cs->site_aa_rh              += cs->aa_rh[lu] * cs->area_fraction[lu];
      cs->site_dndt               += cs->dndt[lu]*cs->area_fraction[lu];
#ifdef ED
      cs->site_gpp                += cs->gpp[lu] * cs->area_fraction[lu];
      cs->site_npp2               += cs->npp2[lu] * cs->area_fraction[lu];
      cs->site_aa_gpp             += cs->aa_gpp[lu] * cs->area_fraction[lu];
      cs->site_aa_npp2            += cs->aa_npp2[lu] * cs->area_fraction[lu];
      if (data->time_period == 0) {
         cs->nep2[lu] = cs->total_c[lu] - cs->last_total_c[lu];
         cs->last_total_c[lu] = cs->total_c[lu];
      }
      /* water */
      cs->site_total_water        += cs->water[lu] * cs->area_fraction[lu]; 
      cs->site_total_theta        += cs->theta[lu] * cs->area_fraction[lu];
      cs->site_total_water_uptake += cs->total_water_uptake[lu] * cs->area_fraction[lu];         
      cs->site_total_water_demand += cs->total_water_demand[lu] * cs->area_fraction[lu];         
      cs->site_total_perc         += cs->perc[lu] * cs->area_fraction[lu];         
      cs->site_total_soil_evap    += cs->soil_evap[lu] * cs->area_fraction[lu];     
      /* height */
      cs->site_avg_height         += cs->lu_avg_height[lu] * cs->area_fraction[lu];
#endif
   }
   if (data->time_period == 0) {
      cs->site_nep2 = cs->site_total_c - cs->last_site_total_c;
      cs->last_site_total_c = cs->site_total_c;
   }
#ifdef ED
   cs->site_total_soil_N = cs->site_mineralized_soil_N + cs->site_fast_soil_N
      + cs->site_structural_soil_C * 1.0 / data->c2n_structural 
      + cs->site_slow_soil_C * 1.0 / data->c2n_slow 
      + cs->site_passive_soil_C * 1.0 / data->c2n_passive; 

   /* update patch profiles */
   patch* cp = cs->oldest_patch[LU_NTRL];
   while (cp != NULL) {
      species_patch_size_profile(&cp, N_DBH_BINS, data);
      cp = cp->younger;
   } /* end loop over patches */
   species_site_size_profile(&cs, N_DBH_BINS, data);
#endif

#if COUPLED
   data->glm_data[1].smb[cs->sdata->globY_][cs->sdata->globX_] = cs->total_ag_biomass[LU_SCND];
#endif
}


////////////////////////////////////////////////////////////////////////////////
//! update_site_landuse
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void update_site_landuse(site** siteptr, size_t lu, UserData* data) {

   site* cs = *siteptr;

   /* disturbance rates */
   for (size_t i=0; i<NUM_TRACKS; i++) {
      cs->disturbance_rate[i][lu]  = 0.0;
   }
   /* biomass */
   cs->total_ag_biomass[lu]        = 0.0;
   cs->total_biomass[lu]           = 0.0;
   /* soil pools */
   cs->total_soil_c[lu]            = 0.0;
   cs->litter[lu]                  = 0.0;
   cs->fast_soil_C[lu]             = 0.0;
   cs->structural_soil_C[lu]       = 0.0;
#ifdef ED
   cs->slow_soil_C[lu]             = 0.0;
   cs->passive_soil_C[lu]          = 0.0;
   cs->mineralized_soil_N[lu]      = 0.0;
   cs->fast_soil_N[lu]             = 0.0;
   cs->structural_soil_L[lu]       = 0.0;
#endif
   cs->aa_lai[lu]                  = 0.0; 
   /* total c */
   cs->total_c[lu]                 = 0.0;
   /* c fluxes */
   cs->npp[lu]                     = 0.0;
   cs->nep[lu]                     = 0.0;
   cs->rh[lu]                      = 0.0;
   cs->aa_npp[lu]                  = 0.0; 
   cs->aa_nep[lu]                  = 0.0; 
   cs->aa_rh[lu]                   = 0.0;
   cs->dndt[lu]                    = 0.0;
#ifdef ED
   cs->npp2[lu]                    = 0.0;
   cs->gpp[lu]                     = 0.0;
   cs->aa_gpp[lu]                  = 0.0;
   cs->aa_npp2[lu]                 = 0.0; 
   /* water */
   cs->water[lu]                   = 0.0;
   cs->theta[lu]                   = 0.0;
   cs->total_water_uptake[lu]      = 0.0;
   cs->total_water_demand[lu]      = 0.0;
   cs->perc[lu]                    = 0.0;
   cs->soil_evap[lu]               = 0.0;        
   /* height */
   cs->lu_avg_height[lu]           = 0.0;
   cs->basal_area[lu]              = 0.0;
   cs->lai[lu]                     = 0.0;
   for (size_t i=0; i<NSPECIES; i++) {
      cs->total_spp_biomass[i][lu] = 0.0;
      cs->total_spp_babove[i][lu]  = 0.0;
      cs->basal_area_spp[i][lu]    = 0.0;
   }
#endif
   
   patch* cp = cs->youngest_patch[lu];
   while (cp != NULL) {
      update_patch(&cp, data);

      if (cs->area_fraction[lu] > data->min_area_fraction) {
         double frac = cp->area / (cs->area_fraction[lu] * data->area);

         /* disturbance rates */
         for (size_t i=0; i< NUM_TRACKS; i++) {
            cs->disturbance_rate[i][lu] += cp->disturbance_rate[i] * frac;
         }
         /* biomass */
         cs->total_ag_biomass[lu]   += cp->total_ag_biomass * frac;
         cs->total_biomass[lu]      += cp->total_biomass * frac;

         /* soil pools */
         cs->total_soil_c[lu]       += cp->total_soil_c * frac;
         cs->litter[lu]             += cp->litter * frac;
         cs->fast_soil_C[lu]        += cp->fast_soil_C * frac;
         cs->structural_soil_C[lu]  += cp->structural_soil_C * frac;
#ifdef ED
         cs->slow_soil_C[lu]        += cp->slow_soil_C * frac;
         cs->passive_soil_C[lu]     += cp->passive_soil_C * frac;
         cs->mineralized_soil_N[lu] += cp->mineralized_soil_N * frac;
         cs->fast_soil_N[lu]        += cp->fast_soil_N * frac;
         cs->structural_soil_L[lu]  += cp->structural_soil_L * frac;
#endif
         cs->aa_lai[lu]             += cp->aa_lai * frac;
         /* total_c */
         cs->total_c[lu]            += cp->total_c * frac;
         /* c fluxes */
         cs->npp[lu]                += cp->npp * frac;
         cs->nep[lu]                += cp->nep * frac;
         cs->rh[lu]                 += cp->rh * frac;
         cs->aa_npp[lu]             += cp->aa_npp * frac;
         cs->aa_nep[lu]             += cp->aa_nep * frac;
         cs->aa_rh[lu]              += cp->aa_rh * frac;
         cs->dndt[lu]               += cp->dndt * frac;
#ifdef ED
         cs->gpp[lu]                += cp->gpp * frac;
         cs->npp2[lu]               += cp->npp2 * frac;
         cs->aa_gpp[lu]             += cp->aa_gpp * frac;
         cs->aa_npp2[lu]            += cp->aa_npp2 * frac;
         /* water */
         cs->water[lu]              += cp->water * frac;
         cs->theta[lu]              += cp->theta * frac;
         /* TODO: why are these different? -justin */
         cs->total_water_uptake[lu] += cp->total_water_uptake / (cs->area_fraction[lu]*data->area);
         cs->total_water_demand[lu] += cp->total_water_demand / (cs->area_fraction[lu]*data->area);
         cs->perc[lu]               += cp->perc * frac;
         cs->soil_evap[lu]          += cp->soil_evap * frac;
         /* height */
         /* insulate against empty patches */
         if (cp->tallest != NULL)
            cs->lu_avg_height[lu]   += cp->tallest->hite * frac;

         cs->basal_area[lu]         += cp->basal_area * frac;
         cs->lai[lu]                += cp->lai * frac;
         for (size_t i=0; i<NSPECIES; i++) {
            cs->total_spp_biomass[i][lu] += cp->total_spp_biomass[i] * frac;
            cs->total_spp_babove[i][lu]  += cp->total_spp_babove[i] * frac;
            cs->basal_area_spp[i][lu]    += cp->basal_area_spp[i] * frac;
         }
#endif /* ED */
      }
      cp = cp->older;
   }
}


#ifdef ED
////////////////////////////////////////////////////////////////////////////////
//! species_site_size_profile
//! binned site size profiles
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void species_site_size_profile (site** pcurrents, unsigned int nbins, 
                                UserData* data) {
   /* sum patch profiles to get species distributions at site level */

   site* cs = *pcurrents;
   for(size_t i=0; i<NSPECIES; i++) {  
      for(size_t j=0; j<N_DBH_BINS; j++) { 
         cs->spp_density_profile[i][j]    = 0.0;
         cs->spp_basal_area_profile[i][j] = 0.0;
         cs->spp_agb_profile[i][j]        = 0.0;
      }
   }

   patch* cp = cs->oldest_patch[LU_NTRL];
   while (cp != NULL) {
      for (size_t i=0; i<NSPECIES; i++) {  
         for (size_t j=0; j<N_DBH_BINS; j++) { 
            cs->spp_density_profile[i][j] += cp->spp_density_profile[i][j] * cp->area / data->area;
            cs->spp_basal_area_profile[i][j] += cp->spp_basal_area_profile[i][j] * cp->area / data->area;
            cs->spp_agb_profile[i][j] += cp->spp_agb_profile[i][j] * cp->area / data->area; 
         }
      }
      cp = cp->younger; 
   } /* end loop over patches */

   /* size profile summed over all species */
   for (size_t j=0; j<N_DBH_BINS; j++) { 
      cs->density_profile[j]    = 0.0;
      cs->basal_area_profile[j] = 0.0;
      cs->agb_profile[j]        = 0.0;
      for (size_t i=0; i<NSPECIES; i++) {
         cs->density_profile[j] += cs->spp_density_profile[i][j];
         cs->basal_area_profile[j] += cs->spp_basal_area_profile[i][j];
         cs->agb_profile[j] += cs->spp_agb_profile[i][j];
      }
   }   
}
#endif /* ED */

