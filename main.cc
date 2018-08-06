#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <cerrno>
#include <sys/stat.h>
#include "netcdf.h"

#include "edmodels.h"
#if GCD
#include <dispatch/dispatch.h> 
#elif TBB
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#endif
#ifdef USEMPI
#include "mpi.h"
#include "edmpi.h"
#endif

#include "site.h"
#include "restart.h"
#include "read_site_data.h"
#include "print_output.h"
#include "readconfiguration.h"

time_t seconds;           /* time variable for rnd seeding */
long intdum;              /* random no. seed */

/**************************** Function Prototypes *****************************/
void setup_dirs (UserData& data, char* basename);
void model (UserData& data);
/******************************************************************************/

using namespace std;

#if TBB
using namespace tbb;

class UpdateSites {
   site** const my_site_arr;
   unsigned int t;
   double t1;
   double t2;
   UserData *data;
 public:
   void operator() ( const blocked_range<size_t>& r ) const {
      site** site_arr = my_site_arr;
      for (size_t i=r.begin(); i!=r.end(); ++i) {
         community_dynamics(t, t1, t2, &(site_arr[i]), data); 
         update_site(&(site_arr[i]), data); 
      }
   }
   UpdateSites (site* site_arr[], unsigned int t, double t1, double t2, UserData *data) 
      : my_site_arr(site_arr), t(t), t1(t1), t2(t2), data(data)
   {}
};
#endif


////////////////////////////////////////////////////////////////////////////////
//! ed_initialize
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
UserData* ed_initialize (char* expName, const char* cfgFile) {
   /* Below added to supply yearly lists */
   FILE *namefile;

   /* allocate data storage structure */
   //UserData* data = (UserData*) malloc(sizeof *data);  /* Allocate data memory */
   struct UserData* data = new UserData;
   if (data == NULL) {
      fprintf(stderr, "main: out of memory - can't allocate UserData\n");
      exit(1);;
   }

#ifdef USEMPI
   init_mpi(*data);
#endif

   init_data(cfgFile, data);

   setup_dirs(*data, expName);

#if COUPLED
   data->sitelist_copy = NULL;
#endif

   /* initialize site structures */
   if(data->do_yearly_mech) {
   data->mechanism_year = INI_Year+data->start_time;  /*  Added to allow initial mech year  */
      if(data->m_string) {
         /* Added to read list of years we want */
         namefile = fopen("/Network/Xgrid/data/MSTMIP/model_driver/cru_ncep/file_lists/fl1.txt","r");
         fscanf(namefile,"%s",data->mech_year_string);
      }
   }

   read_input_data_layers(data);

   loadGlobalEnvironmentData(data);
   loadPREMECH(data);
    if (!data->fire_off && data->fire_gfed)
    {
        //read_gfed_bf(data);
        load_GFED(data);
    }
    
#if LANDUSE
   loabGlobalLUData(data);
    loadCropCalendar(data);
#endif


   /* setup netcdf output */
   data->outputter = new Outputter(data);
   registerOutputVars(data->outputter);

   string restartDir (data->outdir);
#if 1-LOCAL_MACHINE
   if (data->new_restart_write) {
      data->restartWriter = new Restart(restartDir + "/RESTART/");
      //data->restartWriter = new Restart("/tmp/");
   }
   if (data->restart && data->new_restart_read) {
      data->restartReader = new Restart(data->restart_dir);
   }
#endif
   site* first_site = NULL;
   init_sites(&first_site, data);

   data->first_site = first_site;

   if (first_site == NULL) {
      fprintf(stderr, "error no valid sites \n"); 
      // TODO: if an mpi processor ends up here, bad things happen
      exit(0);
   }

   if(data->print_output_files) {
      print_initial(first_site, data);
   }

#if TBB
   task_scheduler_init init(500);
#endif

#if GCD || TBB 
   // TODO: this should replace site list
   data->site_arr = (struct site**) malloc (data->number_of_sites 
                                            * sizeof(struct site*));
                                    
   site *siteptr = data->first_site;
   for (size_t i=0; i<data->number_of_sites; i++) {
      data->site_arr[i] = siteptr;
      siteptr = siteptr->next_site;
   }
#endif

#ifdef COUPLED
   data->lastTotalC = 0.0;
   site *cs = first_site;
   while (cs != NULL) {
      data->lastTotalC += cs->site_total_c * cs->sdata->grid_cell_area * T_PER_KG * GT_PER_T;
      cs = cs->next_site;
   }
#endif

   return data;
}

////////////////////////////////////////////////////////////////////////////////
//! ed_finalize
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void ed_finalize(UserData& data) {
   printf("Problematic Sites:\n");
   int count1 = 0, count2 = 0;
   site* current_site = data.first_site;
   while (current_site!=NULL){
      if (current_site->skip_site){
        printf("%s\n", current_site->sdata->name_); 
         count1++;
      }
      count2++;
      current_site = current_site->next_site;
   }
   printf("Skipped %d out of %d sites\n", count1, count2);
   printf("*** Program Complete ***\n");

   // Free up all used memory
   free_user_data(&data);

#ifdef USEMPI
   MPI::Finalize();
#endif
}

#ifdef MAIN
////////////////////////////////////////////////////////////////////////////////
//! main
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
int main (int ac, char *av[]) {

   if (ac == 1) { /* user did not supply experiment name */
      fprintf(stderr, "Usage: %s experiment-name [config-file-name]\n", av[0]);
      return 1;
   }

   UserData* data = NULL;
   if (ac > 2) { // config file was specified on command line
      data = ed_initialize(av[1], av[2]);
   } else {
      data = ed_initialize(av[1], NULL);
   }
   model(*data);
   ed_finalize(*data);
}
#endif


////////////////////////////////////////////////////////////////////////////////
//! setup_dirs
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void setup_dirs(UserData& data, char* basename) {

   if(data.output_base_path) {
      /* if data.output_base_path defined and arg is not absolute path *
      * prepent data.output_base_path to basename */
      if (strchr(basename,'/') != basename) {
         char tmp[STR_LEN];
         strcpy(tmp, data.output_base_path);
         strcat(tmp, basename);
         strcpy(basename, tmp);
      }
   }
   
   int rv = mkdir(basename, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   if ( (rv != 0) && (errno != EEXIST) ) {
      printf("FAILED to create output dir: %s\n",basename);
      exit(-1);
   }

   char restart_dir[STR_LEN];
   strcpy(restart_dir,basename);
   strcat(restart_dir,"/RESTART");
   rv = mkdir(restart_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
   if ( (rv != 0) && (errno != EEXIST) ) {
      printf("FAILED to create restart dir: %s\n",restart_dir);
      exit(-1);
   }

   strcpy(data.outdir, basename);
   strcpy(data.base_filename, basename);

   char* tok = strtok(basename, "/");
   while(tok != NULL) {
      strcpy(data.expname,tok);
      tok = strtok(NULL, "/");
   }

   strcat(data.base_filename, "/");
   strcat(data.base_filename, data.expname);

#ifdef USEMPI
   sprintf(data.base_filename, "%s_%ld", data.base_filename, data.mpi_rank);
#endif
}


////////////////////////////////////////////////////////////////////////////////
//! model
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void model (UserData& data) {

#if USEMPI
   MPI::COMM_WORLD.Barrier();
   unsigned long n_sites;
   MPI::COMM_WORLD.Reduce(&data.number_of_sites, &n_sites, 1, MPI::UNSIGNED_LONG, MPI::SUM, 0);
   if (data.mpi_rank == 0) {
      printf("TOTAL NUMBER OF SITES: %ld\n", n_sites);
   }
   MPI::COMM_WORLD.Barrier();
#endif

   
   printf("****** running model \n"); 

   unsigned int tsteps = ((int)(data.tmax * N_SUB)) + 1;
   for (unsigned int t=data.start_time; t<tsteps; t++) { /* absolute time offset */
      if(data.print_output_files) {
          if (tsteps-t<35*N_SUB+1)
          {
              print_region_files(t,&data.first_site,&data);
          }
      }

      double t1 = t * TIMESTEP;
      double t2 = (t + 1) * TIMESTEP;
      data.year = (int) floor(t * 1.0 / N_CLIMATE);

      if (t % N_CLIMATE == 0) {
#ifdef USEMPI
         mpi_collect_data (data);
         if (data.mpi_rank == 0) {
            printf(" Year: %d\n", data.year);
         }
#else
         printf(" Year: %d\n", data.year);
#endif         
      }

      if ( (t > 0) && (t%data.print_ss_freq == 0) ) {
         if (data.old_restart_write) {
            print_system_states(t, data.first_site, &data);
         } else if (data.new_restart_write) {
            data.restartWriter->storeStates(data.first_site, data.year);
         }
      }

      data.time_period = ((int) rint(t1 * N_CLIMATE)) % N_CLIMATE;
#if USEMPI
      if (data.mpi_rank == 0) {
         printf("TIME PERIOD: %d\n", data.time_period);
      }
#else
      printf("TIME PERIOD: %d\n", data.time_period );
#endif         

      if(data.do_hurricane) {
         if (t%12 == 0) {
            if ( data.hurricanetology ) 
               data.hurricane_year = -1;
            else if ( data.hurricane_ramp ) {
               data.hurricane_year = (int)( rand() / ( ( (double)RAND_MAX + 1 ) / data.n_hurricane_years ) ); 
               printf("Hurricane Year to use: %d\n", data.hurricane_year);
            }
            else
               data.hurricane_year = data.year - data.hurricane_start_year;
         }
      }

      // do_yearly_mech is deprecated in favor of FTS
#if 1
      FILE *namefile;
      if(data.do_yearly_mech)
      {
          if(data.m_int)
          {
             if (t > 0 && t%12 == 0)
             {
                 size_t i=0;
                 data.mechanism_year = INI_Year+t1;
                 if (data.mechanism_year<1850)
                 {
                     data.MERRA2_LUT=0;
                 }
                 else
                 {
                     data.MERRA2_LUT=0;
                 }
                 printf("Mechanism_year_to use: %d\n" , data.mechanism_year);
            
                 if (data.mechanism_year>1980) //loading data for everyear if COUPLE_MERRA2
                 {
                     for (; i< data.num_Vm0;i++)
                     {
                         if (data.mech_c3_file_ncid[i]>0)    ncclose(data.mech_c3_file_ncid[i]);
                         if (data.mech_c4_file_ncid[i]>0)   ncclose(data.mech_c4_file_ncid[i]);
                         
                         data.mech_c3_file_ncid[i] =0;
                         data.mech_c4_file_ncid[i] =0;
                         
                     }
                     if (data.climate_file_ncid>0)  ncclose(data.climate_file_ncid);
                     data.climate_file_ncid =0;

                     printf("Start load Envir\n");
                     loadGlobalEnvironmentData(&data);

                     printf("Start load pre_mech\n");
                     loadPREMECH(&data);
                     
                     if (!data.fire_off && data.fire_gfed)
                     {
                         printf("Start load GBED Burned Fraction data\n");
                         //read_gfed_bf(data);  // Instead, using load_GFED which is defined in read_site_data.cc
                         load_GFED(&data);
                     }

                     site* siteptr = data.first_site;
                     while (siteptr != NULL)
                     {
                         /* Now we have to read the site data again */
                         siteptr->sdata->readSiteData(data);
                         siteptr = siteptr->next_site;
                     }
                 }
             }
          }

          if(data.m_string) {
             if (t > 0 && t%12 == 0) {
                 size_t i=0;
                 for (; i< data.num_Vm0;i++) {
                     ncclose(data.mech_c3_file_ncid[i]);
                     ncclose(data.mech_c4_file_ncid[i]);
                     ncclose(data.climate_file_ncid);
                     
                     
                     data.mech_c3_file_ncid[i] =0;
                     data.mech_c4_file_ncid[i] =0;
                     data.climate_file_ncid =0;
                 }
                //ncclose(data.mech_c3_file_ncid);
                //ncclose(data.mech_c4_file_ncid);
                //ncclose(data.climate_file_ncid);
                //data.mech_c3_file_ncid =0;
                //data.mech_c4_file_ncid =0;
                //data.climate_file_ncid =0;
                fscanf(namefile,"%s",data.mech_year_string);
                printf("Mechanism_year_to use: %s\n" , data.mech_year_string);
                if (strlen(data.mech_year_string)!= 4){
                   fclose(namefile);
                   namefile = fopen("/Network/Xgrid/data/MSTMIP/model_driver/cru_ncep/file_lists/fl1.txt","r");
                   fscanf(namefile,"%s",data.mech_year_string);
                   printf("Mechanism_year_to use: %s\n" , data.mech_year_string);
                }
                 
         
                 freeGlobalEnvironmentData(&data);
                 loadGlobalEnvironmentData(&data);
                 
                site* siteptr = data.first_site;
                while (siteptr != NULL) {
                   /* Now we have to read the site data again */
                   siteptr->sdata->readSiteData(data); 
                   siteptr = siteptr->next_site;
                }    
             }
          }
      }
#endif
#if GCD 
      dispatch_apply(data.number_of_sites, dispatch_get_global_queue(0,0), ^(size_t i) { 
         community_dynamics(t, t1, t2, &(site_arr[i]), &data); 
         update_site(&(site_arr[i]), &data); 
      }); 
#elif TBB
      parallel_for(blocked_range<size_t>(0,data.number_of_sites,3),
                   UpdateSites(data.site_arr, t, t1, t2, &data) );
#endif
      site* siteptr = data.first_site;
      while (siteptr != NULL) {
#if !GCD && !TBB
          
#if CHECK_C_CONSERVE
          ///CarbonConserve
          update_site(&siteptr, &data);
          site* currents = siteptr;
          /// Variables starting with all means total of the entire site. before means results befor implementing a process
          double all_tb_before =currents->site_total_biomass, all_sc_before=currents->site_total_soil_c, all_tc_before =currents->site_total_c;
          double all_f_harv_before = currents->site_forest_harvest, all_p_harv_before = currents->site_pasture_harvest, all_c_harv_before = currents->site_crop_harvest;
          double all_tb_after = 0.0, all_sc_after = 0.0, all_tc_after = 0.0, all_fire_emission_after = 0.0;
          double all_f_harv_after = 0.0, all_p_harv_after = 0.0, all_c_harv_after = 0.0;
          double esti_dt_tc = 0.0, actual_dt_tc= 0.0;
#endif
          community_dynamics(t, t1, t2, &siteptr, &data);
          update_site(&siteptr, &data);
          
#if CHECK_C_CONSERVE
          ///CarbonConserve
          all_tb_after = currents->site_total_biomass;
          all_sc_after = currents->site_total_soil_c;
          all_tc_after = currents->site_total_c;
          all_fire_emission_after = currents->site_fire_emission;
          all_f_harv_after = currents->site_forest_harvest;
          all_p_harv_after = currents->site_pasture_harvest;
          all_c_harv_after = currents->site_crop_harvest;
          actual_dt_tc = all_tc_after- all_tc_before;
          esti_dt_tc = (currents->site_npp_avg-currents->site_rh_avg-currents->site_fire_emission-all_f_harv_after-all_p_harv_after-all_c_harv_after)*data.deltat;
          if (abs(actual_dt_tc-esti_dt_tc)>1e-9)
          {
              printf("Carbon leakage in community_dynamics: imbalance  %.15f actual_dt_tc %.15f esti_dt_tc %.15f\n",actual_dt_tc-esti_dt_tc,actual_dt_tc,esti_dt_tc);
              printf("                                    : site_tc_bf %.15f site_sc_bf   %.15f site_tb_bf %.15f site_f_harc_bf %.15f site_p_harv_bf %.15f site_c_harv_bf %.15f\n",all_tc_before,all_sc_before,all_tb_before,all_f_harv_before,all_p_harv_before,all_c_harv_before);
              printf("                                    : site_tc_af %.15f site_sc_af   %.15f site_tb_af %.15f site_f_harv_af %.15f site_p_harv_af %.15f site_c_harv_af %.15f\n",all_tc_after,all_sc_after,all_tb_after,all_f_harv_after,all_p_harv_after,all_c_harv_after);
              printf("                                    : site_npp  %.15f site_rh   %.15f\n",currents->site_npp_avg,currents->site_rh_avg);
              printf("                                    : site_lat %.3f site_lon %.3f\n",currents->sdata->lat_,currents->sdata->lon_);
              printf(" --------------------------------------------------------------------------------------\n");
          }
#endif
          
#endif
         
         if(data.print_output_files) {
            print_soi_files(t, &siteptr, &data);
         }
         siteptr = siteptr->next_site;
      }
   }
    
}


// TODO: this should just become the default and "model" loop rewritten
#ifdef COUPLED
////////////////////////////////////////////////////////////////////////////////
//! ed_step
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void ed_step (int year, UserData& data) {

   printf("****** ed step: year %d \n", year); 

   unsigned int t = (year - data.start_year) * N_CLIMATE;
   for (int tt=0; tt<N_CLIMATE; tt++) {
      data.time_period = tt;
      double t1 = (t+tt) * TIMESTEP;
      double t2 = (t+tt+1) * TIMESTEP;
      data.year = year;

      printf("TIME PERIOD: %d\n", data.time_period);

      if ( (t > 0) && (t%PRINT_SS_FREQ == 0) ) {
         if (OLD_RESTART_WRITE) {
            print_system_states(t, data.first_site, &data);
         } else if (NEW_RESTART_WRITE) {
            data.restartWriter->storeStates(data.first_site, year);
         }
      }

#if DO_HURRICANE
      if (tt%N_CLIMATE == 0) {
         if ( HURRICANETOLOGY ) 
            data.hurricane_year = -1;
         else if ( HURRICANE_RAMP ) {
            data.hurricane_year = (int)( rand() / ( ( (double)RAND_MAX + 1 ) / N_HURRICANE_YEARS ) ); 
            printf("Hurricane Year to use: %d\n", data.hurricane_year);
         }
         else
            data.hurricane_year = data.year;
      }
#endif

#if GCD 
      dispatch_apply(data.number_of_sites, dispatch_get_global_queue(0,0), ^(size_t i) { 
         community_dynamics(t, t1, t2, &(data.site_arr[i]), &data); 
         update_site(&(data.site_arr[i]), &data); 
      }); 
#elif TBB
      parallel_for(blocked_range<size_t>(0,data.number_of_sites,100), 
                   UpdateSites(data.site_arr, t, t1, t2, &data) );
#endif
      site* siteptr = data.first_site;
      while (siteptr != NULL) {
#if !GCD && !TBB
         community_dynamics(t, t1, t2, &siteptr, &data);
         update_site(&siteptr, &data);
#endif
#if PRINT_OUTPUT_FILES
         print_soi_files(t, &siteptr, &data);
#endif
         siteptr = siteptr->next_site;
      }

#if PRINT_OUTPUT_FILES          
      print_region_files(t,&(data.first_site),&data);
#endif
      t++;
   }
   printf("finished ed_step\n");
}
#endif // COUPLED
