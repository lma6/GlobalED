#ifndef EDM_READ_SITE_DATA_H_
#define EDM_READ_SITE_DATA_H_

#include "edmodels.h"

////////////////////////////////////////
//    SiteData contains site-related 
//    environmental data
////////////////////////////////////////
struct SiteData {
   size_t x_;
   size_t y_;

   size_t globX_;
   size_t globY_;

   double lat_;
   double lon_;
   
   char name_[STR_LEN];

   bool soi;

   double grid_cell_area;            ///< land area of grid cell  (m^2)
   double grid_cell_area_total;      ///< total area of grid cell (m^2)

#ifdef ED
   double temp[N_CLIMATE];           ///< monthly mean temperature (deg C)
   double precip[N_CLIMATE];         ///< monthly rate of rainfall (mm/yr)
   double soil_temp[N_CLIMATE];      ///< monthly mean soil temp (deg C)
   double pet[N_CLIMATE];            ///< potential evapotranspiration (mm/yr) 
#endif // ED
   double precip_average;            ///< avg. yrly precip (mm/yr) 
   double temp_average;              ///< average annual temperature (deg C)
   double soil_temp_average;         ///< annual average soil temp (deg C)
   double pet_average;

   double soil_depth;                ///< depth of soil (mm)
   double soil_evap_conductivity;    ///< soil conductivity for evap water loss

   double dryness_index[N_CLIMATE];  ///< used in fire model */
   double dryness_index_avg;         ///< used in fire model 

   double loss_fraction[NUM_TRACKS]; ///< loss fractions following disturbance

#ifdef MIAMI_LU
   double miami_npp;                 ///< annual net primary production (Kg / (m^2 yr)) 
#endif

#ifdef ED
   double N_conc_in_rain;            ///< concentration of N in rain water 
   double L_top;                     ///< Light at the top of the canopy 
   double Rn_top;                    ///< Net Radiative flux at top of canopy

   double theta_max;                 ///< volume of water in saturated soil (mm/mm^3)
   double k_sat;                     ///< conductivity of saturated soil (yr^-1) 
   double tau;                       ///< conductivity exponent 
    
   // potential photosynthesis and evap 
   double light_levels[PT][NUM_Vm0s][N_LIGHT];
   double tf[PT][NUM_Vm0s][N_CLIMATE];           ///< value of temp function (used in resp calc) (Dimensionless???)

#if FTS
   bool readFTSdata (UserData& data);
   double Input_Temperature[N_CLIMATE_INPUT*CLIMATE_INPUT_INTERVALS];
   double Input_Specific_Humidity[N_CLIMATE_INPUT*CLIMATE_INPUT_INTERVALS];
   double Input_Par[N_CLIMATE_INPUT*CLIMATE_INPUT_INTERVALS];
#else
   double An[PT][NUM_Vm0s][N_CLIMATE][N_LIGHT];  ///< pot. net photosynthesis (gC/(m2 mo)) 
   double E[PT][NUM_Vm0s][N_CLIMATE][N_LIGHT];   ///< potential transpitation (gW/(m2 mo))
   double Anb[PT][NUM_Vm0s][N_CLIMATE][N_LIGHT]; ///< pot. psyn when shut (g/(m2 mo))
   double Eb[PT][NUM_Vm0s][N_CLIMATE][N_LIGHT];  ///< pot. transp when shut (gW/(m2 mo)) 
#endif
   



   // N.A. stuff 
   double first_frost;
   double last_frost;
   double growing_season_length;
#endif // ED 

#if LANDUSE
   double beta[N_LANDUSE_TYPES][N_LANDUSE_TYPES-1][N_LANDUSE_YEARS];

   double vbh[N_VBH_TYPES][N_LANDUSE_YEARS];
   double sbh[N_SBH_TYPES][N_LANDUSE_YEARS];
#endif // LANDUSE 
  
   double avg_hurricane_disturbance_rate;
   double *hurricane_disturbance_rate;
    
#if 0 // Disabled temporarily until fixed
   double ph[N_LANDUSE_YEARS]; ///< prob_harvest(A) (1/yr) 
#endif
    
#if 0 // Disabled temporarily until fixed (if worth it)
   // business as usual 
   double beta_vs_bau;         ///< transistion rates (fr_total_area/yr)
   double beta_vp_bau;            
   double beta_vc_bau;
   double beta_sp_bau;
   double beta_sc_bau;
   double beta_cs_bau;
   double beta_ps_bau;
   double beta_gc_bau;
   double beta_gp_bau;
   double beta_pg_bau;
   double beta_cg_bau;
   double beta_cp_bau;
   double beta_pc_bau;
#endif

   SiteData (size_t y, size_t x, UserData& data);
   bool readSiteData (UserData& data);

 private:

   bool readEnvironmentalData (UserData& data);
#ifdef ED
   bool readMechanismLUT (UserData& data);
   void calcSiteDrynessIndex (UserData& data);
   double calcPETMonthly (size_t month);
#endif
   void calcPETAverage ();
};


bool is_soi(double lat,double lon, UserData& data); ///< flag sites of interest
size_t read_input_data_layers (UserData* data);
int read_hurricane_disturbance (site** siteptr, UserData* data);

#endif // EDM_READ_SITE_DATA_H_ 
