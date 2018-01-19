 #ifndef EDM_DOMAIN_H_
#define EDM_DOMAIN_H_

#include <cstdlib>
#include <vector>
#include <string>

#define MODEL_CONFIG_FILE "models.cfg"
#define PARAMS "params"
#define MODEL_IO     "io"
#define PFTS   "pfts"
#define DOES_COMPILE 0

#define NCERR(s, e) { printf("Error: %s %s\n", s, nc_strerror(e)); exit(62); }
#define ABS(a)    (((a) < 0) ? -(a) : (a))

// THREADING: CHOOSE ONE OR NEITHER OF THE FOLLOWING TWO
// Buggy with landuse
#define GCD 0 ///< Grand Central Dispatch. Works on Mac only
#define TBB 1 ///< Intel Thread Building Blocks

////////////////////////////////////////
//    LAND USE
////////////////////////////////////////
#define LANDUSE 0 ///< Flag to turn on land use dynamics
#define FASTLOAD 1
#define COUPLE_FAR 1
#define INI_Year 1500
#define N_LAI 6
#define WT_Abg_PROFILE 1
const int LAI_INTERVAL[]={0,1.5,5, 10, 20, 30}; //The elements number should be same as N_LAI

#define FTS 0

////////////////////////////////////////
//    INTEGRATION                
////////////////////////////////////////
#if defined ED
#define N_CLIMATE 12              ///< 12 is monthly
#define N_CLIMATE_INPUT 365
#define CLIMATE_INPUT_INTERVALS 4 ///< 4 is 6 hrly, 24 is hrly...
#define N_SUB N_CLIMATE
#define COHORT_FREQ N_CLIMATE/12  ///< Freq of cohort_dynamics in N_SUB units (mnthly)
#define PATCH_FREQ N_CLIMATE
#define LANDUSE_FREQ N_CLIMATE    ///< Freq of land use dynamics in N_SUB units
#elif defined MIAMI_LU
#define N_CLIMATE 1               ///< 12 is monthly 
#define N_SUB 1
#define PATCH_FREQ 1
#define LANDUSE_FREQ 1            ///< Freq of land use dynamics in N_SUB units
#endif // ED V. MIAMI_LU
#define TIMESTEP 1.0/N_SUB        ///< Integration timestep (yrs)

////////////////////////////////////////
//    BIOLOGY/BIOGEOCHEMISTRY     
////////////////////////////////////////
#define NUM_TRACKS 2      ///< Number of disturbance tracks
#define NSPECIES 7        ///< Number of species
#define PT 2              ///< Number of different plant physiologies
#define NCMPT 6           ///< Number carbon compartments in allocation array

////////////////////////////////////////
//    FUSION/FISSION              
////////////////////////////////////////
#define N_DBH_BINS 20     ///< No. of dbh bins used when comparing patches

////////////////////////////////////////
//    MISC                        
////////////////////////////////////////
#define MAX_DECL_ANGLE 23.4667
#define NUM_Vm0s 3                           ///< Number of Vm0 bins

#define LIGHT_RESOLUTION 20                  ///< Resolution of light in lut
#define MIN_LIGHT 6
#define N_LIGHT LIGHT_RESOLUTION*MIN_LIGHT+1 ///< Number light bins in lut

#define KM2_PER_M2 0.000001
#define T_PER_KG 0.001
#define GT_PER_T 0.000000001

#define STR_LEN 512            ///< Was 256

// Landuse
#if LANDUSE
#define N_LANDUSE_YEARS 506    ///< Number of years of land use history data
#define N_LANDUSE_TYPES 4
#define N_VBH_TYPES 2
#define N_SBH_TYPES 3
#else
#define N_LANDUSE_YEARS 1
#define N_LANDUSE_TYPES 1
#endif

// Constants for types of landuse tiles
#define LU_NTRL 0
#define LU_SCND 1 ///< For now, SCND has to be second for sbh/vbh arrays
#define LU_CROP 2
#define LU_PAST 3

////////////////////////////////////////
//    PRINTING                    
////////////////////////////////////////
#ifdef ED
#define PRINTFREQ 1 ///< Timesteps bewteen output, 1=1month
#elif defined MIAMI_LU
#define PRINTFREQ  1        ///< Timesteps between output, 1=1yr
#endif

#if FTS
#define TEMP_TO_INDEX(T)  1.*(T+75.)
#define RAD_TO_INDEX(R)   R/1.
#define HUM_TO_INDEX(H)   H/0.001
#endif


// Forward declarations
class Outputter;
class Restart;
struct site;
namespace libconfig {
   class Config;
}

//! site of interest
struct soi{ 
   char name[STR_LEN];
   double lat;
   double lon;
   soi* next_soi; ///< Pointer to next soi
};


//////////////////////////////////////////////////////////
// Typedef : UserData contains constants and data vectors
//////////////////////////////////////////////////////////
struct UserData {
   // @TODO: only one of these is necessary
   site* first_site;
   site** site_arr;

   Outputter* outputter;
   Restart* restartWriter;
   Restart* restartReader;
   
   const char *model_name; ///< Which model are we running: ED or MLU?
   int allometry_type; 
   
   //! Objects to read the configuration file
   libconfig::Config* model_cfg;
   libconfig::Config* io_cfg_default;
   libconfig::Config* io_cfg_alternate;
   libconfig::Config* params_cfg_default;
   libconfig::Config* params_cfg_alternate;
   libconfig::Config* pfts_cfg_default;
   libconfig::Config* pfts_cfg_alternate;
   
   std::vector<std::string> list_of_pfts;
   size_t nspecies;                        ///< Number of PFTs, can be changed by user
   size_t num_Vm0;                         ///< Number of Vm0 bins (>= 1)
   double Vm0_max[NSPECIES];               ///< Maximum Vm0 corresponding to each PFT
   std::vector<double> Vm0_bins;           ///< List of Vm0s
   const char *Vm0_basepath;               ///< Path containing multiple mech files
   std::vector<std::string> list_c3_files; ///< C3 mech file corresponding to each Vm0 bin
   std::vector<std::string> list_c4_files; ///< C4 mech file corresponding to each Vm0 bin
   const char *title[NSPECIES];            ///< Name of PFT 
   size_t n_years_to_simulate;
   int is_grass[NSPECIES];                ///< Is the PFT a grass?
   int is_tropical[NSPECIES];             ///< Is the PFT tropical?
   int is_drought_deciduous[NSPECIES];    ///< Is the PFT drought deciduous?
   int is_cold_deciduous[NSPECIES];       ///< Is the PFT cold deciduous?
   double initial_density[NSPECIES];      ///< Initial plant density (indivs/m2)
   
   double ref_hgt[NSPECIES];
   double min_hgt[NSPECIES];
   double b1Ht[NSPECIES];
   double b2Ht[NSPECIES];
   double b1Bs[NSPECIES];
   double b2Bs[NSPECIES];
   double b1Bl[NSPECIES];
   double b2Bl[NSPECIES];
   
   const char *output_base_path;
   const char *gridspec;
   const char *soil_file;
   const char *lu_file;
   const char *lu_init_c_file;
   const char *gfedbf_file;

   int is_site;        ///< Set to 1 to select region, 0 for site
   const char *region; ///< E.g.: GLOBAL, AMAZONIA, US, EUS ...
   double latmin;
   double latmax;
   double lonmin;
   double lonmax;
   
   int landuse_bau;
   int landuse_stop;

   int fire_suppression_stop;
   int fire_suppression;
   int fire_off;
   int fire_gfed;
 
   const char *which_mech_to_use;
   const char *climate_file;
   const char *climate_file_avg;
   const char *mech_c3_file;
   const char *mech_c4_file;
   std::vector<std::string> mech_c3_file_avg;  ///< C3 mech file corresponding to each Vm0 bin when no yearly mech available
   std::vector<std::string> mech_c4_file_avg;  //< C4 mech file corresponding to each Vm0 bin when no yearly mech available
#if FTS
   const char *QAIR_FILE;
   const char *TAIR_FILE;
   const char *SW_FILE;
   const char *C3_FILE;
   const char *C4_FILE;
#endif
    
#if COUPLE_FAR
    const char *PREMECH;
    const char *PREMECH_avg;
    
#endif
    
   int single_year;
   int do_yearly_mech;
   int m_int;
   int m_string;
   
   int hurricane;
   const char *hurricane_file;
   int do_hurricane;
   int hurricanetology;
   int hurricane_ramp;
   const char *hurricane_ramp_exp;
   int hurricane_ramp_years;
   int n_hurricane_years; 
   int hurricane_start_year;

   ////////////////////////////////////////
   //    INTEGRATION                
   ////////////////////////////////////////
   double height_threshold_delta;
   double tmax;        ///< Number of years to simulated
   int stiff_light;    ///< 1: Yes to stiff integration of light levels
   int patch_dynamics; ///< Patch dynamics flag, 1=yes to patch dynamics
   int substeps; 
   
   int restart;
   int old_restart_write;
   int old_restart_read;
   const char *old_restart_exp_name;

   int new_restart_write;
   int new_restart_read;
   const char *restart_dir;

   ////////////////////////////////////////
   //    BIOLOGY/BIOGEOCHEMISTRY     
   ////////////////////////////////////////
   double canopy_damage;     ///< Fraction of canopy damaged by small scale dist
   int open_cycles;          ///< Open biogeochemical cycles, 1=yes,0=no
   int water_competition;    ///< Competition for water flag, 1=yes
   int n_competition;        ///< Competition for N flag, 1=yes
   int n_decomp_limitation;  ///< N control of structuralc decomp, 1=yes
   int internal_recruitment; ///< Internal recruitment flag 1=yes
   int external_recruitment; ///< External fixed recruitment flag 1=yes
   int ncmpt;                ///< Number carbon compartments in allocation array
   double hgtmin;            ///< Lowest height value assigned to data->hgt_min
   
   double area;              ///< Set in pde to reasonable value, say 10000.0, for
                             ///< numerics, actual site area is read in, is huge,
                             ///< and can be used for area dependent calcs later
   
   double c2b;               ///< carbon to biomass conversion
   double growth_resp;       ///< 0.333 fraction of npp lost as growth respiration
   double mass_of_water;     ///< @TODO roughly- fix!, (kg/mm/m2)

   ////////////////////////////////////////
   //    FUSION/FISSION              
   ////////////////////////////////////////
   int patch_termination;       ///< terminate small patches and adj area 
   int patch_fusion;            ///< patch fusion flag, 1=yes to patch fusion
   int cohort_fusion;           ///< cohort fusion flag, 1=yes to cohort fusion
   int cohort_fission;          ///< cohort fission flag, 1=yes to cohort fission
   int cohort_termination;      ///< terminate small cohorts 
   double f_area;               ///< 0.01 min area of patch as a fraction of total area
   double btol;                 ///< 0.0001 termination tol. for cohort biomass (kg m^-2)
   double profile_tol;          ///< 0.2 fractional tolerance for patch spp hgt profiles
   double dbhmax;               ///< max dbh value used in hgt profile comparison
   unsigned int n_dbh_bins;     ///< no. of dbh bins used when comparing patches
   double ntol;                 ///< min plant density for hgt bin to be used in height profile comparisons
   double fusetol;              ///< min fractional difference in dbh bet cohorts
   double lai_tol;              ///< maximum LAI allowed for a cohort
   double smallest_new_patch_f; ///< 0.00001 as fraction of data->area
   double min_change_in_area;
   double min_change_in_area_patch;
   double min_area_fraction;
   double min_patch_area;
   double min_landuse_area_fraction;
   double min_landuse_change_fraction;   
   
   ////////////////////////////////////////
   //    MISC                        
   ////////////////////////////////////////
   int do_downreg;        ///< Set to 1 to downregulate Vm0 value based on LAI and day length
   int additional_mort;   ///< Bool to add mortality to all PFTs (den_ind_mort in ED_pft.defaults.cfg)
   int mort_s_hemi;       ///< Use different value for mortality in Southern hemisphere (den_ind_mort_s_hemi in ED_pft.defaults.cfg)
   int hgt_lim_to_repro;  ///< Bool to have minimum height for reproduction (repro_ht_thresh in ED_pft.defaults.cfg)
   double tropic_n_limit; ///< Northern limits to the tropics
   double tropic_s_limit; ///< Southern limits to the tropics
   int rarify_factor;     ///< To coarsen regional run, skip every # of sites

   double cohort_shading; ///< degree of within cohort shading
   double self_shading;   ///< 0.5, degree of self shading
   int n_init_patches; 
   double cell_area;      ///< Used to be data->area/data->n_init_patches

   ////////////////////////////////////////
   //    PRINTING                    
   ////////////////////////////////////////
   int print_output_files;   ///< 1 = yes, do you want all of those output files printed?
                             ///< 1 to print full time dependent files,
                             ///< 0 for small file for most recent time period only,                           
                             ///< Warning: long file could get huge

   int print_system_state;   ///< flag to print system state files
   int print_ss_freq;        ///< in NSUB units
   
   // These two (cd_file, fp_file) have to be off if using GCD or Intel TBB !!!
   int cd_file;              
   int fp_file;       
 
   int long_patch_file;
   int long_cd_file;
   int long_fp_file;
   int long_cohort_file;
        
   ////////////////////////////////////////
   //    ACCOUNTING                  
   ////////////////////////////////////////
   // Minimum sizes used in calc of patch site chars
   double min_dbh_class;         ///< Min dbh (cm) class for dbh calc
   double min_hgt_class;         ///< Min h(m) class used for dbh calc
   
   site ***map;
   unsigned int number_of_sites; ///< Number of sites modeled
   int number_of_sites_hydro;    ///< Number of hydro grid cells modeled, HYDRO
   char base_filename[STR_LEN];  ///< Output file base name
   char outdir[STR_LEN];
   char expname[STR_LEN];
   int climate_file_ncid;        ///< Stores the handle to the climate file to avoid re-opening
   int soil_file_ncid;
   int lu_file_ncid;             ///< Stores the handle to the landuse file to avoid re-opening
   int premech_file_ncid;
   int qair_file_ncid;
   int sw_file_ncid;
   int tair_file_ncid;

   // These are ED only
   int mech_c3_file_ncid[NUM_Vm0s];
   int mech_c4_file_ncid[NUM_Vm0s];
   int mechanism_year;           ///<stores the mechanism year to use
   char mech_year_string[256];

   // Integration structures parameters and coefficients
   double deltat;
   unsigned int time_period;     ///< Within year timestep being integrated
   unsigned int year;            ///< Year of simulation
   int start_time;               ///< Time to start simulation when restarting from landuse

   int hurricane_year;           ///< Year of hurricane record to use
   double *hurricane_ramp_factor;

   // Misc parameters LANDUSE
   double crop_residue;          ///< Fraction of crop harvest left in litter and soil
   double grazing_residue;
   double grazing_intensity;
   double fraction_harvest_left_on_site;
   double fraction_balive_2_fast; ///< Fraction of balive loaded into fast pool
   double forest_definition;      ///< in kg/m2

   // model parameters
   int patch_freq;            ///< Freq of spawning patches
   double gpp_max;            ///< Max gpp per unit carbon kg kg^-1 m^-2 yr^-1
   double k_lite;             ///< Light saturation coefficient
   double L_extinct;           ///< Light extinction per unit leaf biomass
   double Rn_extinct;          ///< Net Radiative flux extinction per unit leaf biomass
   double m2;                  ///< Rate at which mortality declines with bl (per year - units don't match description)
   double m3;                  ///< Mortality at zero bl (dimensionless???)
   double m1;                  ///< Slope of wood density related mortality 1/((g/cm^3)*yr)
   double sd_mort;             ///< Seedling mortality
   double wup_max;             ///< Max water uptake per unit root carbon kg kg^-1 yr^-1
   double water1;              ///< Water uptake/stress parameter 1/(kg(roots)*yr)
   double nitrogen1;           ///< Plant nitrogen uptake/stress parameter
   double nitrogen2;           ///< Microbial nitrogen uptake/stress parameter
   double seed_rain[NSPECIES]; ///< Open recruitment */
   double sapwood_fraction;    ///< Fraction of stem that is living (sapwood)
   double theta_microbes;      ///< Exponent used to shut down nitrogen uptake by microbes
   double theta_plants;        ///< Exponent used to shut down water and nitrogen uptake by plants
 
   double An_max[NSPECIES];                ///< (g/m2/mo)
   double Vm_max[NSPECIES];
   double leaf_life_span[NSPECIES];        ///< (months)
   double specific_leaf_area[NSPECIES];    ///< (m2/kg)
   double Fm_leaf[NSPECIES];               ///< Leaf litter metabolic fraction
   double Fm_sd[NSPECIES];                 ///< Leaf litter metabolic fraction
   double Fm_stem;                         ///< Stem litter metabolic fraction
   double Fm_root;                         ///< Root litter metabolic fraction

   // Disturbance params
   double treefall_max_disturbance_rate_temp;
   double treefall_max_disturbance_rate_trop;
   double fire_max_disturbance_rate;
   double smoke_fraction;             ///< Fraction of biomass lost as smoke
   double fp1;                        ///< Fuel parameter used in fire dist rate calculation
   double selective_harvest_rate;
   double max_patch_age;              ///< Maximum patch age
#ifdef ED
   double treefall_hite_threshold;    ///< Disturbance threshold height (m)
   double treefall_s_ltht[NSPECIES];  ///< Survivorship from disturbance below hite threshold
   double treefall_s_gtht[NSPECIES];  ///< Survivorship from disturbance above hite threshold
   double fire_hite_threshold;
#endif

   // Nitrogen Model coefficients
#ifdef ED
   double NC_perc_coeff;              ///< coefficient for hydrologic loss of N and C
   double fraction_of_GPP_to_Nfixers; ///< payment for N fixers
   int Nfixer[NSPECIES];              ///< 1=nfixer, =non_nfixer
#endif

#ifdef MIAMI_LU
   double agf_biomass; ///< Above ground fraction biomass
   double K1;          ///< Fast decomp rate
   double K2;          ///< Slow decomp rate
   double wood_allocation;
#endif

#ifdef ED
   int phenology[NSPECIES];          ///< Phenology flag:
                                     ///< 0: Evergreen,
                                     ///< 1: Drought-deciduous,
                                     ///< 2: Cold-deciduous
   
   int pt[NSPECIES];                 ///< Photosynthetic physiology of species 1=c3 2=c4
   double m[NSPECIES];               ///< fraction dispersal is global
   double repro_ht_thresh[NSPECIES]; ///< Minimum height threshold for reproduction
   double q[NSPECIES];               ///< ratio of root to leaf biomass (dimensionless)
   double qsw[NSPECIES];             ///< sapwood area per unit leaf biomass (m2/kg)
   double agf_bs;                    ///< fraction of total structural that is above ground
   double rho[NSPECIES];             ///< wood density values g cm^-3 
   double rho_max1;                  ///< maximum wood density,tropics 
   double rho_max2;                  ///< maximum wood density,temperate 
   double alpha[NSPECIES][NCMPT];    ///< decay rates of plant carbon pools (per year per kgC)
   double beta[NSPECIES][NCMPT];     ///< resp rates of plant carbon pools (per year per kgC)
 
   double max_dbh[NSPECIES];         ///< size at which height growth stops (cm)
   double r_fract[NSPECIES];         ///< fraction of excess c going to seed repro
   double c_fract[NSPECIES];         ///< fraction of excess c going to clonal repro
   double hgt_min[NSPECIES];  
   double bl_min[NSPECIES];  
   double bs_min[NSPECIES];  
   double den_ind_mort[NSPECIES];
   double den_ind_mort_s_hemi[NSPECIES];
 
   // CENTRUY 4BOX
   double l2n_stem;
   double l2n_root;
   double l2n_leaf;
   double l2n_sd;
   double c2n_leaf[NSPECIES];    ///< (gC/gN)
   double c2n_recruit[NSPECIES]; ///< (gC/gN)
   double c2n_sd[NSPECIES];      ///< (gC/gN)
   double c2n_stem;              ///< (gC/gN)
   double c2n_root;              ///< (gC/gN)
   double c2n_structural;        ///< (gC/gN)
   double c2n_slow;              ///< (gC/gN)
   double c2n_passive;           ///< (gC/gN)

#endif // ED

   int num_sois;
   soi *first_soi;

#ifdef USEMPI
   size_t mpi_rank;
   size_t mpi_nproc;
#endif

   size_t n_lat;
   size_t n_lon;
   size_t start_lat;
   size_t start_lon;
   double * lats;
   double * lons;

   // These will hold dynamically allocated 2d arrays
   // @TODO: currently they are not being freed -justin
   double ** grid_cell_area;
   double ** grid_cell_area_total;
   double ** wtr_ice_f;
   //int mask[180][360];

   double ** init_height;
   double ** max_pot_height;

   // Initial lu
   double ** init_c;
   double ** init_p;
   double ** init_v;
   double ** init_csc;
   double ** init_psc;
   double ** init_pb;
   // Yannick GFED   
   double ** gfed_bf; ///< declaring pointer to burned fraction data   

   /*Mechanism look up table for FTS ONLY*/
   
#if FTS
   double An[2][136][30][1300];
   double Anb[2][136][30][1300];
   double E[2][136][30][1300];
   double Eb[2][136][30][1300];
#endif
    
#if COUPLE_FAR
    double global_tmp[288][360][720];
    double global_hum[288][360][720];
    double global_swd[288][360][720];
#endif
    
#if FASTLOAD
#if LANDUSE
    float gfl[N_LANDUSE_TYPES][N_LANDUSE_TYPES-1][N_LANDUSE_YEARS][360][720];
    float gfvh[N_VBH_TYPES][N_LANDUSE_YEARS][360][720];
    float gfsh[N_SBH_TYPES][N_LANDUSE_YEARS][360][720];
#endif
    float **soil_depth;
    float **theta_max;
    float **k_sat;
    float **tau;
    
    float ***climate_temp;
    float ***climate_precip;
    float ***climate_soil;
    float *light_levels;
    
#ifndef COUPLE_FAR
    float *****tf;
    float ******An;
    float ******Anb;
    float ******E;
    float ******Eb;
#endif
    
#endif

#ifdef COUPLED
   double lastTotalC;
   double backupTotalC;
   struct new_data* glm_data; ///< This is defined in glm_coupler.h in the glm src
   int start_year;
   site* sitelist_copy;
#endif

};

UserData* ed_initialize(char *name, const char* cfgFile);
void ed_step (int year, UserData& data);
void ed_finalize(UserData& data);

void** malloc_2d (size_t nrows, size_t ncols, int elementsize);
void init_data(const char* cfgFile, UserData* data);
void init_mech_table (UserData *data);

#endif // EDM_DOMAIN_H_
