#include <cstdio>
#include <cmath>
#include <cstring>
#include "netcdf.h"

#include "edmodels.h"
#include "site.h"
#include "miami.h"
#if LANDUSE
#include "landuse.h"
#endif
#if USEMPI
#include "edmpi.h"
#endif
#include "fire.h"
#include "read_site_data.h"


size_t read_gridspec (UserData* data);
void read_sois (UserData* data);
////////////////////////////////////////////////////////////////////////////////
//! read_input_data_layers
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
size_t read_input_data_layers (UserData* data) {

   read_sois(data);
   size_t nPotentialSites = read_gridspec(data);
#if LANDUSE
#ifndef COUPLED
   read_initial_landuse_fractions(data);
#endif
   read_initial_landuse_c(data);
#endif

   printf("gfed bool= %d\n", data->fire_gfed);
   printf("fre_off bool= %d\n", data->fire_off);
   if(data->fire_gfed){
      printf("In read_input_data_layer(), fire_gfed bool is true... \n");
      read_gfed_bf(data);
   }	
   
   return nPotentialSites;
}


////////////////////////////////////////////////////////////////////////////////
//! read_sois
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void read_sois (UserData* data) {
   soi *new_soi, *current_soi;
   int i;
   double lat, lon;
   char name[STR_LEN];
   char filename[STR_LEN];
   FILE *infile;

  
   strcpy(filename, "sois.txt");
   infile = fopen(filename,"r");

   printf("read_sois... %s\n", filename);

   if(fscanf(infile,"%*s %d\n",&(data->num_sois)) != 1){
      printf("rmi: Error reading from file: %p\n", infile);
      fclose(infile);
      return;
   } 

   printf("num_sois= %d\n", data->num_sois);

   if(data->num_sois > 0){
      data->first_soi=(soi *) malloc (sizeof(soi));
      fscanf(infile, "%s %lf %lf\n", name, &lat, &lon);
      data->first_soi->lat = lat;
      data->first_soi->lon = lon;
      strcpy(data->first_soi->name, name);
      data->first_soi->next_soi = NULL;

      if(data->num_sois > 1){
         current_soi=data->first_soi;
         for(i=1;i<data->num_sois;i++){
            fscanf(infile, "%s %lf %lf\n", name, &lat, &lon);
            new_soi=(soi*) malloc (sizeof(soi));
            new_soi->lat = lat;
            new_soi->lon = lon;
            strcpy(new_soi->name, name);
            new_soi->next_soi = NULL;
            current_soi->next_soi = new_soi;
            current_soi = new_soi;
         }
      }      
   }
   fclose(infile);
 
   current_soi = data->first_soi;
   while (current_soi != NULL) {
      printf("soi = %s lat= %lf lon= %lf\n",
             current_soi->name, current_soi->lat, current_soi->lon);
      current_soi = current_soi->next_soi;
   }

   printf("read_sois COMPLETE\n");
   return;
}


////////////////////////////////////////////////////////////////////////////////
//! is_soi
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
bool is_soi (double lat, double lon, UserData& data) { /* flag site of interest */

   double dlat = fabs(data.lats[0] - data.lats[1]);
   double dlon = fabs(data.lons[0] - data.lons[1]);
   double lat_min = lat - dlat/2.0;
   double lat_max = lat + dlat/2.0;
   double lon_min = lon - dlon/2.0;
   double lon_max = lon + dlon/2.0;


   bool flag = false;
   soi* cur_soi = data.first_soi;
   while (cur_soi != NULL) {
      if( (cur_soi->lat >= lat_min) && (cur_soi->lat < lat_max)
         && (cur_soi->lon >= lon_min) && (cur_soi->lon < lon_max) )
         flag = true;
      cur_soi = cur_soi->next_soi;
   }
   return flag;   
}


////////////////////////////////////////////////////////////////////////////////
//! get_grid_bounds
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void get_grid_bounds (int ncid, UserData* data) {
   int rv, latid, lonid, varid;
   double *tmp;
   size_t i, start, end, count;
   int in_region_flag;

   /* TODO: this allocates memory that is never freed.
      shouldn't be a problem because there is only ever 
      one 'data' object created, but should fix anyway -- justin */
   
   if ((rv = nc_inq_dimid(ncid, "lat", &latid)))
      NCERR("lat", rv);
   if ((rv = nc_inq_dimlen(ncid, latid, &count)))
      NCERR("lat", rv);
   tmp = (double *)malloc(count * sizeof(double));
   if ((rv = nc_inq_varid(ncid, "lat", &varid)))
      NCERR("lat", rv);
   if ((rv = nc_get_var_double(ncid, varid, &tmp[0])))
      NCERR("lat", rv);

   in_region_flag = 0;
   for (i=0; i<count; i++) {
      if ((!in_region_flag) && (tmp[i] <= data->latmax) && (tmp[i] >= data->latmin)) {
         start = i; 
         in_region_flag = 1;
      }
      if ((in_region_flag) && ((tmp[i] > data->latmax) || (tmp[i] < data->latmin))) {
         end = i; 
         break;
      } else if (i == count - 1) {
         end = i + 1;
      }
   }
   
   data->start_lat = start;
   data->n_lat = end - start;
   data->lats = (double *)malloc(data->n_lat * sizeof(double));
   for (i=start; i<end; i++) {
      data->lats[i-start] = tmp[i];
   }
   free(tmp);

   if ((rv = nc_inq_dimid(ncid, "lon", &lonid)))
      NCERR("lon", rv);
   if ((rv = nc_inq_dimlen(ncid, lonid, &count)))
      NCERR("lon", rv);
   tmp = (double *)malloc(count * sizeof(double));
   if ((rv = nc_inq_varid(ncid, "lon", &varid)))
      NCERR("lon", rv);
   if ((rv = nc_get_var_double(ncid, varid, &tmp[0])))
      NCERR("lon", rv);
   
   in_region_flag = 0;
   for (i=0; i<count; i++) {
      if ((!in_region_flag) && (tmp[i] <= data->lonmax) && (tmp[i] >= data->lonmin)) {
         start = i; 
         in_region_flag = 1;
      }
      if ((in_region_flag) && ((tmp[i] > data->lonmax) || (tmp[i] < data->lonmin))) {
         end = i; 
         break;
      } else if (i == count - 1) {
         end = i + 1;
      }
   }
   
   data->start_lon = start;
   data->n_lon = end - start;
   data->lons = (double *)malloc(data->n_lon * sizeof(double));
   for (i=start; i<end; i++) {
      data->lons[i-start] = tmp[i];
   }
   free(tmp);
}

////////////////////////////////////////////////////////////////////////////////
//! malloc_2d
//! This is a kludge to dynamically allocate a contiguous 2d array
//! that will still be able to be accessed with arr[i][j]
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void** malloc_2d (size_t nrows, size_t ncols, int elementsize) {
   size_t i;
   void ** ptr;
   if ( (ptr = (void**)malloc(nrows * sizeof(void *))) == NULL ) {
      fprintf(stderr, "malloc_2d: out of memory\n");
      exit(1);
   }
   if ( (ptr[0] = malloc(nrows * ncols * elementsize)) == NULL ) {
      fprintf(stderr, "malloc_2d: out of memory\n");
      exit(1);
   }
   /* this is dangerous. I am assuming that char* has the same size as whatever
      pointer I am using at all times */
   for (i=1; i<nrows; i++) 
      ptr[i] = (char*)ptr[0] + i * ncols * elementsize;
   return ptr;
}

////////////////////////////////////////////////////////////////////////////////
//! read_gridspec
//! This function reads land/sea mask, ice and water fractions, grid
//! cell area, and country codes from the grid spec file 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
size_t read_gridspec (UserData* data) {
   int rv, ncid, varid;
   unsigned int i, j;
   if ((rv = nc_open(data->gridspec, NC_NOWRITE, &ncid)))
      NCERR(data->gridspec, rv);

   size_t index[2], count[2]; 
   get_grid_bounds(ncid, data);

   index[0] = data->start_lat;
   index[1] = data->start_lon;

   count[0] = data->n_lat;
   count[1] = data->n_lon;

   // This is the only layer that is for whole region when using mpi
   // All others are local to processor
   data->wtr_ice_f = (double **)malloc_2d(data->n_lat, data->n_lon, sizeof(double));
   printf("read_water_and_ice_fractions...\n");
   if ((rv = nc_inq_varid(ncid, "wtr_ice_frac", &varid)))
      NCERR("wtr_ice_frac", rv);
   if ((rv = nc_get_vara_double(ncid, varid, index, count, &data->wtr_ice_f[0][0])))
      NCERR("wtr_ice_frac", rv);

   size_t offset = 0;
#ifdef USEMPI
   size_t old_start_lon = data->start_lon;
   narrow_lonbounds_to_proc(*data);
   offset = data->start_lon - old_start_lon;

   index[1] = data->start_lon;
   count[1] = data->n_lon;
#endif

   /* allocate arrays */
   data->map = (site ***)malloc_2d(data->n_lat, data->n_lon, sizeof(site*));
   data->grid_cell_area = (double **)malloc_2d(data->n_lat, data->n_lon, sizeof(double));
   data->grid_cell_area_total = (double **)malloc_2d(data->n_lat, data->n_lon, sizeof(double));

#if 0
   data->mask = (unsigned char **)malloc_2d(data->n_lat, data->n_lon, sizeof(unsigned char));
   data->gcode = (int **)malloc_2d(data->n_lat, data->n_lon, sizeof(int));

   printf("read_mask...\n");
   if ((rv = nc_inq_varid(ncid, "mask", &varid)))
      NCERR("mask", rv);
   if ((rv = nc_get_vara_uchar(ncid, varid, index, count, &(data->mask[0][0]))))
      NCERR("mask", rv);
#endif


   printf("read_grid_cell_area...\n");
   if ((rv = nc_inq_varid(ncid, "grid_cell_area", &varid)))
      NCERR("gridcellarea", rv);
   if ((rv = nc_get_vara_double(ncid, varid, index, count, &data->grid_cell_area_total[0][0])))
      NCERR("gridcellarea", rv);
   size_t nPotentialSites = 0;
   for (i=0; i<data->n_lat; i++) {
      for (j=0; j<data->n_lon; j++) {
         data->grid_cell_area_total[i][j] *= 1000000.0; /*km2 to m2*/
         data->grid_cell_area[i][j] = data->grid_cell_area_total[i][j] 
            * (1.0000 - data->wtr_ice_f[i][j+offset]);
         if (data->grid_cell_area[i][j] < 0) {
            data->grid_cell_area[i][j] = 0.0;
         } else {
            nPotentialSites++;
         }
      }
   }

#if 0
   printf("read_country_codes...\n");
   if ((rv = nc_inq_varid(ncid, "countrycode", &varid))) {
      NCERR("countrycode", rv);
   }
   if ((rv = nc_get_vara_int(ncid, varid, index, count, &data->gcode[0][0]))) {
      NCERR("countrycode", rv);
   }
#endif

   if ((rv = nc_close(ncid))) {
      NCERR(data->gridspec, rv);  
   }
   return nPotentialSites;
}


////////////////////////////////////////////////////////////////////////////////
//! SiteData
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
SiteData::SiteData (size_t y, size_t x, UserData& data) {
   y_ = y;
   x_ = x;

   globY_ = y + data.start_lat;
   globX_ = x + data.start_lon;

   lat_ = data.lats[y];
   lon_ = data.lons[x];

   sprintf(name_, "lat%.2flon%.2f", lat_, lon_);

   soi = is_soi(lat_, lon_, data);

   // site area
   grid_cell_area = data.grid_cell_area[y][x];
   grid_cell_area_total = data.grid_cell_area_total[y][x];

   // TODO: should the rest be here? Same for all sites

   // max soil evaporation per mm of soil moisture per yr
   soil_evap_conductivity = 0.0;

   loss_fraction[0] = 0.0; // nothing lost during treefall
   loss_fraction[1] = data.smoke_fraction;

#ifdef ED
   if (data.open_cycles) {
      /* assume a basin wide avg of 1000 mm/yr and a basin wide avg *
       * deposisiton of N of 0.001 kg N /m2 / y, the latter is      *
       * reported in a ref in Vitousek et al 1986                   */
      N_conc_in_rain = 0.0000001; // kg N / m2 * mm
   }
   L_top  = 1.0;  // Light at the top of the canopy
   Rn_top = 1.0;  // Net Radn flx at the top of the canopy
#endif
}


////////////////////////////////////////////////////////////////////////////////
//! readSiteData
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
bool SiteData::readSiteData (UserData& data) {

   if ( ! readEnvironmentalData(data) ) {
      return false;
   }

   calcPETAverage ();

#ifdef ED
   for (size_t i=0; i<12; i++) {
      pet[i] = calcPETMonthly (i);
   }
   calcSiteDrynessIndex(data);
   
   // read in physiology or FTS data
#if FTS
   if (! readFTSdata(data) ) {
      return false;
   }
#else
   if (! readMechanismLUT(data) ) {
         return false;
   }
#endif
    

#elif defined MIAMI_LU
   // calculate miami npp
   miami_npp = miami(precip_average, temp_average);
#endif

   return true;
}


#ifdef MIAMI_LU
////////////////////////////////////////////////////////////////////////////////
//! readEnvironmentalData
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
bool SiteData::readEnvironmentalData (UserData& data) {
   int rv, ncid, varid;

   size_t index[2];
   index[0] = globX_;
   index[1] = globY_;

   if (data.climate_file_ncid == 0) {
      if ((rv = nc_open(data.climate_file, NC_NOWRITE, &ncid))) {
         NCERR(data.climate_file, rv);
      }
      data.climate_file_ncid = ncid;
   } else {
      ncid = data.climate_file_ncid;
   }

   // precip 
   if ((rv = nc_inq_varid(ncid, "annual_precipitation", &varid))) {
      NCERR("precip", rv);
   }
   if ((rv = nc_get_var1_double(ncid, varid, index, &precip_average))) {
      NCERR("precip", rv);
   }
   if (precip_average < 0) {
      //fprintf(stderr, "No precip data found for site lat %f lon %f\n", cs->lat, cs->lon);
      return false;
   } 

   if ((rv = nc_inq_varid(ncid, "annual_temperature", &varid))) {
      NCERR("temp", rv);
   }
   if ((rv = nc_get_var1_double(ncid, varid, index, &temp_average))) {
      NCERR("temp", rv);
   }
   if (temp_average == -9999.0) { // TODO: better test?
      //fprintf(stderr, "No temp data found for site lat %f lon %f\n", cs->lat, cs->lon);
      return false;
   } else {
      // TODO: check for need to convert?
      temp_average -= 273.2; // convert to celsius
   }

   if ((rv = nc_inq_varid(ncid, "annual_soil_temp", &varid))) {
      // No soil_temp? default to use air temp
      if ((rv = nc_inq_varid(ncid, "annual_temperature", &varid))) {
         NCERR("soil_temp", rv);
      }
   }
   if ((rv = nc_get_var1_double(ncid, varid, index, &soil_temp_average)))
      NCERR("soil_temp", rv);
   if (soil_temp_average == -9999.0) { // TODO: better test?
      //fprintf(stderr, "No soil temp data found for site lat %f lon %f\n", cs->lat, cs->lon);
      return false;
   } else {
      // TODO: check for need to convert?
      soil_temp_average -= 273.2; // convert to celsius
   }

   return true;
}
#endif /* MIAMI_LU */


#ifdef ED
////////////////////////////////////////////////////////////////////////////////
//! readEnvironmentalData
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
bool SiteData::readEnvironmentalData (UserData& data) {
   int rv, ncid, varid;

   size_t index1[2] = { globY_, globX_ };
   size_t index2[3] = { 0, globY_, globX_ };
   size_t count[3]  = { N_CLIMATE, 1, 1 };

   if (data.soil_file_ncid == 0) {
      if ((rv = nc_open(data.soil_file, NC_NOWRITE, &ncid))) {
         NCERR(data.soil_file, rv);
      }
      data.soil_file_ncid = ncid;
   } else {
      ncid = data.soil_file_ncid;
   }
  
   // soil chars 
   if ((rv = nc_inq_varid(ncid, "soil_depth", &varid))) {
      NCERR("soil_depth", rv);
   }
   if ((rv = nc_get_var1_double(ncid, varid, index1, &soil_depth))) {
      NCERR("soil_depth", rv);
   }
   //soil_depth *= 10.0; // convert from cm to mm

   if ((rv = nc_inq_varid(ncid, "soil_theta_max", &varid))) {
      NCERR("theta_max", rv);
   }
   if ((rv = nc_get_var1_double(ncid, varid, index1, &theta_max))) {
      NCERR("theta_max", rv);
   }

   if ((rv = nc_inq_varid(ncid, "soil_k_sat", &varid))) {
      NCERR("k_sat", rv);
   }
   if ((rv = nc_get_var1_double(ncid, varid, index1, &k_sat))) {
      NCERR("k_sat", rv);
   }

   if ((rv = nc_inq_varid(ncid, "soil_tau", &varid))) {
      NCERR("tau", rv);
   }
   if ((rv = nc_get_var1_double(ncid, varid, index1, &tau))) {
      NCERR("tau", rv);
   }
    
   if ( (soil_depth <= 0) || (theta_max == -9999.0) 
        || (k_sat == -9999.0) || (tau == -9999.0) ) {
      //fprintf(stderr, "No soil char data found for site lat %f lon %f\n", cs->lat, cs->lon);
      return false;
   } 

   // Added below for yearly climate 
   // TODO: this should not be done here. needs to be done once, not every site
   char nc[4] = ".nc"  ;
   char base[256] = "";  
   char convert[256]; 
   char climatename[256];   
   if (data.do_yearly_mech) {
      if (data.m_int) {
         sprintf(convert, "%s%d", base, data.mechanism_year);  
         strcpy(climatename, data.climate_file);
         strcat(climatename, convert);
         strcat(climatename, nc);
      }
    
      if(data.m_string) {
         strcpy(climatename,data.climate_file);
         strcat(climatename,data.mech_year_string);
         strcat(climatename, nc);
      }
   } else if(data.single_year) {
      strcpy(climatename, data.climate_file);
   }

   if (data.climate_file_ncid == 0) {
      if ((rv = nc_open(climatename, NC_NOWRITE, &ncid))) {
         NCERR(climatename, rv);
      }
      data.climate_file_ncid = ncid;
   } else {
      ncid = data.climate_file_ncid;
   }
   
   double climate_temp[12], climate_precip[12], climate_soil[12];

   // precip 
   if ((rv = nc_inq_varid(ncid, "precipitation", &varid))) {
      NCERR("precip", rv);
   }
   if ((rv = nc_get_vara_double(ncid, varid, index2, count, &climate_precip[0]))) {
      NCERR("precip", rv);
   }
   if (climate_precip[0] < 0.0) {
      //fprintf(stderr, "No precip data found for site lat %f lon %f\n", cs->lat, cs->lon);
      return false;
   } 
   
   if ((rv = nc_inq_varid(ncid, "temperature", &varid))) {
      NCERR("temp", rv);
   }
   if ((rv = nc_get_vara_double(ncid, varid, index2, count, &climate_temp[0]))) {
      NCERR("temp", rv);
   }
   if (climate_temp[0] == -9999.0) { // TODO: better test?
      //fprintf(stderr, "No temp data found for site lat %f lon %f\n", cs->lat, cs->lon);
      return false;
   } 
   
   if ((rv = nc_inq_varid(ncid, "soil_temp", &varid))) {
      // if no soil_temp, default to air temp
      if ((rv = nc_inq_varid(ncid, "temperature", &varid))) {
         NCERR("soil_temp", rv);
      }
   }
   if ((rv = nc_get_vara_double(ncid, varid, index2, count, &climate_soil[0]))) {
      NCERR("soil_temp", rv);
   }
   if (climate_soil[0] == -9999.0) { // TODO: better test?
      //fprintf(stderr, "No soil temp data found for site lat %f lon %f\n", cs->lat, cs->lon);
      return false;
   } 
   
   precip_average    = 0.0;
   temp_average      = 0.0;
   soil_temp_average = 0.0;
   for (size_t i=0; i<N_CLIMATE; i++) {
      precip[i] = climate_precip[i*12/N_CLIMATE]*12.; // to make units avg mm/yr 
      precip_average += precip[i] / N_CLIMATE;
      temp[i] = climate_temp[i*12/N_CLIMATE] - 273.15; 
      temp_average += temp[i] / N_CLIMATE;
      soil_temp[i] = climate_soil[i*12/N_CLIMATE] - 273.15; // convert Kelvin to C 
      soil_temp_average += soil_temp[i] / N_CLIMATE;
   }
   return true;
}

#if FTS
////////////////////////////////////////////////////////////////////////////////
//! readFTSData
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////

bool SiteData::readFTSdata (UserData& data) {
   
   size_t lengthp;
   int t_id, rv, ncid, varid;
   size_t index3[3] = { 0, globY_, globX_ };
   size_t count3[3]  = {N_CLIMATE_INPUT*CLIMATE_INPUT_INTERVALS, 1, 1 };
   char qairname[256], tairname[256], swname[256];
   
   // Read in the appropriate TEMP/HUM/PAR files 
   strcpy(qairname, data.QAIR_FILE);
   if ((rv = nc_open(qairname, NC_NOWRITE, &ncid))){ 
      NCERR(qairname, rv); 
   }
   else{
      data.qair_file_ncid = ncid;
   }
   if ((rv=nc_inq_dimid(ncid, "time", &t_id))) NCERR("time", rv);
   if ((rv=nc_inq_dimlen(ncid, t_id, &lengthp))) NCERR("time_len", rv);
   if (lengthp!=N_CLIMATE_INPUT*CLIMATE_INPUT_INTERVALS) {printf("Input file time length not correct - qair\n"); exit(1);}
   if ((rv = nc_inq_varid(ncid, "qair", &varid))) NCERR("qair", rv);
   if ((rv = nc_get_vara_double(ncid, varid, index3, count3, &Input_Specific_Humidity[0]))) NCERR("qair", rv);
   rv = nc_close(ncid);
   
   strcpy(tairname, data.TAIR_FILE);
   if ((rv = nc_open(tairname, NC_NOWRITE, &ncid))){ 
      NCERR(tairname, rv); 
   }
   else{
      data.tair_file_ncid = ncid;
   }
   if ((rv=nc_inq_dimid(ncid, "time", &t_id))) NCERR("time", rv);
   if ((rv=nc_inq_dimlen(ncid, t_id, &lengthp))) NCERR("time_len", rv);
   if (lengthp!=N_CLIMATE_INPUT*CLIMATE_INPUT_INTERVALS) {printf("Input file time length not correct - qair\n"); exit(1);}
   if ((rv = nc_inq_varid(ncid, "tair", &varid))) NCERR("tair", rv);
   if ((rv = nc_get_vara_double(ncid, varid, index3, count3, &Input_Temperature[0]))) NCERR("tair", rv);
   rv = nc_close(ncid);
   
   strcpy(swname, data.SW_FILE);
   if ((rv = nc_open(swname, NC_NOWRITE, &ncid))){ 
      NCERR(swname, rv); 
   }
   else{
      data.sw_file_ncid = ncid;
   }
   if ((rv=nc_inq_dimid(ncid, "time", &t_id))) NCERR("time", rv);
   if ((rv=nc_inq_dimlen(ncid, t_id, &lengthp))) NCERR("time_len", rv);
   if (lengthp!=N_CLIMATE_INPUT*CLIMATE_INPUT_INTERVALS) {printf("Input file time length not correct - qair\n"); exit(1);}
   if ((rv = nc_inq_varid(ncid, "swdown", &varid))) NCERR("swdown", rv);
   if ((rv = nc_get_vara_double(ncid, varid, index3, count3, &Input_Par[0]))) NCERR("swdown", rv);
   rv = nc_close(ncid);
   
   //Convert inputs to right units
   int i;
   for (i=0;i<N_CLIMATE_INPUT*CLIMATE_INPUT_INTERVALS;i++){
      Input_Temperature[i]-=273.2;
      Input_Par[i]/=2.;
      Input_Specific_Humidity[i]*=28.96 / 18.02;
   }
   
   // Monthly Averaging done the same way as mech processing 
   if (1){ 
      int ndays = 365;
      float dpp = ndays/12.0;
      float dc = 0.;
      //Two grids to keep the leftover fraction of each month
      float averages[12] ={0};
      float averages2[12] = {0};
      int start_day = 0;
      int j,k,s,t;
      t = 0;
      float check = 0;
      //loop through days and apply appropriate fraction for that day/hour
      for (i=0;i<ndays;i++){
         dc +=1 ;
         if (dc >= dpp){
            float frac = 1 - (dc-dpp);
            for (k=0;k<4;k++){
               if (t%2==0){
                  averages[k]+=Input_Specific_Humidity[4*i+k]*frac;
                  averages[4+k]+=Input_Temperature[4*i+k]*frac;
                  averages[8+k]+=Input_Par[4*i+k]*frac;
                  averages2[k]+=Input_Specific_Humidity[4*i+k]*(1-frac);
                  averages2[4+k]+=Input_Temperature[4*i+k]*(1-frac);
                  averages2[8+k]+=Input_Par[4*i+k]*(1-frac);
                  
               }
               else{
                  averages[k]+=Input_Specific_Humidity[4*i+k]*(1-frac);
                  averages[4+k]+=Input_Temperature[4*i+k]*(1-frac);
                  averages[8+k]+=Input_Par[4*i+k]*(1-frac);
                  averages2[k]+=Input_Specific_Humidity[4*i+k]*(frac);
                  averages2[4+k]+=Input_Temperature[4*i+k]*(frac);
                  averages2[8+k]+=Input_Par[4*i+k]*(frac);
               }
            }
            dc = (1-frac);
            for (j=start_day;j<=i;j++){
               for (k=0;k<4;k++){
                  if (t%2==0){
                     Input_Specific_Humidity[4*j+k] = averages[k]/(dpp);
                     Input_Temperature[4*j+k] = averages[4+k]/(dpp);
                     Input_Par[4*j+k] = averages[8+k]/(dpp);
                  }
                  else{
                     Input_Specific_Humidity[4*j+k] = averages2[k]/(dpp);
                     Input_Temperature[4*j+k] = averages2[4+k]/(dpp);
                     Input_Par[4*j+k] = averages2[8+k]/(dpp);
                  }
               }
            }
            for (s=0;s<12;s++){
               if (t%2==0){
                  averages[s]=0;
               }
               else{
                  averages2[s]=0;
               }
            }
            t +=1;
            start_day = i+1;
         }
         else{
            for (k=0;k<4;k++){
               if (t%2 ==0){
                  averages[k]+=Input_Specific_Humidity[4*i+k];
                  averages[4+k]+=Input_Temperature[4*i+k];
                  averages[8+k]+=Input_Par[4*i+k];
               }
               else{
                  averages2[k]+=Input_Specific_Humidity[4*i+k];
                  averages2[4+k]+=Input_Temperature[4*i+k];
                  averages2[8+k]+=Input_Par[4*i+k];
               }
            }
         }
      }
   }
   // End averaging  
   return true;
}
#else
////////////////////////////////////////////////////////////////////////////////
//! readMechanismLUT
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
bool SiteData::readMechanismLUT (UserData& data) {
   int rv, ncid, varid;
   size_t index1[3] = { globY_, globX_, 0 };
   size_t count1[3] = { 1, 1, N_CLIMATE };

   size_t index2[4] = { globY_, globX_, 0, 0 };
   size_t count2[4] = { 1, 1, N_CLIMATE, N_LIGHT };

   // TODO: this shouldn't be here. Do once, not for each site.
   char nc[] = ".nc"  ;
   char base[] = "";  
   char convert[256];
   char c3name[256]; 
   char c4name[256];  
   if (data.do_yearly_mech) {
      // NOTE!!! do_yearly_mech does not work with the multiple Vm0 bins mechanism (currently)
      // Please use FTS
       if (data.m_int) {
          sprintf(convert, "%s%d", base, data.mechanism_year);  
          strcpy(c3name, data.mech_c3_file);
          strcpy(c4name, data.mech_c4_file);
     
          strcat(c3name, convert);
          strcat(c4name, convert);   
    
          strcat(c3name, nc);
          strcat(c4name, nc); 
       } else if (data.m_string) {
          strcpy(c3name, data.mech_c3_file);
          strcpy(c4name, data.mech_c4_file);
    
          strcat(c3name, data.mech_year_string);
          strcat(c4name, data.mech_year_string);   
    
          strcat(c3name, nc);
          strcat(c4name, nc); 
       }
   } else if (data.single_year) {
      size_t i = 0;
      for(; i < data.num_Vm0; i++) {
          if (data.num_Vm0 > 1) {
            strcpy(c3name, data.Vm0_basepath);
            strcat(c3name, data.list_c3_files.at(i).c_str());
            strcpy(c4name, data.Vm0_basepath);
            strcat(c4name, data.list_c4_files.at(i).c_str());
         } else {
            strcpy(c3name, data.mech_c3_file);
            strcpy(c4name, data.mech_c4_file);
         }

         if (data.mech_c3_file_ncid[i] == 0) {
            if ((rv = nc_open(c3name, NC_NOWRITE, &ncid))) {
               NCERR(c3name, rv);
            } else {
               data.mech_c3_file_ncid[i] = ncid;
            }
         }
         if (data.mech_c4_file_ncid[i] == 0) {
            if ((rv = nc_open(c4name, NC_NOWRITE, &ncid))) {
               NCERR(c4name, rv);
            } else {
               data.mech_c4_file_ncid[i] = ncid;
            }
         }
         
         for (size_t pt=0; pt<PT; pt++) {
            if (pt == 0) {
               ncid = data.mech_c3_file_ncid[i];
            } else {
               ncid = data.mech_c4_file_ncid[i];
            }
            // light levels 
            // TODO: is this worth storing for each site? Same for everywhere.
            if ((rv = nc_inq_varid(ncid, "shade", &varid))) {
               NCERR("shade", rv);
            }
            if ((rv = nc_get_var_double(ncid, varid, &light_levels[pt][i][0]))) {
               NCERR("shade", rv);
            }

            // temp function 
            if ((rv = nc_inq_varid(ncid, "tf", &varid))) {
               NCERR("tf", rv);
            }
            if ((rv = nc_get_vara_double(ncid, varid, index1, count1, &tf[pt][i][0]))) {
               NCERR("tf", rv);
            }

            // An 
            if ((rv = nc_inq_varid(ncid, "An", &varid))) {
               NCERR("An", rv);
            }
            if ((rv = nc_get_vara_double(ncid, varid, index2, count2, &An[pt][i][0][0]))) {
               NCERR("An", rv);
            }

            // Anb 
            if ((rv = nc_inq_varid(ncid, "Anb", &varid))) {
               NCERR("Anb", rv);
            }
            if ((rv = nc_get_vara_double(ncid, varid, index2, count2, &Anb[pt][i][0][0]))) {
               NCERR("Anb", rv);
            }

            // E 
            if ((rv = nc_inq_varid(ncid, "E", &varid))) {
               NCERR("E", rv);
            }
            if ((rv = nc_get_vara_double(ncid, varid, index2, count2, &E[pt][i][0][0]))) {
               NCERR("E", rv);
            }

            // Eb 
            if ((rv = nc_inq_varid(ncid, "Eb", &varid))) {
               NCERR("Eb", rv);
            }
            if ((rv = nc_get_vara_double(ncid, varid, index2, count2, &Eb[pt][i][0][0]))) {
               NCERR("Eb", rv);
            }
         }
      }
   }
   // TODO: check for bad inputs
   
   return true;
}
#endif

////////////////////////////////////////////////////////////////////////////////
//! calcPETMonthly
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double SiteData::calcPETMonthly (size_t month) {
   /* Based on Thornthwaite 1948, as presented in Bras, R. L. 1990. Hydrology: *
    * An introduction to hydrologic science. Addison-Wesley. NY. p.224.        */

   /* function used in soil biogeochemistry      */
   /* This function returns monthly pet in mm/yr */

   /* adjustment factor for daylength and month length- lat dependent, this   *
    * value is the average of the values given in table 5.9 on page 225 of    *
    * Bras 1990                                                               */

   double b[12];
   if ( (lat_ <= 5.0) && (lat_ >= -5.0) ) {
      b[0]=1.04; b[1]=0.94; b[2]=1.04; b[3]=1.01; b[4]=1.04; b[5]=1.01;
      b[6]=1.04; b[7]=1.04; b[8]=1.01; b[9]=1.04; b[10]=1.01; b[11]=1.04;
   }
   else if ( (lat_ <= 15.0) && (lat_ >= -15.0) ) {
      b[0]=1.00; b[1]=0.91; b[2]=1.03; b[3]=1.03; b[4]=1.08; b[5]=1.06;
      b[6]=1.08; b[7]=1.07; b[8]=1.02; b[9]=1.02; b[10]=0.98; b[11]=0.99;
   }
   else if ( (lat_ <= 25.0) && (lat_ >= -25.0) ) {
      b[0]=0.95; b[1]=0.90; b[2]=1.03; b[3]=1.05; b[4]=1.13; b[5]=1.11;
      b[6]=1.14; b[7]=1.11; b[8]=1.02; b[9]=1.00; b[10]=0.93; b[11]=0.94;
   }
   else if ( (lat_ <= 32.5) && (lat_ >= -32.5) ) {
      b[0]=0.90; b[1]=0.87; b[2]=1.03; b[3]=1.08; b[4]=1.18; b[5]=1.17;
      b[6]=1.20; b[7]=1.14; b[8]=1.03; b[9]=0.98; b[10]=0.89; b[11]=0.88;
   }
   else if ( (lat_ <= 37.5) && (lat_ >= -37.5) ) {
      b[0]=0.87; b[1]=0.85; b[2]=1.03; b[3]=1.09; b[4]=1.21; b[5]=1.21;
      b[6]=1.23; b[7]=1.16; b[8]=1.03; b[9]=0.97; b[10]=0.86; b[11]=0.85;
   }
   else if ( (lat_ <= 42.5) && (lat_ >= -42.5) ) {
      b[0]=0.84; b[1]=0.83; b[2]=1.03; b[3]=1.11; b[4]=1.24; b[5]=1.25;
      b[6]=1.27; b[7]=1.18; b[8]=1.04; b[9]=0.96; b[10]=0.83; b[11]=0.81;
   }
   else if ( (lat_ <= 47.5) && (lat_ >= -47.5) ) {
      b[0]=0.80; b[1]=0.81; b[2]=1.02; b[3]=1.13; b[4]=1.28; b[5]=1.29;
      b[6]=1.31; b[7]=1.21; b[8]=1.04; b[9]=0.94; b[10]=0.79; b[11]=0.75;
   }
   else if( (lat_ <= 52.5) && (lat_ >= -52.5) ) {
      b[0]=0.74; b[1]=0.78; b[2]=1.02; b[3]=1.15; b[4]=1.33; b[5]=1.36;
      b[6]=1.37; b[7]=1.25; b[8]=1.06; b[9]=0.92; b[10]=0.76; b[11]=0.70;
   }
   /*WARNING! These next latitudes are off the table! values are wrong*/
   else if ( (lat_ <= 90.0) && (lat_ >= -90.0) ) {
      b[0]=0.74; b[1]=0.78; b[2]=1.02; b[3]=1.15; b[4]=1.33; b[5]=1.36;
      b[6]=1.37; b[7]=1.25; b[8]=1.06; b[9]=0.92; b[10]=0.76; b[11]=0.70;
   }

   double I = 0.0; /* annual heat index */
   for(size_t i=0; i<12; i++){
      double T = temp[i];
      /*printf("T %f \t pow(T/5.0,1.51) %f\n",T,pow(T/5.0,1.51));*/

      /* KLUDGE prevent PET calc from crashing,        *
       * cant raise negative numner to the 1.51 power! */
      if (T < 0.0001) T = 0.0001;

      I += pow(T / 5.0, 1.51);
      /*printf("T %f \t pow(T/5.0,1.51) %f\n",T,pow(T/5.0,1.51));*/
   }

   double a = 67.5 * pow(10.0, -8.0) * pow(I, 3.0) - 77.1 * pow(10.0, -6.0)
            * pow(I, 2.0) + 0.0179 * I + 0.492;

   double T = temp[month];
   /* KLUDGE prevent PET calc from crashing;        *
    * cant raise negative number to the 1.51 power! */
   if (T < 0.0001) T = 0.0001;

   /* monthly potential evapotranspiration (cm/month) */
   double pet = 1.62 * b[month] * pow(10.0 * T / I, a);
   pet *= 10.0; /*convert cm to mm*/
   pet *= 12.0; /*convert mm/month to mm/yr*/

#if 0
   printf("pet monthly: month %d b= %f\n",month,b[month]);
   for(i=0;i<12;i++) printf("pet monthly: month %d temp %f pet= %f \n",month,T, pet);
#endif

   return(pet);
}


////////////////////////////////////////////////////////////////////////////////
//! calcSiteDrynessIndex
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void SiteData::calcSiteDrynessIndex (UserData& data) {

   size_t month = 0;
   double di = 0.0;
   if (pet[month] / precip[month] > 1.0) {
      di += 30.0 * (pet[month] - precip[month]);
      int stop = 0;
      if (stop == 0) {
         for (size_t month2=11; month2>month; month2--) {
            if ( (pet[month2] / precip[month2] > 1.0) && (stop == 0) ) {
               di += 30.0 * (pet[month2] - precip[month2]);
            } else {
               stop = 1; // wet month breaks consectuive string
            }
         }
      }

      dryness_index[month] = di;
      //printf("month %d di_month %f\n",month,cs->sdata->dryness_index[month]);
   }

   for (size_t month=1; month<12; month++) {
      if (pet[month] / precip[month] > 1.0) {
         di = dryness_index[month-1] + 30.0 * (pet[month] - precip[month]);
      } else {
         di = 0.0;
      }

      dryness_index[month] = di;
      //printf("month %d di_month %f\n",month,cs->sdata->dryness_index[month]);
   }

   // average
   dryness_index_avg = 0.0;
   for (size_t month=0; month<12; month++) {
      dryness_index_avg += dryness_index[month] / 12.0;
   }

   //printf("di avg %f\n",cs->sdata->dryness_index_avg);
}
#endif /* ED */


////////////////////////////////////////////////////////////////////////////////
//! calcPETAverage
//! Based on Thornthwaite 1948, as presented in Bras, R. L. 1990. Hydrology
//! An introduction to hydrologic science. Addison-Wesley. NY. p.224.
//! function used in soil biogeochemistry, this function returns annual average pet in mm/yr
//! this function is the function of the average adjustment factor for daylength 
//! and month length- lat dependent, this value is the average of the values given in 
//! table 5.9 on page 225 of Bras 1990
//! 
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void SiteData::calcPETAverage () {
   double b = 1.02;

   double T = temp_average;
   if (T < 0.0001) T = 0.0001; // cant take negative number to the 1.51 power!
   double I = 12 * pow(T / 5.0, 1.51); // annual heat index
   double a = 67.5 * pow(10.0, -8.0) * pow(I, 3.0) - 77.1 * pow(10.0, -6.0)
   * pow(I, 2.0) + 0.0179 * I + 0.492;

   // monthly potential evapotranspiration (cm/month)
   // last two terms convert cm to mm and mm/month to mm/yr
   pet_average = 1.62 * b * pow(10.0 * T / I, a) * 10.0 * 12.0;
   //printf("hydrology: pet yrly average: pet= %f \n",pet_average);
}



////////////////////////////////////////////////////////////////////////////////
//! read_hurricane_disturbance
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
int read_hurricane_disturbance(site** current_site, UserData* data) {
   int rv, ncid, varid;
   size_t index1[2], index2[3], count[3];

   site* cs = *current_site;

   index1[0] = cs->sdata->lat_ - 25.5;
   index1[1] = cs->sdata->lon_ + 101.5;

   index2[0] = 0;
   index2[1] = cs->sdata->lat_ - 25.5;
   index2[2] = cs->sdata->lon_ + 101.5;

   count[0] = data->n_hurricane_years;
   count[1] = 1;
   count[2] = 1;

   cs->sdata->hurricane_disturbance_rate = (double *)malloc(data->n_hurricane_years * sizeof(double));
   if ((rv = nc_open(data->hurricane_file, NC_NOWRITE, &ncid)))
      NCERR(data->hurricane_file, rv);

   if ((rv = nc_inq_varid(ncid, "average", &varid)))
      NCERR("average", rv);
   if ((rv = nc_get_var1_double(ncid, varid, index1, 
                               &(cs->sdata->avg_hurricane_disturbance_rate))))
      NCERR("average", rv);

   if ((rv = nc_inq_varid(ncid, "damage", &varid)))
      NCERR("damage", rv);
   if ((rv = nc_get_vara_double(ncid, varid, index2, count, 
                               &(cs->sdata->hurricane_disturbance_rate[0]))))
      NCERR("damage", rv);

   if (cs->sdata->avg_hurricane_disturbance_rate <= 0.0) {
      //fprintf(stderr, "No hurricane data found for site lat %f lon %f\n", cs->lat, cs->lon);
      return 0;
   } 

#if 0 /*temp set avg to 0 for landuse runs from potential -- needs it for years before hurr rec*/
   cs->sdata->avg_hurricane_disturbance_rate = 0;
#endif
   if ((rv = nc_close(ncid)))
      NCERR(data->hurricane_file, rv);  

   return 1;
}

/******************************************************************************/
