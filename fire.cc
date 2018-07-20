#include <cstdio>
#include <cstring>
#include <cmath>
#include "netcdf.h"
#include "miami.h"
#include "read_site_data.h"

#include "edmodels.h"
#include "site.h"
#include "patch.h"
#include "read_site_data.h"
#ifdef ED
#include "cohort.h"
#endif

///////////////////////////////////////////////////////////////////////////////
//! read_gfed_bf
//! Yannick GFED: function to read burned fraction from GFED
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
/* Yannick GFED: function to read burned fraction from GFED */
//void read_gfed_bf (UserData* data) {
//   int rv, ncid, varid;
//   size_t i,j;
//   size_t index[2], count[2];
//
//   printf("read_gfed_bf...\n");
//
//   index[0] = data->start_lat;
//   index[1] = data->start_lon;
//
//   count[0] = data->n_lat;
//   count[1] = data->n_lon;
//
//   data->gfed_bf = (double **)malloc_2d(data->n_lat, data->n_lon, sizeof(double));
//   printf("%s\n",data->gfedbf_file);
//
//   if ((rv = nc_open(data->gfedbf_file, NC_NOWRITE, &ncid))) {
//      NCERR(data->gfedbf_file, rv);
//   }
//
//   if ((rv = nc_inq_varid(ncid, "burnedfraction", &varid))) {
//      NCERR("burnedfraction", rv);
//   }
//
//   if ((rv = nc_get_vara_double(ncid, varid, index, count, &(data->gfed_bf[0][0])))) {
//      NCERR("burnedfraction", rv);
//   }
//
//   /* set all missing values to 0 */
//   for (i=0; i<data->n_lat; i++) {
//      for (j=0; j<data->n_lat; j++) {
//         if (data->gfed_bf[i][j] < 0.0) {
//            data->gfed_bf[i][j] = 0.0;
//         }
//      }
//   }
//
//}


////////////////////////////////////////////////////////////////////////////////
//! fire
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double fire (int t, patch** patchptr, UserData* data) {
   /* for all linked patches */

   if(data->fire_off) {
      return 0.0;
   }

   patch* cp =* patchptr;
   /* no fires on cropland */
   if (cp->landuse == LU_CROP) {
      return 0.0; 
   }

   double fireterm = 0.0;
   /* Yannick: returns gfed burned fraction if option activated (gfed_bf) */
    // Lei: add time dimension to gfed_bf for loading monthly GFED data
   if(data->gfed_bf) {
       if (PATCH_FREQ==12)
       {
           for (size_t mon=0;mon<12.0;mon++)
           {
               fireterm += data->gfed_bf[data->time_period][cp->siteptr->sdata->globY_][cp->siteptr->sdata->globX_];
           }
           fireterm /=12.0;   // Lei- If raw GFED burned fraction is multiplied by 12 in load_GFED(), here to take a average than sum. If no, here take sum to get yearly total disturbance rate
       }
       else
           fireterm = data->gfed_bf[data->time_period][cp->siteptr->sdata->globY_][cp->siteptr->sdata->globX_];  // Lei- the unit of data->gfed_bf after loading from load_GFED should be yearly based.
      
      return(fireterm);
   }

#if defined ED 
   for ( size_t i=0; i<N_SUB; i++ ) {
      fireterm += cp->lambda1[i] / ( 1.0 * N_SUB );
   }

#elif MIAMI_LU
   site* cs = cp->siteptr;

   cp->fuel = 0.0; 
   cp->ignition_rate = 0.0;
   cp->fuel = cp->total_biomass;

   if (cp->fuel > 0.0) {
      if (cs->sdata->precip_average < 400 + 40 * cs->sdata->temp_average) {
         fireterm = cp->fuel * (400 + 40 * cs->sdata->temp_average - cs->sdata->precip_average);
      }
   }
#endif

   if (fireterm > data->fire_max_disturbance_rate) {
      fireterm = data->fire_max_disturbance_rate; 
   }

#if DOES_COMPILE
   if(data->fire_suppression) {
       double old_fire;
       /* ramp fire down 1870 to 1970 to 2% of original value; ref Houghton */

       /* using values from ed. miami-lu was using 220*N_SUB in place of 170*N_SUB */

       if( (data->fire_suppression_stop && ((t > 170*N_SUB) && (t < 300*N_SUB)) ) 
               || (!data->fire_suppression_stop && (t > 170*N_SUB))) {
  
          old_fire = fireterm;

          /* using values from ed. miami-lu was using (t-220*N_SUB)*1.0/(50*N_SUB) */
          if ( !strcmp(data->region,"NORTH_AMERICA") || !strcmp(data->region,"US")  ) {
            fireterm -=  fireterm*(t-170*NSUB)*1.0/(100.0*N_SUB);
            if(fireterm < old_fire * 0.02) 
               fireterm = old_fire * 0.02;
          }
          else if (!strcmp(data->region,"SOUTH_AMERICA")) {
            fireterm -=  fireterm*(t-170*NSUB)*1.0/(100.0*N_SUB);
            if(fireterm < old_fire * 0.09) 
               fireterm = old_fire * 0.09;
          }
          else if (!strcmp(data->region,"EUROPE")) {
            fireterm -=  fireterm*(t-170*NSUB)*1.0/(100.0*N_SUB);
            if(fireterm < old_fire * 0.04) 
               fireterm = old_fire * 0.04;
          }
          else if (!strcmp(data->region,"AFRICA")) {
            fireterm -=  fireterm*(t-170*NSUB)*1.0/(100.0*N_SUB);
            if(fireterm < old_fire * 0.43) 
               fireterm = old_fire * 0.43;
          }
          else if (!strcmp(data->region,"NASIA")) {
            fireterm -=  fireterm*(t-170*NSUB)*1.0/(100.0*N_SUB);
            if(fireterm < old_fire * 0.22) 
               fireterm = old_fire * 0.22;
          }
          else if (!strcmp(data->region,"SASIA")) {
            fireterm -=  fireterm*(t-170*NSUB)*1.0/(100.0*N_SUB);
            if(fireterm < old_fire * 0.22) 
               fireterm = old_fire * 0.22;
          }
          else if (!strcmp(data->region,"AUSTRALIA")) {
            fireterm -=  fireterm*(t-170*NSUB)*1.0/(100.0*N_SUB);
            if(fireterm < old_fire * 0.43) 
               fireterm = old_fire * 0.43;
          }
       }
   }/* FIRE_SUPPRESSION */
#endif
   return(fireterm);  
}


#ifdef ED
////////////////////////////////////////////////////////////////////////////////
//! update_fuel
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void update_fuel (int t, patch** patchptr, UserData* data) {
   /* for all linked patches */
  
   patch* youngest_patch = *patchptr;
   site* cs = youngest_patch->siteptr;

   size_t landuse = youngest_patch->landuse;

   /* init */
   cs->fuel[landuse] = 0.0;
   cs->ignition_rate[landuse] = 0.0;
  
   patch* cp = youngest_patch;
   while (cp != NULL) {
      /* init patch fuel and ignition rate */
      cp->fuel = 0.0;
      cp->ignition_rate = 0.0;

      cohort* cc = cp->shortest;
      while(cc != NULL) {
         if (cc->hite < data->fire_hite_threshold) {
            cp->fuel += cc->nindivs * cc->babove;
         }
         cc = cc->taller;
      } /* end loop over cohorts */
    
      /* divide fuel by patch area */
      cp->fuel /= cp->area;
    
      /* add patch total to site total */
      cs->fuel[landuse] += cp->fuel * cp->area;
    
      cp = cp->older;
   }  /* end loop over patches */

   /* divide by site area */
   cs->fuel[landuse] /= (cs->area_fraction[landuse] * data->area);
   cs->ignition_rate[landuse] = cs->fuel[landuse] * data->fp1 
      * pow(cs->sdata->dryness_index_avg / 30000.0, 10.0);

   /* calculate fire dist rate and add into array of within yr values */
   cs->lambda1[t%N_SUB][landuse] = 0.0;
   cp = youngest_patch;
   while (cp != NULL) {
      cp->lambda1[t%N_SUB] = cs->ignition_rate[landuse];
      cs->lambda1[t%N_SUB][landuse] += cp->lambda1[t%N_SUB] * cp->area;
      cp=cp->older;
   }  /* end loop over patches */
   cs->lambda1[t%N_SUB][landuse] /= cs->area_fraction[landuse] * data->area;  
}
#endif /* ED */

