#include <cmath>
#include <cstdio>
#include <stdio.h>
#include <string.h>
#include <cstdlib>
#include "netcdf.h"
#include "edmodels.h"
#include "site.h"
#include "patch.h"
#include "read_site_data.h"
#ifdef ED
#include "cohort.h"
#endif
#include "disturbance.h"

#if LANDUSE
void landuse_transition (site** siteptr, UserData* data, 
                         int donor_lu, int target_lu, double beta);
void cut_forest (site** siteptr, UserData* data, 
                 int donor_lu, double beta, int bh_type);
void harvest_croplands (site** siteptr, UserData* data);
void graze_pastures (site** siteptr, UserData* data);

////////////////////////////////////////////////////////////////////////////////
//! read_initial_landuse_c
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void read_initial_landuse_c (UserData* data) {
   int rv, ncid, varid;
   size_t i,j;
   size_t index[2], count[2];

   printf("read_initial_landuse_c...\n");

   index[0] = data->start_lat;
   index[1] = data->start_lon;

   count[0] = data->n_lat;
   count[1] = data->n_lon;

   data->init_csc = (double **)malloc_2d(data->n_lat, data->n_lon, sizeof(double));
   data->init_psc = (double **)malloc_2d(data->n_lat, data->n_lon, sizeof(double));
   data->init_pb = (double **)malloc_2d(data->n_lat, data->n_lon, sizeof(double));

   if ((rv = nc_open(data->lu_init_c_file, NC_NOWRITE, &ncid)))
      NCERR(data->lu_init_c_file, rv);

   if ((rv = nc_inq_varid(ncid, "crop_sc", &varid)))
      NCERR("crop_sc", rv);
   if ((rv = nc_get_vara_double(ncid, varid, index, count, &(data->init_csc[0][0]))))
      NCERR("crop_sc", rv);

   if ((rv = nc_inq_varid(ncid, "past_sc", &varid)))
      NCERR("past_sc", rv);
   if ((rv = nc_get_vara_double(ncid, varid, index, count, &(data->init_psc[0][0]))))
      NCERR("past_sc", rv);

   if ((rv = nc_inq_varid(ncid, "past_b", &varid)))
      NCERR("past_b", rv);
   if ((rv = nc_get_vara_double(ncid, varid, index, count, &(data->init_pb[0][0]))))
      NCERR("past_b", rv);
  
   /* set all missing values to 0 */
   for (i=0; i<data->n_lat; i++)
      for (j=0; j<data->n_lon; j++) {
         if (data->init_csc[i][j] < 0.0) data->init_csc[i][j] = 0.0;
         if (data->init_psc[i][j] < 0.0) data->init_psc[i][j] = 0.0;
         if (data->init_pb[i][j]  < 0.0) data->init_pb[i][j]  = 0.0;
      }
}


////////////////////////////////////////////////////////////////////////////////
//! read_initial_landuse_fractions
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void read_initial_landuse_fractions (UserData* data) {
   int rv, ncid, varid;
   size_t i, j;
   size_t start[3], count[3];

   printf("read_initial_landuse_fractions...\n");


   data->init_c = (double **)malloc_2d(data->n_lat, data->n_lon, sizeof(double));
   data->init_p = (double **)malloc_2d(data->n_lat, data->n_lon, sizeof(double));
   data->init_v = (double **)malloc_2d(data->n_lat, data->n_lon, sizeof(double));

   if ((rv = nc_open(data->lu_file, NC_NOWRITE, &ncid)))
      NCERR(data->lu_file, rv);

   start[0] = 0;
   start[1] = data->start_lat;
   start[2] = data->start_lon;

   count[0] = 1;
   count[1] = data->n_lat;
   count[2] = data->n_lon;

   if ((rv = nc_inq_varid(ncid, "gcrop", &varid)))
      NCERR("gcrop", rv);
   if ((rv = nc_get_vara_double(ncid, varid, start, count, &data->init_c[0][0])))
      NCERR("gcrop", rv);

   if ((rv = nc_inq_varid(ncid, "gpast", &varid)))
      NCERR("gpast", rv);
   if ((rv = nc_get_vara_double(ncid, varid, start, count, &data->init_p[0][0])))
      NCERR("gpast", rv);

   if ((rv = nc_inq_varid(ncid, "gothr", &varid)))
      NCERR("gothr", rv);
   if ((rv = nc_get_vara_double(ncid, varid, start, count, &data->init_v[0][0])))
      NCERR("gothr", rv);

   if ((rv = nc_close(ncid)))
      NCERR(data->lu_file, rv);


   /*re-normalize from fraction of grid-cell area to fraction of the land area*/  
   for (i=0; i<data->n_lat; i++) {
      for (j=0; j<data->n_lon; j++) {
         data->init_c[i][j] /= (data->init_c[i][j] + data->init_p[i][j] + data->init_v[i][j]);
         data->init_p[i][j] /= (data->init_c[i][j] + data->init_p[i][j] + data->init_v[i][j]);
         data->init_v[i][j] /= (data->init_c[i][j] + data->init_p[i][j] + data->init_v[i][j]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
//! lu2charname
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
char lu2charname (int lu) {
   char name;

   switch (lu) {
      case LU_NTRL: 
         name = 'v';
         break;
      case LU_SCND:
         name = 's';
         break;
      case LU_CROP:
         name = 'c';
         break;
      case LU_PAST:
         name = 'p';
         break;
      default:
         fprintf(stderr, "unkown landuse type: %d\n", lu);
         return 0;
   }
   return name;
}


////////////////////////////////////////////////////////////////////////////////
//! read_transition_rates
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
int read_transition_rates (site** current_site, UserData* data) {
   int rv, ncid, varid, dlu, tlu, i;
   char dn, tn;
   char varname[10];
   size_t index[3], count[3];
   double factor;

   site* cs = *current_site;

   index[0] = 0;
   index[1] = cs->sdata->globY_;
   index[2] = cs->sdata->globX_;

   count[0] = N_LANDUSE_YEARS;
   count[1] = 1;
   count[2] = 1;

   if (data->lu_file_ncid == 0) {
      if ((rv = nc_open(data->lu_file, NC_NOWRITE, &ncid)))
         NCERR(data->lu_file, rv);
      data->lu_file_ncid = ncid;
   } else {
      ncid = data->lu_file_ncid;
   }

   for (dlu=0; dlu<N_LANDUSE_TYPES; dlu++) {
      for (tlu=1; tlu<N_LANDUSE_TYPES; tlu++) {
         /* skip transitions self->self and v->s (v->s dealt with in sbh/vbh) */
         if ((dlu != tlu) && !( (dlu == LU_NTRL) && (tlu == LU_SCND) ) ) { 
            if ((dn = lu2charname(dlu)) && (tn = lu2charname(tlu))) {
               sprintf(varname, "gfl%c%c", dn, tn);
               if ((rv = nc_inq_varid(ncid, varname, &varid)))
                  NCERR(varname, rv);
               if ((rv = nc_get_vara_double(ncid, varid, index, count, 
                                            &(cs->sdata->beta[dlu][tlu-1][0]))))
                  NCERR(varname, rv);
            }
         }
      }
   }

   for (dlu=0; dlu<N_VBH_TYPES; dlu++) {
      sprintf(varname, "gfvh%d", dlu+1);
      if ((rv = nc_inq_varid(ncid, varname, &varid)))
         NCERR(varname, rv);

      if ((rv = nc_get_vara_double(ncid, varid, index, count, 
                                   &(cs->sdata->vbh[dlu][0]))))
         NCERR(varname, rv);
   }

   for (dlu=0; dlu<N_SBH_TYPES; dlu++) {
      sprintf(varname, "gfsh%d", dlu+1);
      if ((rv = nc_inq_varid(ncid, varname, &varid)))
         NCERR(varname, rv);

      if ((rv = nc_get_vara_double(ncid, varid, index, count, 
                                   &(cs->sdata->sbh[dlu][0]))))
         NCERR(varname, rv);
   }

   /*re-normalize from fraction of grid cell area to fraction of the land area*/
   if (cs->sdata->grid_cell_area_total != cs->sdata->grid_cell_area) {
      factor = cs->sdata->grid_cell_area_total / cs->sdata->grid_cell_area;
      for (dlu=0; dlu<N_LANDUSE_TYPES; dlu++)
         for (tlu=0; tlu<N_LANDUSE_TYPES-1; tlu++) 
            for (i=0; i<N_LANDUSE_YEARS; i++) 
               cs->sdata->beta[dlu][tlu][i] *= factor;
      for (dlu=0; dlu<N_VBH_TYPES; dlu++)
         for (i=0; i<N_LANDUSE_YEARS; i++) 
            cs->sdata->vbh[dlu][i] *= factor;
      for (dlu=0; dlu<N_SBH_TYPES; dlu++)
         for (i=0; i<N_LANDUSE_YEARS; i++) 
            cs->sdata->sbh[dlu][i] *= factor;      
   } 

#if 0
   /* TODO: this needs to be adjusted to use the multi-D lutype indexed arrays - justin */
   /* average last ten years transistion rates for bau */
   cs->sdata->beta_vs_bau = cs->sdata->beta_vc_bau = 
      cs->sdata->beta_vp_bau = cs->sdata->beta_sp_bau = 
      cs->sdata->beta_sc_bau = cs->sdata->beta_cs_bau = 
      cs->sdata->beta_ps_bau = cs->sdata->beta_gc_bau = 
      cs->sdata->beta_gp_bau = cs->sdata->beta_pg_bau = 
      cs->sdata->beta_cg_bau = cs->sdata->beta_cp_bau = 
      cs->sdata->beta_pc_bau = 0.00;

   for(j=N_LANDUSE_YEARS-11; j<N_LANDUSE_YEARS; j++){
      cs->sdata->beta_vs_bau += cs->sdata->beta_vs[j] / 10.0; 
      cs->sdata->beta_vc_bau += cs->sdata->beta_vc[j] / 10.0; 
      cs->sdata->beta_vp_bau += cs->sdata->beta_vp[j] / 10.0; 
      cs->sdata->beta_sp_bau += cs->sdata->beta_sp[j] / 10.0; 
      cs->sdata->beta_sc_bau += cs->sdata->beta_sc[j] / 10.0; 
      cs->sdata->beta_cs_bau += cs->sdata->beta_cs[j] / 10.0; 
      cs->sdata->beta_ps_bau += cs->sdata->beta_ps[j] / 10.0; 
      cs->sdata->beta_gc_bau += cs->sdata->beta_gc[j] / 10.0; 
      cs->sdata->beta_gp_bau += cs->sdata->beta_gp[j] / 10.0; 
      cs->sdata->beta_pg_bau += cs->sdata->beta_pg[j] / 10.0; 
      cs->sdata->beta_cg_bau += cs->sdata->beta_cg[j] / 10.0;
      cs->sdata->beta_cp_bau += cs->sdata->beta_cp[j] / 10.0; 
      cs->sdata->beta_pc_bau += cs->sdata->beta_pc[j] / 10.0; 
   }    
#endif  
   return 1;
}


////////////////////////////////////////////////////////////////////////////////
//! init_landuse_patches
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void init_landuse_patches (site** siteptr, UserData* data) {

   site* cs = *siteptr;

   /*init*/
   int track = 0;
   double age   = 0.0;
   double fsc   = 0.01; 
#ifdef ED
   double ssc   = 0.01;
   double psc   = 0.01;
   double msn   = 0.01;
   double fsn   = 0.01;
   double stsl  = 0.01;
   /* initialize patch water to eqm value in abs of vegetation */
   double water = (cs->sdata->soil_depth * cs->sdata->theta_max) 
                * pow(cs->sdata->precip_average / cs->sdata->k_sat, 
                      1.0 / (2.0 * cs->sdata->tau + 3.0));
#endif /* ED */
    // These two will be initialized from file values
   double area, stsc;

   for (size_t lu=1; lu<N_LANDUSE_TYPES; lu++) {
      cs->area_fraction[lu] = 0.0;
   }

   /* CROP */
   patch* newp = NULL;
#if COUPLED
   double factor = cs->sdata->grid_cell_area_total / cs->sdata->grid_cell_area;
   if (data->glm_data[0].c[cs->sdata->globY_][cs->sdata->globX_] > MIN_LANDUSE_AREA_FRACTION) {
      cs->area_fraction[LU_CROP] = data->glm_data[0].c[cs->sdata->globY_][cs->sdata->globX_] * factor;
#else
   if (data->init_c[cs->sdata->y_][cs->sdata->x_] > data->min_landuse_area_fraction) {
      cs->area_fraction[LU_CROP] = data->init_c[cs->sdata->y_][cs->sdata->x_];
#endif
      area = cs->area_fraction[LU_CROP] * data->area;
      stsc = data->init_csc[cs->sdata->y_][cs->sdata->x_];
#if defined ED
      create_patch(&cs, &newp, LU_CROP, track, age, area, 
                   water, fsc, stsc, stsl, ssc, psc, msn, fsn, data);
      init_cohorts(&newp, data);
#elif defined MIAMI_LU
      double tb = 0.0;
      create_patch(&cs, &newp, LU_CROP, track, age, area, 
                   fsc, stsc, tb, data);
#endif
      cs->youngest_patch[LU_CROP] = newp;
      cs->oldest_patch[LU_CROP] = newp;
   } 

   /* PASTURE */
   newp = NULL;
#if COUPLED
   if (data->glm_data[0].p[cs->sdata->globY_][cs->sdata->globX_] > MIN_LANDUSE_AREA_FRACTION) {
      cs->area_fraction[LU_PAST] = data->glm_data[0].p[cs->sdata->globY_][cs->sdata->globX_] * factor;
#else
   if (data->init_p[cs->sdata->y_][cs->sdata->x_] > data->min_landuse_area_fraction) {
      cs->area_fraction[LU_PAST] = data->init_p[cs->sdata->y_][cs->sdata->x_];
#endif
      area = cs->area_fraction[LU_PAST] * data->area;
      stsc = data->init_psc[cs->sdata->y_][cs->sdata->x_];
#if defined ED
      create_patch(&cs, &newp, LU_PAST, track, age, area, 
                   water, fsc, stsc, stsl, ssc, psc, msn, fsn, data);
      init_cohorts(&newp, data);
#elif defined MIAMI_LU
      double tb = data->init_pb[cs->sdata->y_][cs->sdata->x_];
      create_patch(&cs, &newp, LU_PAST, track, age, area, 
                   fsc, stsc, tb, data);
#endif
      cs->youngest_patch[LU_PAST] = newp;
      cs->oldest_patch[LU_PAST] = newp;
   } 

   /* NATURAL */
   /* adjust areas */
   cs->area_fraction[LU_NTRL] = 1.0 - cs->area_fraction[LU_CROP] - cs->area_fraction[LU_PAST];
   if (cs->area_fraction[LU_NTRL] < data->min_landuse_area_fraction) {
      cs->area_fraction[LU_NTRL] = data->min_landuse_area_fraction;
   }
   patch* cp = cs->youngest_patch[LU_NTRL];
   while (cp != NULL) {
      cp->area *= cs->area_fraction[LU_NTRL];
#if defined ED
      /* REDUCE DONOR COHORT DENSITY */
      cohort* cc = cp->shortest;
      while (cc != NULL) {
         cc->nindivs *= cs->area_fraction[LU_NTRL];
         cc = cc->taller;
      }
#endif
      cp = cp->older;
   } 
}


////////////////////////////////////////////////////////////////////////////////
//! landuse_dynamics
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void landuse_dynamics (unsigned int t, site** siteptr, UserData* data) {

   site* currents = *siteptr;

   if (t == 0) {
      if (currents->total_ag_biomass[0] > data->forest_definition) {
         currents->maintain_pasture_flag = 1;
      } else {
         currents->maintain_pasture_flag = 0;
      }
   }
   if ( (t % LANDUSE_FREQ == 0) && (t > 0) ) {
      if ( (currents->total_ag_biomass[0] > data->forest_definition) ) {
         currents->forest_harvest_flag = 1;
      } else {
         currents->forest_harvest_flag = 0;
      }

      /************************************************/
      /****            aging of crops              ****/
      /************************************************/
      /* currently crops not passed to patch dynamics, so do patch aging here */
      patch* currentp = currents->oldest_patch[LU_CROP];
      while (currentp != NULL) {
         currentp->age += data->deltat * LANDUSE_FREQ; /* patch aging */
         currentp = currentp->younger;
      }

      /************************************************/
      /*****         create new patches           *****/
      /************************************************/
      for (int lu=1; lu<N_LANDUSE_TYPES; lu++) {
#if defined ED
         create_patch(&currents, &(currents->new_patch[lu]), lu, 
                      0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      data);
#elif defined MIAMI_LU
         create_patch(&currents, &(currents->new_patch[lu]), lu, 
                      0, 0.0, 0.0, 0.0, 0.0, 0.0, data);
#endif
      }

#ifndef COUPLED
      size_t lu_year;
      if (data->year < N_LANDUSE_YEARS) {
         lu_year = data->year;
      } else {
         /* TODO: probably not the right assumption for all transitions. -justin */
         lu_year = N_LANDUSE_YEARS-1;
      }

      /************************************************/
      /*****            lu transistions            ****/
      /************************************************/
      for (int dlu=0; dlu<N_LANDUSE_TYPES; dlu++) {
         for (int tlu=1; tlu<N_LANDUSE_TYPES; tlu++) {
            if (tlu != dlu) {
               /* v2s will be dealt with in harvesting. skip */
               if ( !((dlu == LU_NTRL) && (tlu == LU_SCND)) ) {
                  double beta = currents->sdata->beta[dlu][tlu-1][lu_year];
#if 1
                   if (data->year < N_LANDUSE_YEARS) {
                       printf("LU year exceeds 2005 \n");
                       beta=0;
                   }
#endif
                  landuse_transition(&currents, data, dlu, tlu, beta);
               }
            }
         }
      }

#if 1
      for (int i=0; i<N_VBH_TYPES; i++) {
         cut_forest(&currents, data, LU_NTRL, currents->sdata->vbh[i][lu_year], i);
      }
      for (int i=0; i<N_SBH_TYPES; i++) {
         cut_forest(&currents, data, LU_SCND, currents->sdata->sbh[i][lu_year], i);
      }
#endif

#else // COUPLED
      size_t y = currents->sdata->globY_;
      size_t x = currents->sdata->globX_;
      double factor = currents->sdata->grid_cell_area_total / currents->sdata->grid_cell_area;

      landuse_transition(&currents, data, LU_NTRL, LU_CROP, data->glm_data[0].flowvc[y][x] * factor);
      landuse_transition(&currents, data, LU_NTRL, LU_PAST, data->glm_data[0].flowvp[y][x] * factor);
      landuse_transition(&currents, data, LU_SCND, LU_CROP, data->glm_data[0].flowsc[y][x] * factor);
      landuse_transition(&currents, data, LU_SCND, LU_PAST, data->glm_data[0].flowsp[y][x] * factor);
      landuse_transition(&currents, data, LU_CROP, LU_SCND, data->glm_data[0].flowcs[y][x] * factor);
      landuse_transition(&currents, data, LU_CROP, LU_PAST, data->glm_data[0].flowcp[y][x] * factor);
      landuse_transition(&currents, data, LU_PAST, LU_SCND, data->glm_data[0].flowps[y][x] * factor);
      landuse_transition(&currents, data, LU_PAST, LU_CROP, data->glm_data[0].flowpc[y][x] * factor);
      cut_forest(&currents, data, LU_NTRL, data->glm_data[0].flowvbh[y][x] * factor, 0);
      cut_forest(&currents, data, LU_NTRL, data->glm_data[0].flowvbh2[y][x] * factor, 0);
      cut_forest(&currents, data, LU_SCND, data->glm_data[0].flowsbh[y][x] * factor, 1);
      cut_forest(&currents, data, LU_SCND, data->glm_data[0].flowsbh2[y][x] * factor, 1);
      cut_forest(&currents, data, LU_SCND, data->glm_data[0].flowsbh3[y][x] * factor, 1);
#endif //COUPLED

#if 1
      harvest_croplands(&currents, data);
      graze_pastures(&currents, data);
#endif

      /************************************************/
      /****     insert or delete new patches      *****/
      /************************************************/
      for (size_t lu=1; lu<N_LANDUSE_TYPES; lu++) {
         if (currents->new_patch[lu]->area > data->min_change_in_area) {
#ifdef ED
            /* plant */
            init_cohorts(&(currents->new_patch[lu]), data);
#endif
            patch* target = currents->new_patch[lu];
            currentp = currents->youngest_patch[lu];
            if (currentp == NULL)
               currents->oldest_patch[lu] = target;
            else
               currentp->younger = target;
            target->older = currentp;
            target->younger = NULL;
            currents->youngest_patch[lu] = target;

            if(data->patch_fusion) {
               if(currents->youngest_patch[lu] != currents->oldest_patch[lu])
                  fuse_2_patches(&(currents->youngest_patch[lu]),
                                 &(currents->youngest_patch[lu]->older), 1, data);
            }
         } else {
            /* free the new patch */
#ifdef ED
            cohort* currentc = currents->new_patch[lu]->shortest;
            while (currentc != NULL) {
               cohort* tmpc = currentc->taller;
               free(currentc);
               currentc = tmpc;
            }
#endif
            if (lu == LU_SCND) {
               free(currents->new_patch[lu]->phistory);
            }
            free(currents->new_patch[lu]);
         }
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
//! landuse_transition
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void landuse_transition (site** siteptr, UserData* data,
                         int donor_lu, int target_lu, double beta) {

   site* currents = *siteptr;

   if (beta == 0.0) {
      return;
   }

   /* determine how much donor land is available */
   double total_area = 0.0;
   patch* currentp = currents->youngest_patch[donor_lu];
   while (currentp != NULL) {
      total_area += currentp->area;
      currentp = currentp->older;
   }
   /* TODO: should this be MIN_PATCH_AREA? -justin */
   if (total_area <= 0.0) {
      /*printf("no land to transition %s %d %d->%d %e %e\n",
             currents->name, data->year, donor_lu, target_lu, 
             currents->area_fraction[donor_lu], beta);*/
      return;
   }
   double donor_fraction = total_area / data->area;

   /* convert beta to fraction of donor landuse */
   if (beta > donor_fraction) {
      /*printf("not enough land %s %d %d->%d %e %e %e %e\n",
             currents->name, data->year, donor_lu, target_lu, currents->area_fraction[donor_lu],
             total_area, donor_fraction, beta);*/
      beta = 1.0;
   } else {
      beta /= donor_fraction;
   }
      
   double change_in_area = beta * total_area;

   if ( (change_in_area > data->min_landuse_change_fraction * data->area) 
        && (total_area - change_in_area > data->min_landuse_area_fraction * data->area)) {

      currentp = currents->youngest_patch[donor_lu]; 
      while (currentp != NULL) {
         patch* nextp = currentp->older;
      
         /* COMPUTE CHANGE IN AREA IN EACH DONOR PATCH *
          * ASSUMED PROPORTIONAL TO AREAS              */      
         /*change_in_area = beta * (currents->area_fraction[donor_lu] * data->area)
          * currentp->area / total_area;*/
         change_in_area = beta * currentp->area;
      
         if (change_in_area > data->min_change_in_area_patch) {
            /* INCREASE NEW PATCH AREA */
            currents->new_patch[target_lu]->area += change_in_area;
   
            /* ACCUMULATE LITTER */
            /* abondoned cropland (to pasture or secondary) abondoned pasture *
             * (to secondary) have complete survivorship. no litter to        *
             * accumulate. skip for any transition from crop, or any          *
             * transition from pasture except to crop                         */
            if ( (donor_lu != LU_CROP) && 
                 ( (donor_lu != LU_PAST) || (target_lu == LU_CROP) ) ) {
               accumulate_litter_from_disturbance(&(currents->new_patch[target_lu]),
                                                  &currentp, change_in_area, 9, data);
            }

#ifdef ED
            /* REDUCE DONOR COHORT DENSITY */
            cohort* currentc = currentp->shortest;
            while (currentc != NULL) {
               currentc->nindivs -= currentc->nindivs * change_in_area / currentp->area;
               currentc = currentc->taller;
            }
#endif

            /* AGGREGATE IN PATCH STATE VARIABLES (repro Water, N, etc...) */
            aggregate_in_soil_state_from_disturbance(&(currents->new_patch[target_lu]),
                                                   &currentp, change_in_area, data);
   
            /* update area of donor patch */ 
            currentp->area -= change_in_area;
   
            if (currentp->area < data->min_patch_area) {
               /*printf("landuse transition: warning, patch area %e (%d -> %d) in %s, deleting\n",
                 currentp->area, donor_lu, target_lu, currents->name );*/
               /* TODO: this area should be added back somewhere - justin */
               terminate_patch(&currentp);
            }   
         } else {
            /*printf("patch transition too small %s %d %d->%d %e %e\n",
                   currents->name, data->year, donor_lu, target_lu, 
                   change_in_area, beta);*/
         }
         currentp = nextp;      
      }
      // Divide aggregate variables (soil C, N etc.) by area of new path 
      // to obtain averages
      (currents->new_patch[target_lu])->fast_soil_C /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->structural_soil_C /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->fast_soil_N /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->slow_soil_C /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->passive_soil_C /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->mineralized_soil_N /=(currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->structural_soil_L /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->water /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->theta /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->rh /= (currents->new_patch[target_lu])->area;
   } else {
      /*printf("site transition too small %s %d %d->%d %e %e\n",
             currents->name, data->year, donor_lu, target_lu, 
             change_in_area, beta);*/
   }
      
}

/* cut forest based on fraction */
void cut_forest (site** siteptr, UserData* data, 
                 int donor_lu, double beta, int bh_type) {

   int target_lu = LU_SCND;

   site* currents=*siteptr;

   if (beta == 0.0) {
      return;
   }

   /* determine how much donor land is available */
   double total_area = 0.0;
   patch* currentp = currents->youngest_patch[donor_lu];
   while (currentp != NULL) {
      total_area += currentp->area;
      currentp = currentp->older;
   }
   /* TODO: should this be MIN_PATCH_AREA? -justin */
   if (total_area <= 0.0) {
      /*printf("no land to transition %s %d %d->%d %e %e\n",
             currents->name, data->year, donor_lu, target_lu, 
             currents->area_fraction[donor_lu], beta);*/
      return;
   }
   double donor_fraction = total_area / data->area;

   /* convert beta to fraction of donor landuse */
   if (beta > donor_fraction) {
      /*printf("not enough land %s %d %d->%d %e %e %e %e\n",
             currents->name, data->year, donor_lu, target_lu, currents->area_fraction[donor_lu],
             total_area, donor_fraction, beta);*/
      beta = 1.0;
   } else {
      beta /= donor_fraction;
   }
      
   double change_in_area = beta * total_area;

   if ( (change_in_area > data->min_landuse_change_fraction * data->area) 
        && (total_area - change_in_area > data->min_landuse_area_fraction * data->area)) {

      currentp = currents->youngest_patch[donor_lu]; 
      while (currentp != NULL) {
         patch* nextp = currentp->older;
      
         /*COMPUTE CHANGE IN AREA IN EACH DONOR PATCH*/
         /*ASSUMED PROPORTIONAL TO AREAS*/      
         change_in_area = beta * currentp->area;
      
         if (change_in_area > data->min_change_in_area_patch) {
            /* INCREASE NEW PATCH AREA */
            currents->new_patch[target_lu]->area += change_in_area;
   
            /* ACCUMULATE LITTER */
            accumulate_litter_from_disturbance(&(currents->new_patch[target_lu]),
                                               &currentp, change_in_area, 9, data);

#ifdef ED
            /* REDUCE DONOR COHORT DENSITY */
            cohort* currentc = currentp->shortest;
            while (currentc != NULL) {
               currentc->nindivs -= currentc->nindivs * change_in_area / currentp->area;
               currentc = currentc->taller;
            }
#endif

            /* AGGREGATE IN PATCH STATE VARIABLES (repro Water, N, etc...) */
            aggregate_in_soil_state_from_disturbance(&(currents->new_patch[target_lu]),
                                                   &currentp, change_in_area, data);
   
            /* update area of donor patch */ 
            currentp->area -= change_in_area;
            currents->area_harvested[donor_lu][bh_type] += change_in_area;
            currents->biomass_harvested[donor_lu][bh_type] += change_in_area*currentp->total_ag_biomass;
   
            if (currentp->area < data->min_patch_area) {
               /*printf("landuse transition: warning, patch area %e (%d -> %d) in %s, resetting\n",
                 currentp->area, donor_lu, target_lu, currents->name );*/
               /* TODO: this area should be added back somewhere */
               terminate_patch(&currentp);
            }   
         } else {
            /*printf("patch transition too small %s %d %d->%d %e %e\n",
                   currents->name, data->year, donor_lu, target_lu, 
                   change_in_area, beta);*/
         }
         currentp = nextp;
      }
      
      // Divide aggregate variables (soil C, N etc.) by area of new path 
      // to obtain averages
      (currents->new_patch[target_lu])->fast_soil_C /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->structural_soil_C /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->fast_soil_N /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->slow_soil_C /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->passive_soil_C /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->mineralized_soil_N /=(currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->structural_soil_L /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->water /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->theta /= (currents->new_patch[target_lu])->area;
      (currents->new_patch[target_lu])->rh /= (currents->new_patch[target_lu])->area;
   } else {
      /*printf("site transition too small %s %d %d->%d %e %e\n",
             currents->name, data->year, donor_lu, target_lu, 
             change_in_area, beta);*/
   }
      
}


////////////////////////////////////////////////////////////////////////////////
//! update_landuse
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void update_landuse (site* siteptr, UserData& data) {

   for (size_t lu=0; lu<N_LANDUSE_TYPES; lu++) {
      double area = 0.0;
      patch* cp = siteptr->oldest_patch[lu];
      while (cp != NULL) {
         area += cp->area;
         cp = cp->younger;
      }
      siteptr->area_fraction[lu] = area / data.area;
   }
}


////////////////////////////////////////////////////////////////////////////////
//! harvest_croplands
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void harvest_croplands (site** siteptr, UserData* data) {
  
   site* currents = *siteptr;
   patch* currentp = currents->youngest_patch[LU_CROP];
   while (currentp != NULL) {
      double harvested_structural_C = 0.0;
      double harvested_fast_C       = 0.0;
#if defined ED
      double harvested_fast_N       = 0.0;

      cohort* currentc = currentp->tallest;
      while (currentc != NULL) {
         /*harvest everything*/
         harvested_structural_C = currentc->bdead * currentc->nindivs;
         harvested_fast_C = currentc->balive * currentc->nindivs;
         harvested_fast_N = (1.0 / data->c2n_leaf[currentc->species])
            * currentc->balive * currentc->nindivs;  
         currentc->nindivs = 0.00001;
    
         currentc = currentc->shorter;
      }

      /* harvest everything */
      /* put a fraction of the biomass back on the harvested land */
      currentp->structural_soil_C += data->crop_residue * harvested_structural_C
         / currentp->area + (1.0 - data->fraction_balive_2_fast)
         * data->crop_residue * harvested_fast_C / currentp->area;
      currentp->fast_soil_C += data->fraction_balive_2_fast 
         * data->crop_residue * harvested_fast_C / currentp->area;
      currentp->structural_soil_L += (data->l2n_stem / data->c2n_stem)
         * data->crop_residue * harvested_structural_C / currentp->area 
         + (data->l2n_stem / data->c2n_stem) * (1.0-data->fraction_balive_2_fast)
         * data->crop_residue * harvested_fast_C / currentp->area;
      currentp->fast_soil_N += data->fraction_balive_2_fast 
         * data->crop_residue * harvested_fast_N / currentp->area;
      /*plant new crop*/
      init_cohorts(&currentp, data);
#elif defined MIAMI_LU
      harvested_structural_C = currentp->total_biomass;
      /*currentp->fast_soil_C += data->fraction_balive_2_fast 
          * data->crop_residue * harvested_fast_C; */
      currentp->structural_soil_C += data->crop_residue * harvested_structural_C 
         + (1.0 - data->fraction_balive_2_fast) * data->crop_residue * harvested_fast_C;    

      /*resest total biomass*/
      currentp->total_biomass = 0.0;
      currentp->total_ag_biomass = 0.0;
#endif /* ED V. MIAMI_LU */

      currentp = currentp->older;
   }
}


////////////////////////////////////////////////////////////////////////////////
//! graze_pastures
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void graze_pastures (site** siteptr, UserData* data) {

   site* currents =* siteptr;
   patch* currentp = currents->youngest_patch[LU_PAST];
   while (currentp != NULL) {
      double grazed_structural_C = 0.0;
      double grazed_fast_C       = 0.0;
#if defined ED
      double grazed_fast_N       = 0.0;

      cohort* currentc = currentp->tallest;
      while (currentc != NULL) {
         /* 1. graze grasses */
         if (data->is_grass[currentc->species]) {
            grazed_structural_C = data->grazing_intensity * currentc->bdead * currentc->nindivs;
            grazed_fast_C = data->grazing_intensity * currentc->balive * currentc->nindivs;
            grazed_fast_N = (1.0 / data->c2n_leaf[currentc->species])
               * data->grazing_intensity * currentc->balive * currentc->nindivs;
            currentc->nindivs *= (1.0 - data->grazing_intensity);
         }
         /* 2. remove trees explicitly if primary vegetation forested*/
         if ( (currents->maintain_pasture_flag == 1) && (!data->is_grass[currentc->species]) ) {
            grazed_structural_C = currentc->bdead * currentc->nindivs;
            grazed_fast_C = currentc->balive * currentc->nindivs;
            grazed_fast_N = (1.0 / data->c2n_leaf[currentc->species])
               * currentc->balive * currentc->nindivs;
            currentc->nindivs = 0.000001;
         }

         currentc = currentc->shorter;
      }

      /* put a fraction of the biomass back on the harvested land */
      currentp->fast_soil_C += data->fraction_balive_2_fast
         * data->grazing_residue * grazed_fast_C / currentp->area;
      currentp->structural_soil_C += data->grazing_residue 
         * grazed_structural_C / currentp->area 
         + (1.0 - data->fraction_balive_2_fast) * data->grazing_residue
         * grazed_fast_C / currentp->area;
      currentp->structural_soil_L +=  (data->l2n_stem / data->c2n_stem)
         * data->grazing_residue * grazed_structural_C / currentp->area 
         + (data->l2n_stem / data->c2n_stem) 
         * (1.0 - data->fraction_balive_2_fast) 
         * data->grazing_residue * grazed_fast_C / currentp->area;
      currentp->fast_soil_N += data->fraction_balive_2_fast 
         * data->grazing_residue * grazed_fast_N / currentp->area;

#elif defined MIAMI_LU
      /* remove trees explicitly if primary vegetation forested */
      if (currents->maintain_pasture_flag == 1) 
         grazed_structural_C = currentp->total_biomass;
      else /* graze grasses */
         grazed_structural_C = data->grazing_intensity * currentp->total_biomass;

      /* reduce biomass accordingly*/
      currentp->total_biomass -= grazed_structural_C;
      currentp->total_ag_biomass -= data->agf_biomass * grazed_structural_C;

      /* put a fraction of the biomass back on the harvested land */
      /* currentp->fast_soil_C += data->fraction_balive_2_fast
            * data->grazing_residue * grazed_fast_C; */
      currentp->structural_soil_C += data->grazing_residue * grazed_structural_C 
         + (1.0 - data->fraction_balive_2_fast) * data->grazing_residue * grazed_fast_C;
#endif

      currentp = currentp->older;
   }
}


////////////////////////////////////////////////////////////////////////////////
//! print_landuse
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void print_landuse (unsigned int time, site** siteptr, UserData* data) {

   site* cs = *siteptr;
  
   /* open print files */
   char filename[STR_LEN];
   strcpy(filename, data->base_filename);
   strcat(filename, ".");
   strcat(filename, cs->sdata->name_);
   strcat(filename, ".landuse");

   FILE *outfile;
   if (time == 0) { 
      outfile = fopen(filename, "w");
   } else {
      outfile = fopen(filename, "a");
   }
   /* print to files */
   fprintf(outfile, 
           "%s t= %f natural %f crop %f past %f secondary %f total %f\n",
           cs->sdata->name_, time * TIMESTEP,
           cs->area_fraction[LU_NTRL],
           cs->area_fraction[LU_CROP],
           cs->area_fraction[LU_PAST],
           cs->area_fraction[LU_SCND],
           (cs->area_fraction[LU_NTRL] + cs->area_fraction[LU_CROP]
            + cs->area_fraction[LU_PAST] + cs->area_fraction[LU_SCND]));
  
   fclose(outfile);
}


////////////////////////////////////////////////////////////////////////////////
//! update_mean_age_secondary_lands
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
void update_mean_age_secondary_lands (site** siteptr, UserData* data) {

   site* cs = *siteptr;

   cs->mean_AGE_sec = 0.0;
 
#if 0 // verify this works with coupler years (or fix)
   double area_scnd = 0.0;
   /*for each AGE*/
   for (size_t age=0; age<data->year+1; age++) {
      size_t age_index = data->year - age;  
       
      patch* cp = cs->youngest_patch[LU_SCND];
      while (cp != NULL) {
         cs->mean_AGE_sec += age * (*(cp->phistory + age_index));
         area_scnd += (*(cp->phistory + age_index));
         cp = cp->older;
      }
   }

   if (area_scnd > 0.0) {
      cs->mean_AGE_sec /= area_scnd;
   } else {
      cs->mean_AGE_sec = 0.0;
   }
#endif
}


/**********************************************************************/
#endif /* LANDUSE */
