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
void plant_croplands (site** siteptr, UserData* data);
void graze_pastures (site** siteptr, UserData* data);
void wood_pool_decay(site** siteptr, UserData* data);   /// Function to compute emission from wood product pools, wrote by Lei Ma


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
    int rv, ncid, varid, dlu, tlu, i,ylu;
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

    for (dlu=0; dlu<N_LANDUSE_TYPES; dlu++) {
        for (tlu=1; tlu<N_LANDUSE_TYPES; tlu++) {
            for (ylu=0;ylu<N_LANDUSE_YEARS;ylu++)
            {
                if ((dlu != tlu) && !( (dlu == LU_NTRL) && (tlu == LU_SCND) ) ) {
                    if ((dn = lu2charname(dlu)) && (tn = lu2charname(tlu))) {
                        cs->sdata->beta[dlu][tlu-1][ylu]=data->gfl[dlu][tlu-1][ylu][cs->sdata->globY_][cs->sdata->globX_];
                    }
                }
                
            }
            /* skip transitions self->self and v->s (v->s dealt with in sbh/vbh) */
        }
    }
    
    for (dlu=0; dlu<N_VBH_TYPES; dlu++) {
        for (ylu=0;ylu<N_LANDUSE_YEARS;ylu++)
        {
            cs->sdata->vbh[dlu][0]=data->gfvh[dlu][ylu][cs->sdata->globY_][cs->sdata->globX_];
        }
    }
    for (dlu=0; dlu<N_SBH_TYPES; dlu++) {
        for (ylu=0;ylu<N_LANDUSE_YEARS;ylu++)
        {
            cs->sdata->sbh[dlu][ylu]=data->gfsh[dlu][ylu][cs->sdata->globY_][cs->sdata->globX_];
        }
    }
    
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
       
       ///CarbonConserve
       /// Reset carbon emission variables in patches as zero
       for (int lu=0; lu<N_LANDUSE_TYPES; lu++)
       {
           currents->new_patch[lu] = NULL; ///Here, reset new_patch as NULL since the landuse type of newpath of primary forest is 3, there may be bug in patch_dynamic when creating new patch --Lei
           currentp = currents->youngest_patch[lu];
           while (currentp != NULL)
           {
               currentp->forest_harvested_c = 0.0;
//               currentp->yr1_decay_product_pool = 0.0;
//               currentp->yr10_decay_product_pool = 0.0;
//               currentp->yr100_decay_product_pool = 0.0;
               currentp->past_harvested_c = 0.0;
               currentp->crop_harvested_c = 0.0;
               currentp->product_emission = 0.0;
               currentp = currentp->older;
           }
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
      if (data->mechanism_year <= LANDUSE_END) {
         lu_year = data->mechanism_year-LANDUSE_START-1;
      } else {
         /* TODO: probably not the right assumption for all transitions. -justin */
          ///Carbon Conserve
          /// There is something wrong with LUH transition data as all transition is zero at year 506
          /// So I change the lu_year to the last second year than last year to avoid zero LU transtion. -- Lei
         lu_year = LANDUSE_END-LANDUSE_START-1;
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
//#if 1
//                   if (data->year > N_LANDUSE_YEARS) {
//                       //printf("LU year exceeds 2005 %d %d, using 2005 LU transition instead\n",data->year,N_LANDUSE_YEARS);
//                       beta=currents->sdata->beta[dlu][tlu-1][lu_year];
//                   }
//#endif
                   
#if CHECK_C_CONSERVE
                   double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0;
                   double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0;
                   double all_f_hrC_before = 0.0, all_p_hrC_before = 0.0, all_c_hrC_before = 0.0;
                   double all_f_hrC_after = 0.0, all_p_hrC_after = 0.0, all_c_hrC_after = 0.0;
                   for (int lu=0; lu<N_LANDUSE_TYPES; lu++)
                   {
                       currentp = currents->youngest_patch[lu];
                       cohort* mlcc = NULL;
                       while(currentp!=NULL)
                       {
                           mlcc = currentp->shortest;
                           double tmp_tb = 0.0;
                           while (mlcc!=NULL) {
                               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
                               mlcc = mlcc->taller;
                           }
                           all_tb_before += tmp_tb*currentp->area/data->area;
                           all_sc_before += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
                           all_f_hrC_before += currentp->forest_harvested_c * currentp->area/data->area;
                           all_p_hrC_before += currentp->past_harvested_c * currentp->area/data->area;
                           all_c_hrC_before += currentp->crop_harvested_c * currentp->area/data->area;
//                           printf("ck old_patch_bf LU %d age %.4f area %.5f tc %.15f tb %.15f sc %.15f f_hr %.15f p_hr %.15f c_hr %.15f all_tc %.15f\n",currentp->landuse,currentp->age,currentp->area,
//                                  currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C+tmp_tb,tmp_tb,
//                                  currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C,currentp->forest_harvested_c,currentp->past_harvested_c,currentp->crop_harvested_c,all_sc_before);
                           currentp = currentp->older;
                       }
                       currentp = currents->new_patch[lu];
                       if (currentp!=NULL)
                       {
                           mlcc = currentp->shortest;
                           double tmp_tb = 0.0;
                           while (mlcc!=NULL) {
                               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
                               mlcc = mlcc->taller;
                           }
                           all_tb_before += tmp_tb*currentp->area/data->area;
                           all_sc_before += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
                           all_f_hrC_before += currentp->forest_harvested_c * currentp->area/data->area;
                           all_p_hrC_before += currentp->past_harvested_c * currentp->area/data->area;
                           all_c_hrC_before += currentp->crop_harvested_c * currentp->area/data->area;
                           //                           printf("ck new_patch_bf LU %d age %.4f area %.5f tc %.15f tb %.15f sc %.15f f_hr %.15f p_hr %.15f c_hr %.15f all_tc %.15f\n",currentp->landuse,currentp->age,currentp->area,
                           //                                  currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C+tmp_tb,tmp_tb,
                           //                                  currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C,currentp->forest_harvested_c,currentp->past_harvested_c,currentp->crop_harvested_c,all_sc_before);
                       }
                   }
                   all_tc_before = all_tb_before + all_sc_before + (all_f_hrC_before + all_p_hrC_before + all_c_hrC_before)*data->deltat;
#endif
                   landuse_transition(&currents, data, dlu, tlu, beta);
                   
#if CHECK_C_CONSERVE
                   for (int lu=0; lu<N_LANDUSE_TYPES; lu++)
                   {
                       currentp = currents->youngest_patch[lu];
                       cohort* mlcc = NULL;
                       while(currentp!=NULL)
                       {
                           mlcc = currentp->shortest;
                           double tmp_tb = 0.0;
                           while (mlcc!=NULL) {
                               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
                               mlcc = mlcc->taller;
                           }
                           all_tb_after += tmp_tb*currentp->area/data->area;
                           all_sc_after += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
                           all_f_hrC_after += currentp->forest_harvested_c * currentp->area/data->area;
                           all_p_hrC_after += currentp->past_harvested_c * currentp->area/data->area;
                           all_c_hrC_after += currentp->crop_harvested_c * currentp->area/data->area;
//                           printf("ck old_patch_af LU %d age %.4f area %.5f tc %.15f tb %.15f sc %.15f f_hr %.15f p_hr %.15f c_hr %.15f all_tc %.15f\n",currentp->landuse,currentp->age,currentp->area,
//                                  currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C+tmp_tb,tmp_tb,
//                                  currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C,currentp->forest_harvested_c,currentp->past_harvested_c,currentp->crop_harvested_c,all_sc_after);
                           currentp = currentp->older;
                       }
                       currentp = currents->new_patch[lu];
                       if (currentp!=NULL)
                       {
                           mlcc = currentp->shortest;
                           double tmp_tb = 0.0;
                           while (mlcc!=NULL) {
                               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
                               mlcc = mlcc->taller;
                           }
                           all_tb_after += tmp_tb*currentp->area/data->area;
                           all_sc_after += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
                           all_f_hrC_after += currentp->forest_harvested_c * currentp->area/data->area;
                           all_p_hrC_after += currentp->past_harvested_c * currentp->area/data->area;
                           all_c_hrC_after += currentp->crop_harvested_c * currentp->area/data->area;
                           //                           printf("ck new_patch_af LU %d age %.4f area %.5f tc %.15f tb %.15f sc %.15f f_hr %.15f p_hr %.15f c_hr %.15f all_tc %.15f\n",currentp->landuse,currentp->age,currentp->area,
                           //                                  currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C+tmp_tb,tmp_tb,
                           //                                  currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C,currentp->forest_harvested_c,currentp->past_harvested_c,currentp->crop_harvested_c,all_sc_after);
                       }
                   }
                   all_tc_after = all_tb_after + all_sc_after + (all_f_hrC_after + all_p_hrC_after + all_c_hrC_after)*data->deltat;
                   
                   if (abs(all_tc_after-all_tc_before)>1e-9)
                   {
                       printf("Carbon leakage in landuse_transition  : imbalance   %.15f site_tc_af %.15f site_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
                       printf("                                      : site_sc_bf  %.15f site_tb_bf %.15f \n",all_sc_before,all_tb_before);
                       printf("                                      : site_sc_af  %.15f site_tb_af %.15f \n",all_sc_after,all_tb_after);
                       printf("                                      : site_f_hr_bf  %.15f site_p_hr_bf %.15f site_c_hr_bf %.15f\n",all_f_hrC_before,all_p_hrC_before,all_c_hrC_before);
                       printf("                                      : site_f_hr_af  %.15f site_p_hr_af %.15f site_c_hr_af %.15f\n",all_f_hrC_after,all_p_hrC_after,all_c_hrC_after);
                       printf("                                      : dlu %2d tlu %2d beta %.15f\n",dlu,tlu,beta);
                       printf(" --------------------------------------------------------------------------------------\n");
                   }
#endif
                   
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
      //harvest_croplands(&currents, data);
#if CHECK_C_CONSERVE
      ///CarbonConserve
       double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0;
       double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0;
       double all_c_hrC_bf = 0.0, all_f_hrC_bf = 0.0, all_p_hrC_bf = 0.0;
       double all_c_hrC_af = 0.0, all_f_hrC_af = 0.0, all_p_hrC_af = 0.0;
       currentp = currents->youngest_patch[LU_PAST];
       cohort* mlcc = NULL;
       while(currentp!=NULL)
       {
           mlcc = currentp->shortest;
           double tmp_tb = 0.0;
           while (mlcc!=NULL) {
               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
               mlcc = mlcc->taller;
           }
           all_tb_before += tmp_tb*currentp->area/data->area;
           all_sc_before += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
           all_p_hrC_bf += currentp->past_harvested_c * currentp->area/data->area;
           currentp = currentp->older;
       }
       all_tc_before = all_tb_before + all_sc_before;
#endif
      graze_pastures(&currents, data);
#endif
       
#if CHECK_C_CONSERVE
       ///CarbonConserve
       currentp = currents->youngest_patch[LU_PAST];
       while(currentp!=NULL)
       {
           mlcc = currentp->shortest;
           double tmp_tb = 0.0;
           while (mlcc!=NULL) {
               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
               mlcc = mlcc->taller;
           }
           all_tb_after += tmp_tb*currentp->area/data->area;
           all_sc_after += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
           all_p_hrC_af += currentp->past_harvested_c * currentp->area/data->area;
           currentp = currentp->older;
       }
       all_tc_after = all_tb_after + all_sc_after + all_p_hrC_af*data->deltat;
       
       if (abs(all_tc_after-all_tc_before)>1e-9)
       {
           printf("Carbon leakage in graze_pastures  : imbalance   %.15f past_tc_af %.15f past_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
           printf("                                  : past_sc_bf  %.15f past_tb_bf %.15f \n",all_sc_before,all_tb_before);
           printf("                                  : past_sc_af  %.15f past_tb_af %.15f all_crop_harvest_af %.15f\n",all_sc_after,all_tb_after,all_c_hrC_af);
           printf(" --------------------------------------------------------------------------------------\n");
       }
#endif

      /************************************************/
      /****     insert or delete new patches      *****/
      /************************************************/
      for (size_t lu=1; lu<N_LANDUSE_TYPES; lu++) {
#if CHECK_C_CONSERVE
          ///CarbonConserve
          double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0;
          double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0;
          double all_f_hrC_bf = 0.0, all_p_hrC_bf = 0.0, all_c_hrC_bf = 0.0;
          double all_f_hrC_af = 0.0, all_p_hrC_af = 0.0, all_c_hrC_af = 0.0;
          currentp = currents->youngest_patch[lu];
          cohort* mlcc = NULL;
          while(currentp!=NULL)
          {
              mlcc = currentp->shortest;
              double tmp_tb = 0.0;
              while (mlcc!=NULL) {
                  tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
                  mlcc = mlcc->taller;
              }
              all_tb_before += tmp_tb*currentp->area/data->area;
              all_sc_before += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
              all_f_hrC_bf += currentp->forest_harvested_c * currentp->area/data->area;
              all_p_hrC_bf += currentp->past_harvested_c * currentp->area/data->area;
              all_c_hrC_bf += currentp->crop_harvested_c*currentp->area/data->area;
              currentp = currentp->older;
          }
          currentp = currents->new_patch[lu];
          if (currentp!=NULL)
          {
              mlcc = currentp->shortest;
              double tmp_tb = 0.0;
              while (mlcc!=NULL) {
                  tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
                  mlcc = mlcc->taller;
              }
              all_tb_before += tmp_tb*currentp->area/data->area;
              all_sc_before += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
              all_f_hrC_bf += currentp->forest_harvested_c * currentp->area/data->area;
              all_p_hrC_bf += currentp->past_harvested_c * currentp->area/data->area;
              all_c_hrC_bf += currentp->crop_harvested_c*currentp->area/data->area;
          }
          all_tc_before = all_tb_before + all_sc_before;
#endif
          
         if (currents->new_patch[lu]->area > data->min_change_in_area)
         {
#ifdef ED
            /* plant */
             //Only plant for non-crop types as crop is planted in plant_crops function which is called based on crop calendar.
             if (lu != LU_CROP)
             {
                 init_cohorts(&(currents->new_patch[lu]), data);
             }

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
          
#if CHECK_C_CONSERVE
          ///CarbonConserve
          currentp = currents->youngest_patch[lu];
          while(currentp!=NULL)
          {
              mlcc = currentp->shortest;
              double tmp_tb = 0.0;
              while (mlcc!=NULL) {
                  tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
                  mlcc = mlcc->taller;
              }
              all_tb_after += tmp_tb*currentp->area/data->area;
              all_sc_after += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
              all_f_hrC_af += currentp->forest_harvested_c * currentp->area/data->area;
              all_p_hrC_af += currentp->past_harvested_c * currentp->area/data->area;
              all_c_hrC_af += currentp->crop_harvested_c*currentp->area/data->area;
              currentp = currentp->older;
          }
          all_tc_after = all_tb_after + all_sc_after;
          
          if (abs(all_tc_after-all_tc_before)>1e-9)
          {
              printf("Carbon leakage in ini_patch of LU-%d: imbalance   %.15f site_tc_af %.15f site_tc_bf %.15f\n",lu,all_tc_after-all_tc_before,all_tc_after,all_tc_before);
              printf("                                    : site_sc_bf  %.15f site_tb_bf %.15f all_f_hrC_bf %.15f all_p_hrC_bf %.15f all_c_hrC_bf %.15f\n",all_sc_before,all_tb_before,all_f_hrC_bf,all_p_hrC_bf,all_c_hrC_bf);
              printf("                                    : site_sc_af  %.15f site_tb_af %.15f all_f_hrC_af %.15f all_p_hrC_af %.15f all_c_hrC_af %.15f\n",all_sc_after,all_tb_after,all_f_hrC_af,all_p_hrC_af,all_c_hrC_af);
              printf(" --------------------------------------------------------------------------------------\n");
          }
#endif
      }
       wood_pool_decay(&currents, data);
   }
    if (data->planting_probability[data->time_period][currents->sdata->globY_][currents->sdata->globX_]>0.8)
    {
        plant_croplands(&currents, data);
    }

    if (data->harvest_probability[data->time_period][currents->sdata->globY_][currents->sdata->globX_]>0.8)
    {
        harvest_croplands(&currents, data);
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
       
#if CHECK_C_CONSERVE
       ///CarbonConserve
       double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0, all_f_harv_before = 0.0, all_p_harv_before = 0.0, all_c_harv_before = 0.0;
       double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0, all_f_harv_after = 0.0, all_p_harv_after = 0.0, all_c_harv_after = 0.0;
       double tmp_sc1 = 0.0, tmp_sc2 = 0.0, tmp_sc3 = 0.0;
       currentp = currents->youngest_patch[donor_lu];
       cohort* mlcc = NULL;
       while(currentp!=NULL)
       {
           mlcc = currentp->shortest;
           double tmp_tb = 0.0;
           while (mlcc!=NULL) {
               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
               mlcc = mlcc->taller;
           }
           all_tb_before += tmp_tb*currentp->area/data->area;
           all_sc_before += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
           all_f_harv_before += currentp->forest_harvested_c * currentp->area/data->area;
           all_p_harv_before += currentp->past_harvested_c * currentp->area/data->area;
           all_c_harv_before += currentp->crop_harvested_c * currentp->area/data->area;
           currentp = currentp->older;
       }
       currentp = currents->new_patch[target_lu];
       tmp_sc1 = all_sc_before;
       double tmp_tb = 0.0;
       mlcc = currentp->shortest;
       while (mlcc!=NULL) {
           tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
           mlcc = mlcc->taller;
       }
       all_tb_before += tmp_tb*currentp->area/data->area;
       all_sc_before += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
       all_f_harv_before += currentp->forest_harvested_c * currentp->area/data->area;
       all_p_harv_before += currentp->past_harvested_c * currentp->area/data->area;
       all_c_harv_before += currentp->crop_harvested_c * currentp->area/data->area;
       all_tc_before = all_tb_before + all_sc_before + (all_f_harv_before+all_p_harv_before+all_c_harv_before)*data->deltat;
#endif
       ///CarbonConserve
       double old_fast_soil_C = (currents->new_patch[target_lu])->fast_soil_C;
       double old_structural_soil_C = (currents->new_patch[target_lu])->structural_soil_C;
       double old_fast_soil_N = (currents->new_patch[target_lu])->fast_soil_N;
       double old_slow_soil_C = (currents->new_patch[target_lu])->slow_soil_C;
       double old_passive_soil_C = (currents->new_patch[target_lu])->passive_soil_C;
       double old_mineralized_soil_N = (currents->new_patch[target_lu])->mineralized_soil_N;
       double old_structural_soil_L = (currents->new_patch[target_lu])->structural_soil_L;
       double old_water = (currents->new_patch[target_lu])->water;
       double old_theta = (currents->new_patch[target_lu])->theta;
       double old_rh = (currents->new_patch[target_lu])->rh;
       double old_area = (currents->new_patch[target_lu])->area;
       double old_gpp_avg = (currents->new_patch[target_lu])->gpp_avg;
       double old_npp_avg = (currents->new_patch[target_lu])->npp_avg;
       double old_rh_avg = (currents->new_patch[target_lu])->rh_avg;
       double old_fire_emission = (currents->new_patch[target_lu])->fire_emission;
       double old_product_emission = (currents->new_patch[target_lu])->product_emission;
       double old_forest_harvested_c = (currents->new_patch[target_lu])->forest_harvested_c;
       double old_yr1_decay_product_pool = (currents->new_patch[target_lu])->yr1_decay_product_pool;
       double old_yr10_decay_product_pool = (currents->new_patch[target_lu])->yr10_decay_product_pool;
       double old_yr100_decay_product_pool = (currents->new_patch[target_lu])->yr100_decay_product_pool;
       double old_past_harvested_c = (currents->new_patch[target_lu])->past_harvested_c;
       double old_crop_harvested_c = (currents->new_patch[target_lu])->crop_harvested_c;
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
                 ( (donor_lu != LU_PAST) || (target_lu == LU_CROP) ) )
            {
                accumulate_litter_from_disturbance(&(currents->new_patch[target_lu]),
                                                  &currentp, change_in_area, 9, data);
            }
             else
             {
                 ///CarbonConserve
                 ///The original code does nonthing for transition from abondoned pasture and cropland,
                 ///but just reduce the cohort density in donor patch. I think this does not mean
                 /// "abondoned cropland (to pasture or secondary) abondoned pasture have complete survivorship" as the reduced cohort are not
                 /// put in new patch, and this will lead to carbon leakage.
                 /// The following code put the reduced cohort to new patch using the same way as patch_dynamic did.
                 cohort* currentc = currentp->shortest;
                 while (currentc != NULL) {
                     // make new cohort
                     cohort* newcohort = (cohort*) malloc(sizeof(cohort));
                     if (newcohort == NULL) {
                         printf("pd: out of memory\n");
                         while(getchar() != '\n');
                     }
                     
                     // copy cohort
                     copy_cohort(&currentc,&newcohort);
                     
                     // surviving indivs
                     newcohort->nindivs = currentc->nindivs * change_in_area / currentp->area;
                     newcohort->dndt = newcohort->nindivs*currentc->dndt/currentc->nindivs;
                     // insert
                     insert_cohort(&newcohort, &currents->new_patch[target_lu]->tallest, &currents->new_patch[target_lu]->shortest, data);
                     newcohort->patchptr = currents->new_patch[target_lu];
                     currentc = currentc->taller;
                 } /* end loop over cohorts */
                 for (size_t i=0;i<NSPECIES;i++) {
                     currents->new_patch[target_lu]->repro[i] += currentp->repro[i]*change_in_area/currentp->area;
                 }
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
               terminate_patch(&currentp,data);
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
       /// This below code is problematic -- Lei
//      (currents->new_patch[target_lu])->fast_soil_C /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->structural_soil_C /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->fast_soil_N /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->slow_soil_C /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->passive_soil_C /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->mineralized_soil_N /=(currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->structural_soil_L /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->water /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->theta /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->rh /= (currents->new_patch[target_lu])->area;
       
       ///CarbonConserve
       /// The above code for averaging the lastes soil state is problematic,
       /// Because before accumualting new soil state from currentp, the unit of soil state in new_patch is density-based,
       /// but aggregate_in_soil_state_from_disturbance function will accumulate absolute state (not density but total amount) into newpatch, there to convert the latest soil state variable
       /// into density-base unit again, we should use the following code rather than the above -- Lei
       double new_area =currents->new_patch[target_lu]->area;
       (currents->new_patch[target_lu])->fast_soil_C = (old_fast_soil_C*old_area+(currents->new_patch[target_lu]->fast_soil_C-old_fast_soil_C))/new_area;
       (currents->new_patch[target_lu])->structural_soil_C = (old_structural_soil_C*old_area+(currents->new_patch[target_lu]->structural_soil_C-old_structural_soil_C))/new_area;
       (currents->new_patch[target_lu])->fast_soil_N = (old_fast_soil_N*old_area+(currents->new_patch[target_lu]->fast_soil_N-old_fast_soil_N))/new_area;
       (currents->new_patch[target_lu])->slow_soil_C = (old_slow_soil_C*old_area+(currents->new_patch[target_lu]->slow_soil_C-old_slow_soil_C))/new_area;
       (currents->new_patch[target_lu])->passive_soil_C = (old_passive_soil_C*old_area+(currents->new_patch[target_lu]->passive_soil_C-old_passive_soil_C))/new_area;
       (currents->new_patch[target_lu])->mineralized_soil_N = (old_mineralized_soil_N*old_area+(currents->new_patch[target_lu]->mineralized_soil_N-old_mineralized_soil_N))/new_area;
       (currents->new_patch[target_lu])->structural_soil_L = (old_structural_soil_L*old_area+(currents->new_patch[target_lu]->structural_soil_L-old_structural_soil_L))/new_area;
       (currents->new_patch[target_lu])->water = (old_water*old_area + (currents->new_patch[target_lu]->water-old_water))/new_area;
       (currents->new_patch[target_lu])->theta = (old_theta*old_area+(currents->new_patch[target_lu]->theta-old_theta))/new_area;
       (currents->new_patch[target_lu])->rh = (old_rh*old_area+(currents->new_patch[target_lu]->rh-old_rh))/new_area;
       (currents->new_patch[target_lu])->gpp_avg = (old_gpp_avg*old_area+(currents->new_patch[target_lu]->gpp_avg-old_gpp_avg))/new_area;
       (currents->new_patch[target_lu])->npp_avg = (old_npp_avg*old_area+(currents->new_patch[target_lu]->npp_avg-old_npp_avg))/new_area;
       (currents->new_patch[target_lu])->rh_avg = (old_rh_avg*old_area+(currents->new_patch[target_lu]->rh_avg-old_rh_avg))/new_area;
       (currents->new_patch[target_lu])->fire_emission = (old_fire_emission*old_area+(currents->new_patch[target_lu]->fire_emission-old_fire_emission))/new_area;
       (currents->new_patch[target_lu])->product_emission = (old_product_emission*old_area+(currents->new_patch[target_lu]->product_emission-old_product_emission))/new_area;
       (currents->new_patch[target_lu])->forest_harvested_c = (old_forest_harvested_c*old_area+(currents->new_patch[target_lu]->forest_harvested_c-old_forest_harvested_c))/new_area;
       
       (currents->new_patch[target_lu])->yr1_decay_product_pool = (old_yr1_decay_product_pool*old_area+(currents->new_patch[target_lu]->yr1_decay_product_pool-old_yr1_decay_product_pool))/new_area;
       (currents->new_patch[target_lu])->yr10_decay_product_pool = (old_yr10_decay_product_pool*old_area+(currents->new_patch[target_lu]->yr10_decay_product_pool-old_yr10_decay_product_pool))/new_area;
       (currents->new_patch[target_lu])->yr100_decay_product_pool = (old_yr100_decay_product_pool*old_area+(currents->new_patch[target_lu]->yr100_decay_product_pool-old_yr100_decay_product_pool))/new_area;
       
       (currents->new_patch[target_lu])->past_harvested_c = (old_past_harvested_c*old_area+(currents->new_patch[target_lu]->past_harvested_c-old_past_harvested_c))/new_area;
       (currents->new_patch[target_lu])->crop_harvested_c = (old_crop_harvested_c*old_area+(currents->new_patch[target_lu]->crop_harvested_c-old_crop_harvested_c))/new_area;
       
#if CHECK_C_CONSERVE
       ///CarbonConserve
       currentp = currents->youngest_patch[donor_lu];
       mlcc = NULL;
       while(currentp!=NULL)
       {
           mlcc = currentp->shortest;
           double tmp_tb = 0.0;
           while (mlcc!=NULL) {
               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
               mlcc = mlcc->taller;
           }
           all_tb_after += tmp_tb*currentp->area/data->area;
           all_sc_after += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
           all_f_harv_after += currentp->forest_harvested_c * currentp->area/data->area;
           all_p_harv_after += currentp->past_harvested_c * currentp->area/data->area;
           all_c_harv_after += currentp->crop_harvested_c * currentp->area/data->area;
           currentp = currentp->older;
       }
       tmp_sc3 = all_sc_after;
       currentp = currents->new_patch[target_lu];
       tmp_tb = 0.0;
       mlcc = currentp->shortest;
       while (mlcc!=NULL) {
           tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
           mlcc = mlcc->taller;
       }
       all_tb_after += tmp_tb*currentp->area/data->area;
       all_sc_after += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
       all_f_harv_after += currentp->forest_harvested_c * currentp->area/data->area;
       all_p_harv_after += currentp->past_harvested_c * currentp->area/data->area;
       all_c_harv_after += currentp->crop_harvested_c * currentp->area/data->area;
       all_tc_after = all_tb_after + all_sc_after + (all_f_harv_after+all_p_harv_after+all_c_harv_after)*data->deltat;
       
       if (abs(all_tc_after-all_tc_before)>1e-9)
       {
           printf("Carbon leakage in landuse_transition  : imbalance    %.15f site_tc_af %.15f site_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
           printf("                                      : site_sc_bf  %.15f site_tb_bf %.15f site_f_hrC_bf %.15f site_p_hrC_bf %.15f site_c_hrC_bf %.15f\n",all_sc_before,all_tb_before,all_f_harv_before,all_p_harv_before,all_c_harv_before);
           printf("                                      : site_sc_af  %.15f site_tb_af %.15f site_f_hrC_af %.15f site_p_hrC_af %.15f site_c_hrC_af %.15f\n",all_sc_after,all_tb_after,all_f_harv_after,all_p_harv_after,all_c_harv_after);
           printf(" --------------------------------------------------------------------------------------\n");
       }
#endif
       

       
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
       
#if CHECK_C_CONSERVE
       ///CarbonConserve
       double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0, all_forest_harC_before = 0.0;
       double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0, all_forest_harC_after = 0.0;
       currentp = currents->youngest_patch[donor_lu];
       cohort* mlcc = NULL;
       while(currentp!=NULL)
       {
           mlcc = currentp->shortest;
           double tmp_tb = 0.0;
           while (mlcc!=NULL) {
               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
               mlcc = mlcc->taller;
           }
           all_tb_before += tmp_tb*currentp->area/data->area;
           all_sc_before += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
           all_forest_harC_before += currentp->forest_harvested_c * currentp->area/data->area;
           currentp = currentp->older;
       }
       currentp = currents->new_patch[target_lu];
       double tmp_tb = 0.0;
       mlcc = currentp->shortest;
       while (mlcc!=NULL) {
           tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
           mlcc = mlcc->taller;
       }
       all_tb_before += tmp_tb*currentp->area/data->area;
       all_sc_before += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
       all_forest_harC_before += currentp->forest_harvested_c * currentp->area/data->area;
       all_tc_before = all_tb_before + all_sc_before + all_forest_harC_before;
#endif

       ///CarbonConserve
       double old_fast_soil_C = (currents->new_patch[target_lu])->fast_soil_C;
       double old_structural_soil_C = (currents->new_patch[target_lu])->structural_soil_C;
       double old_fast_soil_N = (currents->new_patch[target_lu])->fast_soil_N;
       double old_slow_soil_C = (currents->new_patch[target_lu])->slow_soil_C;
       double old_passive_soil_C = (currents->new_patch[target_lu])->passive_soil_C;
       double old_mineralized_soil_N = (currents->new_patch[target_lu])->mineralized_soil_N;
       double old_structural_soil_L = (currents->new_patch[target_lu])->structural_soil_L;
       double old_water = (currents->new_patch[target_lu])->water;
       double old_theta = (currents->new_patch[target_lu])->theta;
       double old_rh = (currents->new_patch[target_lu])->rh;
       double old_area = (currents->new_patch[target_lu])->area;
       double old_gpp_avg = (currents->new_patch[target_lu])->gpp_avg;
       double old_npp_avg = (currents->new_patch[target_lu])->npp_avg;
       double old_rh_avg = (currents->new_patch[target_lu])->rh_avg;
       double old_fire_emission = (currents->new_patch[target_lu])->fire_emission;
       double old_product_emission = (currents->new_patch[target_lu])->product_emission;
       double old_forest_harvested_c = (currents->new_patch[target_lu])->forest_harvested_c;
       double old_yr1_decay_product_pool = (currents->new_patch[target_lu])->yr1_decay_product_pool;
       double old_yr10_decay_product_pool = (currents->new_patch[target_lu])->yr10_decay_product_pool;
       double old_yr100_decay_product_pool = (currents->new_patch[target_lu])->yr100_decay_product_pool;
       double old_past_harvested_c = (currents->new_patch[target_lu])->past_harvested_c;
       double old_crop_harvested_c = (currents->new_patch[target_lu])->crop_harvested_c;
       
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
               terminate_patch(&currentp,data);
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
//      (currents->new_patch[target_lu])->fast_soil_C /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->structural_soil_C /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->fast_soil_N /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->slow_soil_C /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->passive_soil_C /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->mineralized_soil_N /=(currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->structural_soil_L /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->water /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->theta /= (currents->new_patch[target_lu])->area;
//      (currents->new_patch[target_lu])->rh /= (currents->new_patch[target_lu])->area;
       
       
       ///CarbonConserve
       ///CarbonConserve
       /// The above code for averaging the lastes soil state is problematic,
       /// Because before accumualting new soil state from currentp, the unit of soil state in new_patch is density-based,
       /// but aggregate_in_soil_state_from_disturbance function will accumulate absolute state (not density but total amount) into newpatch, there to convert the latest soil state variable
       /// into density-base unit again, we should use the following code rather than the above -- Lei
       double new_area =currents->new_patch[target_lu]->area;
       (currents->new_patch[target_lu])->fast_soil_C = (old_fast_soil_C*old_area+(currents->new_patch[target_lu]->fast_soil_C-old_fast_soil_C))/new_area;
       (currents->new_patch[target_lu])->structural_soil_C = (old_structural_soil_C*old_area+(currents->new_patch[target_lu]->structural_soil_C-old_structural_soil_C))/new_area;
       (currents->new_patch[target_lu])->fast_soil_N = (old_fast_soil_N*old_area+(currents->new_patch[target_lu]->fast_soil_N-old_fast_soil_N))/new_area;
       (currents->new_patch[target_lu])->slow_soil_C = (old_slow_soil_C*old_area+(currents->new_patch[target_lu]->slow_soil_C-old_slow_soil_C))/new_area;
       (currents->new_patch[target_lu])->passive_soil_C = (old_passive_soil_C*old_area+(currents->new_patch[target_lu]->passive_soil_C-old_passive_soil_C))/new_area;
       (currents->new_patch[target_lu])->mineralized_soil_N = (old_mineralized_soil_N*old_area+(currents->new_patch[target_lu]->mineralized_soil_N-old_mineralized_soil_N))/new_area;
       (currents->new_patch[target_lu])->structural_soil_L = (old_structural_soil_L*old_area+(currents->new_patch[target_lu]->structural_soil_L-old_structural_soil_L))/new_area;
       (currents->new_patch[target_lu])->water = (old_water*old_area + (currents->new_patch[target_lu]->water-old_water))/new_area;
       (currents->new_patch[target_lu])->theta = (old_theta*old_area+(currents->new_patch[target_lu]->theta-old_theta))/new_area;
       (currents->new_patch[target_lu])->rh = (old_rh*old_area+(currents->new_patch[target_lu]->rh-old_rh))/new_area;
       (currents->new_patch[target_lu])->gpp_avg = (old_gpp_avg*old_area+(currents->new_patch[target_lu]->gpp_avg-old_gpp_avg))/new_area;
       (currents->new_patch[target_lu])->npp_avg = (old_npp_avg*old_area+(currents->new_patch[target_lu]->npp_avg-old_npp_avg))/new_area;
       (currents->new_patch[target_lu])->rh_avg = (old_rh_avg*old_area+(currents->new_patch[target_lu]->rh_avg-old_rh_avg))/new_area;
       (currents->new_patch[target_lu])->fire_emission = (old_fire_emission*old_area+(currents->new_patch[target_lu]->fire_emission-old_fire_emission))/new_area;
       (currents->new_patch[target_lu])->product_emission = (old_product_emission*old_area+(currents->new_patch[target_lu]->product_emission-old_product_emission))/new_area;
       (currents->new_patch[target_lu])->forest_harvested_c = (old_forest_harvested_c*old_area+(currents->new_patch[target_lu]->forest_harvested_c-old_forest_harvested_c))/new_area;
       (currents->new_patch[target_lu])->yr1_decay_product_pool = (old_yr1_decay_product_pool*old_area+(currents->new_patch[target_lu]->yr1_decay_product_pool-old_yr1_decay_product_pool))/new_area;
       (currents->new_patch[target_lu])->yr10_decay_product_pool = (old_yr10_decay_product_pool*old_area+(currents->new_patch[target_lu]->yr10_decay_product_pool-old_yr10_decay_product_pool))/new_area;
       (currents->new_patch[target_lu])->yr100_decay_product_pool = (old_yr100_decay_product_pool*old_area+(currents->new_patch[target_lu]->yr100_decay_product_pool-old_yr100_decay_product_pool))/new_area;
       (currents->new_patch[target_lu])->past_harvested_c = (old_past_harvested_c*old_area+(currents->new_patch[target_lu]->past_harvested_c-old_past_harvested_c))/new_area;
       (currents->new_patch[target_lu])->crop_harvested_c = (old_crop_harvested_c*old_area+(currents->new_patch[target_lu]->crop_harvested_c-old_crop_harvested_c))/new_area;
       
//       (currents->new_patch[target_lu])->gpp_avg /= (currents->new_patch[target_lu])->area;
//       (currents->new_patch[target_lu])->npp_avg /= (currents->new_patch[target_lu])->area;
//       (currents->new_patch[target_lu])->rh_avg /= (currents->new_patch[target_lu])->area;
//       (currents->new_patch[target_lu])->fire_emission /=(currents->new_patch[target_lu])->area;
//       (currents->new_patch[target_lu])->forest_harvested_c /= (currents->new_patch[target_lu])->area;
       
#if CHECK_C_CONSERVE
       ///CarbonConserve
       currentp = currents->youngest_patch[donor_lu];
       mlcc = NULL;
       while(currentp!=NULL)
       {
           mlcc = currentp->shortest;
           double tmp_tb = 0.0;
           while (mlcc!=NULL) {
               tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
               mlcc = mlcc->taller;
           }
           all_tb_after += tmp_tb*currentp->area/data->area;
           all_sc_after += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
           all_forest_harC_after += currentp->forest_harvested_c * currentp->area/data->area;
           currentp = currentp->older;
       }
       currentp = currents->new_patch[target_lu];
       tmp_tb = 0.0;
       mlcc = currentp->shortest;
       while (mlcc!=NULL) {
           tmp_tb += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
           mlcc = mlcc->taller;
       }
       all_tb_after += tmp_tb*currentp->area/data->area;
       all_sc_after += (currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C)*currentp->area/data->area;
       all_forest_harC_after += currentp->forest_harvested_c * currentp->area/data->area;
       all_tc_after = all_tb_after + all_sc_after + (all_forest_harC_after-all_forest_harC_before)*data->deltat;
       
       if (abs(all_tc_after-all_tc_before)>1e-9)
       {
           printf("Carbon leakage in cut_forest  : imbalance    %.15f site_tc_af %.15f site_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
           printf("                              : site_sc_bf  %.15f site_tb_bf %.15f site_f_hrC_bf %.15f\n",all_sc_before,all_tb_before,all_forest_harC_before);
           printf("                              : site_sc_af  %.15f site_tb_af %.15f site_f_hrC_af %.15f\n",all_sc_after,all_tb_after,all_forest_harC_after);
           printf("                              : site_lat %.3f site_lon %.3f\n",currents->sdata->lat_,currents->sdata->lon_);
           printf(" --------------------------------------------------------------------------------------\n");
       }
#endif
       
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
       
#if CHECK_C_CONSERVE
       cohort* mlcc = NULL;
       double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0;
       double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0, all_crop_harvest_af = 0.0;
       double actual_dt_tc = 0.0, esti_dt_tc = 0.0;
       mlcc = currentp->shortest;
       while (mlcc!=NULL) {
           all_tb_before += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
           mlcc = mlcc->taller;
       }
       all_sc_before = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
       all_tc_before = all_tb_before + all_sc_before;
#endif

//      cohort* currentc = currentp->tallest;
//      while (currentc != NULL) {
//         /*harvest everything*/
//         harvested_structural_C = currentc->bdead * currentc->nindivs;
//         harvested_fast_C = currentc->balive * currentc->nindivs;
//         harvested_fast_N = (1.0 / data->c2n_leaf[currentc->species])
//            * currentc->balive * currentc->nindivs;
//         currentc->nindivs = 0.00001;
//
//         currentc = currentc->shorter;
//      }

       ///CarbonConserve
       /// The above commented block is problematic as it does not accumulate the carbon loss from all cohorts but just use the attribute of last cohort
       /// It is corrected in the following
       cohort* currentc = currentp->tallest;
       while (currentc != NULL) {
           /*harvest everything*/
           harvested_structural_C += currentc->bdead * (currentc->nindivs-0.00001);
           harvested_fast_C += currentc->balive * (currentc->nindivs-0.00001);
           harvested_fast_N += (1.0 / data->c2n_leaf[currentc->species])
           * currentc->balive * (currentc->nindivs-0.00001);
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
       
       ///CarbonConserve
       ///Collect additional carbon loss to crop_harvested_c
       currentp->crop_harvested_c += N_CLIMATE*(1-data->crop_residue)*(harvested_structural_C+harvested_fast_C)/currentp->area;
       
      /*plant new crop*/
      //init_cohorts(&currentp, data);
       //Terminate all cohorts and keep the land clear until planting
       terminate_cohorts(&currentp->tallest,&currentp->shortest, data);
#elif defined MIAMI_LU
      harvested_structural_C = currentp->total_biomass;
      /*currentp->fast_soil_C += data->fraction_balive_2_fast 
          * data->crop_residue * harvested_fast_C; */
      currentp->structural_soil_C += data->crop_residue * harvested_structural_C 
         + (1.0 - data->fraction_balive_2_fast) * data->crop_residue * harvested_fast_C;    

      /*resest total biomass*/
      currentp->total_biomass = 0.0;
      currentp->total_ag_biomass = 0.0;
       update_patch (&currentp,data);  //Update patch as all cohorts were already terminated
#endif /* ED V. MIAMI_LU */
       
#if CHECK_C_CONSERVE
       ///CarbonConserve
       mlcc = currentp->shortest;
       while (mlcc!=NULL) {
           all_tb_after += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
           mlcc = mlcc->taller;
       }
       all_sc_after = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
       all_crop_harvest_af = currentp->crop_harvested_c;
       all_tc_after = all_tb_after + all_sc_after + all_crop_harvest_af * data->deltat;
       
       if (abs(all_tc_after-all_tc_before)>1e-9)
       {
           printf("Carbon leakage in harvest_croplands : imbalance   %.15f patch_tc_af %.15f patch_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
           printf("                                    : patch_sc_bf %.15f patch_tb_bf %.15f \n",all_sc_before,all_tb_before);
           printf("                                    : patch_sc_af %.15f patch_tb_af %.15f all_crop_harvest_af %.15f\n",all_sc_after,all_tb_after,all_crop_harvest_af);
           printf(" --------------------------------------------------------------------------------------\n");
       }
#endif

      currentp = currentp->older;
   }
}
////////////////////////////////////////////////////////////////////////////////
//! plant_croplands
//!
//!
//! @param
//! @return
////////////////////////////////////////////////////////////////////////////////
void plant_croplands(site** siteptr, UserData* data)
{
    //plant the crops by initilizing all cohorts if crop calendar indicates
    site* currents = *siteptr;
    patch* currentp = currents->youngest_patch[LU_CROP];
    ///CarbonConserve
    /// As to initiliza cohorts in crop patches will cause carbon leakage, there the added total biomass will be deducted from harvested biomass to compensate the carbon leakage -- Lei
    double total_tb_beforePlanted = 0.0, total_tb_afterPlanted = 0.0;
    cohort* currentc = NULL;
    while (currentp != NULL)
    {
        /*resest total biomass*/
        currentp->total_biomass = 0.0;
        currentp->total_ag_biomass = 0.0;
        
        ///CarbonConserve
        /// Here, calculate total biomass before planting new cohorts
        currentc = currentp->shortest;
        while(currentc!=NULL)
        {
            total_tb_beforePlanted += (currentc->balive+currentc->bdead)*currentc->nindivs/currentp->area;
            currentc = currentc->taller;
        }
        
        init_cohorts(&currentp, data);
        
        /// Here, calculate total biomass after planting new cohorts, the difference is the biomass of planted cohort and it will be deducted from harvested crop biomass -- Lei
        currentc = currentp->shortest;
        while(currentc!=NULL)
        {
            total_tb_afterPlanted += (currentc->balive+currentc->bdead)*currentc->nindivs/currentp->area;
            currentc = currentc->taller;
        }
        currentp->crop_harvested_c -= (total_tb_afterPlanted-total_tb_beforePlanted)*LANDUSE_FREQ;
        
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
       
#if CHECK_C_CONSERVE
       cohort* mlcc = NULL;
       double all_tb_before = 0.0, all_sc_before =0.0, all_tc_before =0.0;
       double all_tb_after = 0.0, all_sc_after =0.0, all_tc_after =0.0, all_crop_harvest_af = 0.0;
       double actual_dt_tc = 0.0, esti_dt_tc = 0.0;
       mlcc = currentp->shortest;
       while (mlcc!=NULL) {
           all_tb_before += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
           mlcc = mlcc->taller;
       }
       all_sc_before = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
       all_tc_before = all_tb_before + all_sc_before;
#endif

//      cohort* currentc = currentp->tallest;
//      while (currentc != NULL) {
//         /* 1. graze grasses */
//         if (data->is_grass[currentc->species]) {
//            grazed_structural_C = data->grazing_intensity * currentc->bdead * currentc->nindivs;
//            grazed_fast_C = data->grazing_intensity * currentc->balive * currentc->nindivs;
//            grazed_fast_N = (1.0 / data->c2n_leaf[currentc->species])
//               * data->grazing_intensity * currentc->balive * currentc->nindivs;
//            currentc->nindivs *= (1.0 - data->grazing_intensity);
//         }
//         /* 2. remove trees explicitly if primary vegetation forested*/
//         if ( (currents->maintain_pasture_flag == 1) && (!data->is_grass[currentc->species]) ) {
//            grazed_structural_C = currentc->bdead * currentc->nindivs;
//            grazed_fast_C = currentc->balive * currentc->nindivs;
//            grazed_fast_N = (1.0 / data->c2n_leaf[currentc->species])
//               * currentc->balive * currentc->nindivs;
//            currentc->nindivs = 0.000001;
//         }
//         currentc = currentc->shorter;
//      }
       
       ///CarbonConserve
       /// The above commmented block is from ED and it is problematic as grazed_structural_c should use += than =
       /// The below block corrent this problem -- Lei
       cohort* currentc = currentp->tallest;
       while (currentc != NULL) {
           /* 1. graze grasses */
           if (data->is_grass[currentc->species]) {
               grazed_structural_C += data->grazing_intensity * currentc->bdead * currentc->nindivs;
               grazed_fast_C += data->grazing_intensity * currentc->balive * currentc->nindivs;
               grazed_fast_N += (1.0 / data->c2n_leaf[currentc->species])
               * data->grazing_intensity * currentc->balive * currentc->nindivs;
               currentc->nindivs *= (1.0 - data->grazing_intensity);
           }
           /* 2. remove trees explicitly if primary vegetation forested*/
           if ( (currents->maintain_pasture_flag == 1) && (!data->is_grass[currentc->species]) ) {
               grazed_structural_C += currentc->bdead * (currentc->nindivs-0.000001);
               grazed_fast_C += currentc->balive * (currentc->nindivs-0.000001);
               grazed_fast_N += (1.0 / data->c2n_leaf[currentc->species])
               * currentc->balive * (currentc->nindivs-0.000001);
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
       
       ///CarbonConserve
       ///Collect carbon loss back to past_harvested_c
       currentp->past_harvested_c += N_CLIMATE*(1.0-data->grazing_residue) * (grazed_fast_C+grazed_structural_C)/currentp->area;
       
#if CHECK_C_CONSERVE
       ///CarbonConserve
       mlcc = currentp->shortest;
       while (mlcc!=NULL) {
           all_tb_after += (mlcc->balive+mlcc->bdead)*mlcc->nindivs/currentp->area;
           mlcc = mlcc->taller;
       }
       all_sc_after = currentp->fast_soil_C+currentp->slow_soil_C+currentp->structural_soil_C+currentp->passive_soil_C;
       all_crop_harvest_af = currentp->past_harvested_c;
       all_tc_after = all_tb_after + all_sc_after + all_crop_harvest_af*data->deltat;
       
       if (abs(all_tc_after-all_tc_before)>1e-9)
       {
           printf("Carbon leakage in graze_pasture  : imbalance    %.15f patch_tc_af %.15f patch_tc_bf %.15f\n",all_tc_after-all_tc_before,all_tc_after,all_tc_before);
           printf("                                 : patch_sc_bf  %.15f patch_tb_bf %.15f \n",all_sc_before,all_tb_before);
           printf("                                 : patch_sc_af  %.15f patch_tb_af %.15f all_crop_harvest_af %.15f\n",all_sc_after,all_tb_after,all_crop_harvest_af);
           printf("                                 : graed_fast_C %.15f grae_stru_C %.15f \n",grazed_fast_C,grazed_structural_C);
           printf(" --------------------------------------------------------------------------------------\n");
       }
#endif


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
       
////////////////////////////////////////////////////////////////////////////////
//! wood_pool_decay
//!
//!
//! @param
//! @return
////////////////////////////////////////////////////////////////////////////////
void wood_pool_decay(site** siteptr, UserData* data)
{
    // This function decay wood product pool and added this to prodcut emission
    site* currents = *siteptr;
    double decay = 0.0;
    for (int lu=0; lu<N_LANDUSE_TYPES; lu++)
    {
        patch* currentp = currents->youngest_patch[lu];
        while(currentp!=NULL)
        {
            currentp->product_emission = 0.0;
            
            decay = currentp->yr1_decay_product_pool * (1.0-exp(LANDUSE_FREQ*TIMESTEP*data->yr1_decay_rate));
            currentp->product_emission += decay;
            currentp->yr1_decay_product_pool -= decay;
            
            decay = currentp->yr10_decay_product_pool * (1.0-exp(LANDUSE_FREQ*TIMESTEP*data->yr10_decay_rate));
            currentp->product_emission += decay;
            currentp->yr10_decay_product_pool -= decay;
            
            decay = currentp->yr100_decay_product_pool * (1.0-exp(LANDUSE_FREQ*TIMESTEP*data->yr100_decay_rate));
            currentp->product_emission += decay;
            currentp->yr100_decay_product_pool -= decay;
            
            currentp = currentp->older;
        }
    }
}

/**********************************************************************/
#endif /* LANDUSE */
