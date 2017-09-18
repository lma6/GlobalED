#include <cmath>
#include <cstring>

#include "edmodels.h"
#include "cohort.h"

/* ALLOMETRY IN UNITS OF CARBON */
/* 10th Dec 98 structural calc from jgs agb allometry */
/* 14th Dec 98 constrined to maximum h *equal to repro hgt */

////////////////////////////////////////////////////////////////////////////////
//! Dbh
//! height(m) diameter(cm) relationships
//! O'Brien et al  - for 56 species at BCI
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double cohort::Dbh (UserData *data) { 
   double dbh;
   double m= 0.64;
   double c= 0.37;


   if ( (!strcmp(data->title[species],"evergreen"))  and data->allometry_type == 0) { /* spruce allometry */    
      dbh = exp( ( log(hite) - 0.04 ) / 0.94 );
   }
   else if(data->is_tropical[species] or data->is_grass[species] or 
           (!strcmp(data->title[species],"cold_decid") and 
           data->allometry_type == 0)) { /* all other species allometry */
      dbh = pow( 10.0, ( ( log10(hite) - c ) / m ) );
   }
   else {
          dbh = log( 1- (hite - data->ref_hgt[species])/data->b1Ht[species])/data->b2Ht[species];         
   }

   return(dbh); 
}

////////////////////////////////////////////////////////////////////////////////
//! Hite
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double cohort::Hite (UserData *data ) {
   double h;
   double m= 0.64;
   double c= 0.37;

   if ((!strcmp(data->title[species],"evergreen")) and data->allometry_type == 0) { /* spruce allometry */
      if (dbh <= data->max_dbh[species] )
         /* canadian forest service report */
         h = exp( 0.94 * log(dbh) + 0.04); 
      else
         h = exp( 0.94 * log(data->max_dbh[species]) + 0.04 );
   } else if(data->is_tropical[species] or data->is_grass[species] or 
           (!strcmp(data->title[species],"cold_decid") 
           and data->allometry_type == 0)) { /* all other species allometry */
      if (dbh <= data->max_dbh[species]) 
         h = pow( 10.0, (log10(dbh) * m + c) );
      else 
         h = pow( 10.0, (log10(data->max_dbh[species]) * m + c) );
   }
   else {
         dbh   = (dbh < data->max_dbh[species])? dbh : data->max_dbh[species];
         h = data->ref_hgt[species] + data->b1Ht[species] * (1 - exp( data->b2Ht[species] * dbh));
   }

   return h; 
}

////////////////////////////////////////////////////////////////////////////////
//! Bleaf
//! calc leaf biomass from height dbh and wood density using allometry of
//! J.G. Saldarriaga et al 1988 - Rio Negro, Journal of Ecology vol 76 p938-958
//! Use sepatate eqns for spruce in NA based on Weight Tables for Maine
//! 
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double cohort::Bleaf (UserData *data ) { 

   double bleaf;
   double a1,c1,d1;
   double a2,b2,c2,d2;
   double f,g,dcrit,p,q,r;
   int spp;
 
   spp = species;

   if ((!strcmp(data->title[spp],"evergreen")) and data->allometry_type == 0) { /*spruce allometry*/
      if (dbh <= data->max_dbh[spp] )
         bleaf = (1.0 / data->c2b) * (1.0 / 2.2) * 
            exp(-0.7980554 + 2.138061 * log(dbh / 2.54 )) + 0.005;
      else 
         bleaf = (1.0 / data->c2b) * (1.0 / 2.2) * 
            exp(-0.7980554 + 2.138061 * log(data->max_dbh[spp] / 2.54)) + 0.005;
    
   } else if(data->is_tropical[species] or data->is_grass[species] or 
           (!strcmp(data->title[spp],"cold_decid") and 
           data->allometry_type == 0)){ /* all other species allometry */
      a1 = -1.981;
      c1 = -0.584;
      d1 = 0.55;
      dcrit = 100.0;
    
      a2 = -4.111;
      b2 = 0.605;
      c2 = 0.848;
      d2 = 0.438;
      f  = 0.64;
      g  = 0.37;
    
      p  = a1 + c1 * g * log(10.0) + d1 * log(data->rho[spp]);
      r  = ( (a2 - a1) + g * log(10.0) * (c2 - c1) + 
             log(data->rho[spp]) * (d2 - d1) ) * ( 1 / log(dcrit) );
      q = 2.0 * b2 + c2 * f + r;  
    
      if(dbh <= data->max_dbh[spp])
         bleaf = (1.0 / data->c2b) * 
            ( exp(p) * pow(dbh, q) + data->bl_min[spp] );
      else 
         bleaf = (1.0 / data->c2b) * 
            ( exp(p) * pow(data->max_dbh[spp], q) + data->bl_min[spp] );
      }
   else {
      dbh   = (dbh < data->max_dbh[species])? dbh : data->max_dbh[species];
      bleaf = (1.0/data->c2b)*data->b1Bl[species] * pow(dbh,data->b2Bl[species]);
   }

   return bleaf;
}

////////////////////////////////////////////////////////////////////////////////
//! Bdead
//! calc stem biomass from height(m) dbh(cm) and wood density(g/cm3) using allometry
//! of J.G. Saldarriaga et al 1988 - Rio Negro, Journal of Ecology vol 76 p938-958
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double cohort::Bdead (UserData *data ) {
   double a1,c1,d1;
   double a2,b2,c2,d2;
   double bdead;
   double dcrit,p,q,r;
   double f,g;
   int spp;

   spp = species;
   
   if ((!strcmp(data->title[spp],"evergreen")) and data->allometry_type == 0) { /* spruce allometry */
      bdead = ( 1.0 / 2.2 ) * ( 1.0 / data->c2b ) * 
         exp( 1.10651 + 2.298388 * log( dbh / 2.54 ) ); 
   } else if(data->is_tropical[species] or data->is_grass[species] or 
           (!strcmp(data->title[spp],"cold_decid") and 
           data->allometry_type == 0)){ /* all other species allometr */
      a1 = -1.981;
      c1 = 0.572;
      d1 = 0.931;
      dcrit = 100.0;
      a2 = -1.086;
      b2 = 0.876;
      c2 = 0.604;
      d2 = 0.871;

      f  = 0.64;
      g  = 0.37;
    
      if ( dbh > data->max_dbh[spp] ) {  
         p = a1 + c1 * log(hite) + d1 * log(data->rho[spp]);
         r = ( (a2 - a1) + (c2 - c1) * log(hite) + 
               log(data->rho[spp]) * (d2 - d1) ) * ( 1 / log(dcrit) );
         q = 2.0 * b2 + r;
      } else {
         p  = a1 + c1 * g * log(10.0) + d1 * log(data->rho[spp]);
         r  = ( (a2 - a1) + g * log(10.0) * (c2 - c1) + 
                log(data->rho[spp]) * (d2 - d1) ) * ( 1 / log(dcrit) );
         q = 2.0 * b2 + c2 * f + r;     
      }
         
      bdead = ( 1.0 / data->c2b ) * ( exp(p) * pow(dbh,q) ) + data->bs_min[spp];
   }
   else {
      dbh   = (dbh < data->max_dbh[species])? dbh : data->max_dbh[species];
      bdead =  (1.0/data->c2b)*data->b1Bs[species] * pow(dbh,data->b2Bs[species]);
   }

   return bdead;
}

////////////////////////////////////////////////////////////////////////////////
//! dHdBd
//! Convert changes in structural biomass bs to changes in height bs is stem biomass 
//! plus strucrual root biomass consistent with Bstem and h-dbh allometries
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double cohort::dHdBd (UserData *data ) {
   double a1,c1,d1;
   double a2,b2,c2,d2;
   double dhdbs;
   double dcrit,q,r;
   double f,g;
   int spp;
   double dhddbh;
   double ddbhdbs;

   spp = species;
   
   if ((!strcmp(data->title[spp],"evergreen")) and data->allometry_type == 0) { /* spruce allometry */
      dhdbs = exp( 0.40898 * log(bdead *2.2 * data->c2b) + 0.4639) * 
         0.40898 / bdead;
   }
   else if(data->is_tropical[species] or data->is_grass[species] or 
           (!strcmp(data->title[spp],"cold_decid") and 
           data->allometry_type == 0)) { /* all other species allometry */
    
      a1 = -1.981;
      c1 = 0.572;
      d1 = 0.931;
      dcrit = 100.0;
      a2 = -1.086;
      b2 = 0.876;
      c2 = 0.604;
      d2 = 0.871;
      f  = 0.64;
      g  = 0.37;
  
      if ( dbh > data->max_dbh[spp] ) {  
         dhdbs = 0.0;
      }
      else {
         r  = ( (a2 - a1) + g * log(10.0) * (c2 - c1) + log(data->rho[spp]) 
                * (d2 - d1) ) * ( 1 / log(dcrit) );
         q = 2.0 * b2 + c2 * f + r;
         dhdbs = ( hite / (bdead) ) * ( 1.0 / q ) * f * 
            log(10.0) * log10( exp(1.0) ); 
      }
   } else {
       if(fabs((data->b1Ht[species]+data->ref_hgt[species]) - hite) > 0.01){
          dhddbh = -data->b1Ht[species] * data->b2Ht[species] * exp(data->b2Ht[species] * 
                  pow((bdead*data->c2b/data->b1Bs[species]),1.0/data->b2Bs[species]));
          ddbhdbs = pow(data->c2b/data->b1Bs[species],1.0/data->b2Bs[species])* 
                  pow(bdead, 1.0/data->b2Bs[species] - 1)*(1/data->b2Bs[species]);

          dhdbs = fabs(dhddbh * ddbhdbs); /* prevent negatives by rounding */
       }else{
          dhdbs = 0.0;                /* return 0 when height is within 1 cm of asymptote */
       }
   }

   return dhdbs;
}

////////////////////////////////////////////////////////////////////////////////
//! dDbhdBd
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double cohort::dDbhdBd (UserData *data ) {
   double ddbhdbs;

   double a1,c1,d1;
   double a2,b2,c2,d2;
   double f,g;
   double dcrit,q,r;
   int spp;

   spp = species;
   
   if ((!strcmp(data->title[spp],"evergreen")) and data->allometry_type == 0) { /* spruce allometry */
      ddbhdbs = 2.54 * exp( (log(bdead * 2.2 * data->c2b) - 1.10651 ) / 
                            2.298388 ) / (2.298388 * bdead);
   }
   else if(data->is_tropical[species] or data->is_grass[species] or 
           (!strcmp(data->title[spp],"cold_decid") and 
           data->allometry_type == 0)) { /* all other species allometry */

      a1 = -1.981;
      c1 = 0.572;
      d1 = 0.931;
      dcrit = 100.0;
      a2 = -1.086;
      b2 = 0.876;
      c2 = 0.604;
      d2 = 0.871;

      f  = 0.64;
      g  = 0.37;

      if (dbh > data->max_dbh[spp]) {  
         r  = ( (a2 - a1) + (c2 - c1) * log(hite) 
                + log(data->rho[spp]) * (d2 - d1) ) * ( 1 / log(dcrit) );
         q = 2.0 * b2 + r;
      }
      else {
         r  = ( (a2 - a1) + g * log(10.0) * (c2 - c1) + 
                log(data->rho[spp]) * (d2 - d1) ) * ( 1 / log(dcrit) );
         q = 2.0 * b2 + c2 * f + r;      
      }
  
      ddbhdbs = ( dbh / (bdead) ) * ( 1.0 / q );        
   }
   else {
        ddbhdbs = pow(data->c2b/data->b1Bs[species],(1.0/data->b2Bs[species]))*pow(bdead, (1.0/data->b2Bs[species] - 1))*(1/data->b2Bs[species]);       
   }
   
   return ddbhdbs;
}

////////////////////////////////////////////////////////////////////////////////
//! dDbhdBl
//! 
//!
//! @param  
//! @return 
////////////////////////////////////////////////////////////////////////////////
double cohort::dDbhdBl ( UserData *data ) {
   double ddbhdbl;
   double a1,c1,d1;
   double a2,b2,c2,d2;
   double dcrit,q,r;
   double f,g;
   int spp;

   spp = species;
 
   if ((!strcmp(data->title[spp],"evergreen")) and data->allometry_type == 0) { /* spruce allometry */
      ddbhdbl = 3.3763239 * pow(bl, -0.5322767);
   }
   else if(data->is_tropical[species] or data->is_grass[species] or 
           (!strcmp(data->title[spp],"cold_decid") and 
           data->allometry_type == 0)) { /* all other species allometry */
    
      a1 = -1.981;
      c1 = -0.584;
      d1 = 0.55;
      dcrit = 100.0;
    
      a2 = -4.111;
      b2 = 0.605;
      c2 = 0.848;
      d2 = 0.438;
    
      f  = 0.64;
      g  = 0.37;
      r  = ( (a2 - a1) + g * log(10.0) * (c2 - c1) + 
             log(data->rho[spp]) * (d2 - d1) ) * ( 1 / log(dcrit) );
      q = 2 * b2 + c2 * f + r;
    
      ddbhdbl = ( dbh / (bl) ) * ( 1.0 / (q) );
      
   }
   else {
       ddbhdbl = pow(data->c2b/data->b1Bl[species],1.0/data->b2Bl[species])*pow(bl, 1.0/data->b2Bl[species] - 1)*(1/data->b2Bl[species]);
   }

   return ddbhdbl;
}

/******************************************************************************/
/********************************  END OF FILE ********************************/
/******************************************************************************/
