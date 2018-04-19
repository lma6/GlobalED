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

bool SiteData::compute_mech(int pt, int spp, double Vm0, int Vm0_bin, int time_period, int light_index, UserData* data)
{
    double farquhar_results[6];
    double shade=0;
    double Vcmax25=Vm0;
    double tmp=0,hum=0,swd=0,windspeed=0,CO2=0,Pa=101.3,Tg=25.0,Ts=0;
    
    shade=light_levels[spp][light_index];
    tf_air[spp][time_period]=0;
    tf_soil[spp][time_period]=0;
    An[spp][time_period][light_index]=0;
    E[spp][time_period][light_index]=0;
    Anb[spp][time_period][light_index]=0;
    Eb[spp][time_period][light_index]=0;

#if CPOUPLE_VcmaxDownreg
    double cumuLAI=log(shade)/(-1.0/(data->cohort_shading*data->L_extinct));   //Derive cumulative LAI that results in the light level
    double Kn = exp(0.00963*data->Vm0_max[spp]* - 2.43);        //Lloyd et al 2011
    cumuLAI/=2.0;    //Suggestion from get_cohort_vm0 function in cohort.cc that only half LAI contributes to shading. --Lei
    Vcmax25=data->Vm0_max[spp]*exp(-Kn*cumuLAI);   //Downregulation from original defined vcmax of each species
    if (shade<1e-5) Vcmax25=0;
    printf("Dwonregulate Vcmax cumuLAI %f lite %f Vcmax25 %f %f\n",cumuLAI,shade,data->Vm0_max[spp],Vcmax25);
#endif
    
    
    //// Test calculation with FLUXNET2015
//    printf("start testing fluxnet\n");
//    char dirpath[500],expname[500],IGBP[100],timestep[10];
//
//    strcpy(dirpath,"/gpfs/data1/hurttgp/gel1/leima/AssignTask/gED/Data/FLUXNET2015/Subset/");
//    FILE* siteinfo;
//    siteinfo=fopen("/gpfs/data1/hurttgp/gel1/leima/AssignTask/gED/Data/FLUXNET2015/SiteInfo.txt","r");
//    double lat=0,lon=0,unit_convertor=0;
//    //////strcpy(expname,"FLX_CA-NS4_FLUXNET2015_FULLSET_2002-2005_1-3_HH");
//    strcpy(timestep,"_HH");  //Timestep of files
//
//    if (!strcmp(timestep,"_HH")) unit_convertor=1e6;    //As unit is different for HH, DD, MM and YY
//    if (!strcmp(timestep,"_DD")) unit_convertor=12*3600*2.0;
//    if (!strcmp(timestep,"_MM")) unit_convertor=12*3600*2.0;
//    if (!strcmp(timestep,"_YY")) unit_convertor=12*3600*2.0;
//
//    if (unit_convertor<0.5)
//    {
//        printf("Error in unit_convertor %f\n",unit_convertor);
//        exit(0);
//    }
//
//    while (fscanf(siteinfo,"%s %s %f %f\n",&expname,&IGBP,&lat,&lon)!=EOF)
//    {
//        char infilename[500],outfilename[500],buf[500];
//
//        //strcpy(expname,"FLX_AR-SLu_FLUXNET2015_FULLSET_2009-2011_1-3");
//        printf("Processing %s\n",expname);
//
//        FILE *infile, *outfile;
//        strcpy(infilename,dirpath);
//        strcat(infilename,expname);
//        strcat(infilename,timestep);
//        strcat(infilename,".txt");
//
//        strcpy(outfilename,dirpath);
//        strcat(outfilename,expname);
//        strcat(outfilename,timestep);
//        strcat(outfilename,"_ED.txt");
//
//        infile = fopen(infilename,"r");
//        outfile = fopen(outfilename,"w");
//        if(infile==NULL)
//            printf("infile is null - %s not exist\n",infilename);
//        if(outfilename==NULL)
//            printf("outfile is null - %s not exist\n",outfilename);
//
//
//        char tair_FLUXs[100],swd_FLUXs[100],VPD_FLUXs[100],Pa_FLUXs[100],Ws_FLUXs[100],CO2_FLUXs[100],GPP1_FLUXs[100],GPP2_FLUXs[100],GPP3_FLUXs[100],GPP4_FLUXs[100];
//        double tair_FLUX=0,swd_FLUX=0,VPD_FLUX=0,Pa_FLUX=0,Ws_FLUX=0,CO2_FLUX=0,GPP1_FLUX=0,GPP2_FLUX=0,GPP3_FLUX=0,GPP4_FLUX=0,Vcmax25_FLUX=0;
//        double results_FLUX[6];
//        int mon=0;
//
//        fscanf(infile, "%s %s %s %s %s %s %s %s %s %s\n",&tair_FLUXs,&swd_FLUXs,&VPD_FLUXs,&Pa_FLUXs,&Ws_FLUXs,&CO2_FLUXs,&GPP1_FLUXs,&GPP2_FLUXs,&GPP3_FLUXs,&GPP4_FLUXs);
//        printf("Header of infile - %s %s %s %s %s %s %s %s %s %s\n",tair_FLUXs,swd_FLUXs,VPD_FLUXs,Pa_FLUXs,Ws_FLUXs,CO2_FLUXs,GPP1_FLUXs,GPP2_FLUXs,GPP3_FLUXs,GPP4_FLUXs);
//        fprintf(outfile, "mon TA_F SW_IN_F VPD_F PA_F WS_F CO2_F_MDS GPP_NT_CUT_REF GPP_NT_VUT_REF GPP_DT_CUT_REF GPP_DT_VUT_REF GPP_ED\n");
//
//        if(lat>23.5 or lat<-23.5)
//        {
//            if (!strcmp(IGBP,"ENF"))    Vcmax25_FLUX=30.0;
//            if (!strcmp(IGBP,"EBF"))    Vcmax25_FLUX=30.0;  //Using grass for EBF
//            if (!strcmp(IGBP,"MF"))    Vcmax25_FLUX=30.0;
//            if (!strcmp(IGBP,"DBF"))    Vcmax25_FLUX=25.0;
//            if (!strcmp(IGBP,"GRA"))    Vcmax25_FLUX=90.0;
//            if (!strcmp(IGBP,"CRO"))    Vcmax25_FLUX=90.0;  //Using grass for open cropland
//            if (!strcmp(IGBP,"OSH"))    Vcmax25_FLUX=90.0;  //Using grass for open shrubland
//            if (!strcmp(IGBP,"CSH"))    Vcmax25_FLUX=90.0;  //Using grass for closed shrubland
//            if (!strcmp(IGBP,"SAV"))    Vcmax25_FLUX=90.0;  //Using grass for closed shrubland
//            if (!strcmp(IGBP,"WSA"))    Vcmax25_FLUX=90.0;  //Using grass for closed shrubland
//            if (!strcmp(IGBP,"WET"))    Vcmax25_FLUX=90.0;  //Using grass for closed shrubland
//            if (!strcmp(IGBP,"SNO"))    Vcmax25_FLUX=90.0;  //Using grass for closed shrubland
//        }
//        else
//        {
//            if (!strcmp(IGBP,"ENF"))    Vcmax25_FLUX=30.0;
//            if (!strcmp(IGBP,"EBF"))    Vcmax25_FLUX=30.0;  //Using grass for EBF
//            if (!strcmp(IGBP,"MF"))     Vcmax25_FLUX=30.0;
//            if (!strcmp(IGBP,"DBF"))    Vcmax25_FLUX=30.0;
//            if (!strcmp(IGBP,"GRA"))    Vcmax25_FLUX=90.0;
//            if (!strcmp(IGBP,"CRO"))    Vcmax25_FLUX=90.0;  //Using grass for open cropland
//            if (!strcmp(IGBP,"OSH"))    Vcmax25_FLUX=90.0;  //Using grass for open shrubland
//            if (!strcmp(IGBP,"CSH"))    Vcmax25_FLUX=90.0;  //Using grass for closed shrubland
//            if (!strcmp(IGBP,"SAV"))    Vcmax25_FLUX=90.0;  //Using grass for closed shrubland
//            if (!strcmp(IGBP,"WSA"))    Vcmax25_FLUX=90.0;  //Using grass for closed shrubland
//            if (!strcmp(IGBP,"WET"))    Vcmax25_FLUX=90.0;  //Using grass for closed shrubland
//            if (!strcmp(IGBP,"SNO"))    Vcmax25_FLUX=90.0;  //Using grass for closed shrubland
//        }
//
//
//        while (fscanf(infile,"%s %s %s %s %s %s %s %s %s %s\n",&tair_FLUXs,&swd_FLUXs,&VPD_FLUXs,&Pa_FLUXs,&Ws_FLUXs,&CO2_FLUXs,&GPP1_FLUXs,&GPP2_FLUXs,&GPP3_FLUXs,&GPP4_FLUXs)!=EOF)
//        {
//            tair_FLUX=atof(tair_FLUXs);
//            swd_FLUX=atof(swd_FLUXs);
//            VPD_FLUX=atof(VPD_FLUXs);
//            Pa_FLUX=atof(Pa_FLUXs);
//            Ws_FLUX=atof(Ws_FLUXs);
//            CO2_FLUX=atof(CO2_FLUXs);
//            GPP1_FLUX=atof(GPP1_FLUXs);
//            GPP2_FLUX=atof(GPP2_FLUXs);
//            GPP3_FLUX=atof(GPP3_FLUXs);
//            GPP4_FLUX=atof(GPP4_FLUXs);
//            if (tair_FLUX<-100)
//            {
//                printf("Error in tair %f\n",tair_FLUX);
//                printf("ErrorLine in Tair tair %s swd %s VPD %s Pa %s WS %s CO2 %s GPP1 %s GPP2 %s GPP3 %s GPP4 %s\n",tair_FLUXs,swd_FLUXs,VPD_FLUXs,Pa_FLUXs,Ws_FLUXs,CO2_FLUXs,GPP1_FLUXs,GPP2_FLUXs,GPP3_FLUXs,GPP4_FLUXs);
//                exit(0);
//            }
//            if (swd_FLUX<-100)
//            {
//                printf("Error in swd_FLUX %f tair is %f\n", swd_FLUX,tair_FLUX);
//                printf("ErrorLine SWD tair %s swd %s VPD %s Pa %s WS %s CO2 %s GPP1 %s GPP2 %s GPP3 %s GPP4 %s\n",tair_FLUXs,swd_FLUXs,VPD_FLUXs,Pa_FLUXs,Ws_FLUXs,CO2_FLUXs,GPP1_FLUXs,GPP2_FLUXs,GPP3_FLUXs,GPP4_FLUXs);
//                exit(0);
//            }
//            if (VPD_FLUX<-100)
//            {
//                printf("Error in VPD_FLUX %f\n", VPD_FLUX);
//                printf("ErrorLine in VPD tair %s swd %s VPD %s Pa %s WS %s CO2 %s GPP1 %s GPP2 %s GPP3 %s GPP4 %s\n",tair_FLUXs,swd_FLUXs,VPD_FLUXs,Pa_FLUXs,Ws_FLUXs,CO2_FLUXs,GPP1_FLUXs,GPP2_FLUXs,GPP3_FLUXs,GPP4_FLUXs);
//                exit(0);
//
//            }
//            if (Pa_FLUX<-100)
//            {
//                /*printf("Error in Pa_FLUX %f and assign to 100.0\n", Pa_FLUX);*/
//                //printf("ErrorLine in Pa tair %s swd %s VPD %s Pa %s WS %s CO2 %s GPP1 %s GPP2 %s GPP3 %s GPP4 %s\n",tair_FLUXs,swd_FLUXs,VPD_FLUXs,Pa_FLUXs,Ws_FLUXs,CO2_FLUXs,GPP1_FLUXs,GPP2_FLUXs,GPP3_FLUXs,GPP4_FLUXs);
//                Pa_FLUX=101.3;
//            }
//            if (Ws_FLUX<-100)
//            {
//                printf("Error in Ws_FLUX %f in %d\n", Ws_FLUX,mon);
//                printf("ErrorLine in Ws tair %s swd %s VPD %s Pa %s WS %s CO2 %s GPP1 %s GPP2 %s GPP3 %s GPP4 %s\n",tair_FLUXs,swd_FLUXs,VPD_FLUXs,Pa_FLUXs,Ws_FLUXs,CO2_FLUXs,GPP1_FLUXs,GPP2_FLUXs,GPP3_FLUXs,GPP4_FLUXs);
//                exit(0);
//
//            }
//            if (CO2_FLUX<-100)
//            {
//                /*printf("Error in CO2_FLUX %f and assign to 350.0\n", CO2_FLUX);*/
//                //printf("ErrorLine in CO2 tair %s swd %s VPD %s Pa %s WS %s CO2 %s GPP1 %s GPP2 %s GPP3 %s GPP4 %s\n",tair_FLUXs,swd_FLUXs,VPD_FLUXs,Pa_FLUXs,Ws_FLUXs,CO2_FLUXs,GPP1_FLUXs,GPP2_FLUXs,GPP3_FLUXs,GPP4_FLUXs);
//                CO2_FLUX=400;
//            }
//            if (Vcmax25_FLUX<0.05)
//            {
//                printf("Error in Vcmax25_FLUX %f\n", Vcmax25_FLUX);
//                exit(0);
//            }
//
//            double Es=(0.611*exp(17.502*tair_FLUX/(240.97+tair_FLUX)))/1e2;
//            double Ea=Es-VPD_FLUX/1000.0;
//
//            Farquhar_couple(0,0,data,tair_FLUX,tair_FLUX,Ea,swd_FLUX,tair_FLUX,CO2_FLUX,Ws_FLUX,Pa_FLUX,1,Vcmax25_FLUX,results_FLUX);
//
////            printf("mon %d tair %f swd %f Ea %f Pa %f ws %f co2 %f gpp1 %f gpp2 %f gpp3 %f gpp4 %f gpp_ED %f Vc %f\n",mon,tair_FLUX,swd_FLUX,Ea/Es,Pa_FLUX,Ws_FLUX,CO2_FLUX,GPP1_FLUX,GPP2_FLUX,GPP3_FLUX,GPP4_FLUX,results_FLUX[1]*unit_convertor,Vcmax25_FLUX);
//
//            mon++;
//            fprintf(outfile, "%d %f %f %f %f %f %f %f %f %f %f %f \n",mon,tair_FLUX,swd_FLUX,VPD_FLUX,Pa_FLUX,Ws_FLUX,CO2_FLUX,GPP1_FLUX,GPP2_FLUX,GPP3_FLUX,GPP4_FLUX,results_FLUX[1]*unit_convertor);
//        }
//        fclose(infile);
//        fclose(outfile);
//        printf("Done\n");
//        //
//
//    }
//    exit(0);

    
//    double temp1=0,hum1=0,swd1=0,co21=0,ws1=0,shade1=1;
//    int lat_i=186,lon_i=249;
//    Vcmax25=15;
//    double farquhar_results1[6],farquhar_results2[6];
//    for(size_t mon=0;mon<12;mon++)
//    {
//        for (size_t hr=mon*24;hr<mon*24+24;hr++)
//        {
//            temp1=data->global_tmp[hr][lat_i][lon_i];
//            hum1=data->global_hum[hr][lat_i][lon_i];
//            swd1=data->global_swd[hr][lat_i][lon_i];
//            ws1=data->global_windspeed[hr][lat_i][lon_i];
//            co21=data->global_CO2[hr][lat_i][lon_i];
//            Farquhar_couple(0,4,data,temp1,temp1,hum1,swd1,temp1,co21*350/390,ws1,Pa,1,Vcmax25,farquhar_results1);
//            //Farquhar_couple(pt,4,data,temp1,temp1,hum1,swd1,temp1,350,ws1,Pa,1,15,farquhar_results2);
//            //farquhar(Vcmax25/1e6,co21/1e6,temp1,temp1,hum1,swd1,shade1,0,farquhar_results1);
//            //farquhar(Vcmax25/1e6,396.65/1e6,temp1,temp1,hum1,swd1,shade1,0,farquhar_results2);
//            printf("mon %d hr %d tair %f hum %f swd %f ws %f co2 %f An1 %f An2 %f\n",mon,hr,temp1,hum1,swd1,ws1,co21,farquhar_results1[1]*1e6,farquhar_results2[1]*1e6);
//        }
//    }
//    exit(0);
    
    
//#if COUPLE_MERRA2_LUT
    if(data->MERRA2_LUT)
    {
    if (1-is_filled_LUT[data->MERRA2_timestamp-1][spp][time_period])  ////if MECH LUT has not been computed before
    {
        //Compute for all co2 and light combination
        double tmp_airtemp=0,tmp_soiltemp=0,tmp_hum=0,tmp_swd=0,tmp_windspeed=0,tmp_co2=0,tmp_shade=0;
        for (size_t lite_i=0;lite_i<N_bins_LUT_LITE;lite_i++)
        {
            double tmp_An=0,tmp_Anb=0,tmp_E=0,tmp_Eb=0,tmp_tf_air=0,tmp_tf_soil=0;
            for (size_t mon=time_period*24;mon<time_period*24+24;mon++)
            {
                tmp_airtemp=data->global_tmp[mon][globY_][globX_];
                tmp_hum=data->global_hum[mon][globY_][globX_];
                tmp_swd=data->global_swd[mon][globY_][globX_];
                tmp_windspeed=data->global_windspeed[mon][globY_][globX_];
                tmp_soiltemp=data->global_soiltmp[mon][globY_][globX_];
                if(data->mechanism_year<1850)
                    tmp_co2=280.0/390.0*data->global_CO2[mon][globY_][globX_];
                else if(data->mechanism_year>=1850 and data->mechanism_year<1950)
                    tmp_co2=(280.0+0.314*(data->mechanism_year-1850))/390.0*data->global_CO2[mon][globY_][globX_];
                else if(data->mechanism_year>=1950 and data->mechanism_year<2001)
                    tmp_co2=(311+1.290*(data->mechanism_year-1950))/390.0*data->global_CO2[mon][globY_][globX_];
                else
                    tmp_co2=data->global_CO2[mon][globY_][globX_];
                
                tmp_shade=light_levels[spp][LITE_IDX[lite_i]];
                
                double tmp_Tg=0;
                for (size_t mon1=time_period*24;mon1<time_period*24+24;mon1++)
                {
                    tmp_Tg+=data->global_tmp[mon1][globY_][globX_];
                }
                tmp_Tg/=24.0;
                
                Farquhar_couple(pt,spp,data,tmp_airtemp,tmp_soiltemp,tmp_hum,tmp_swd,tmp_Tg,tmp_co2,tmp_windspeed,Pa,tmp_shade,Vcmax25,farquhar_results);
                
                tmp_tf_air+=farquhar_results[0];
                tmp_tf_soil+=farquhar_results[5];
                tmp_An+=farquhar_results[1];
                tmp_E+=farquhar_results[2];
                tmp_Anb+=farquhar_results[3];
                tmp_Eb+=farquhar_results[4];
            }
            tf_LUT_air[data->MERRA2_timestamp-1][spp][time_period]=tmp_tf_air/24.0;
            tf_LUT_soil[data->MERRA2_timestamp-1][spp][time_period]=tmp_tf_soil/24.0;
            An_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i]=tmp_An*3600.0*360.0;
            E_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i]=tmp_E*3600.0*540.0;
            Anb_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i]=tmp_Anb*3600.0*360.0;
            Eb_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i]=tmp_Eb*3600.0*540.0;
        }
        //flag for this spp at MERRA2_timestamp and time_period, avoid repeating computation
        is_filled_LUT[data->MERRA2_timestamp-1][spp][time_period]=1;
        //printf("Finished computing of LUT lat %d lon %d mechyear %d timestamp %d mon %d spp %d co2 %f \n",globY_,globX_,data->mechanism_year,data->MERRA2_timestamp-1,time_period,spp,tmp_co2);
        //exit(0);
    }
    if (An_LUT[data->MERRA2_timestamp-1][spp][time_period][0]<-1000)
    {
        printf("Error in An_LUT initilization\n");
        exit(0);
    }
    int lite_i=0,lite_i_1=0; //co2_i_1 and co2_i are left and right points for linear interpolation, same to lite_i and lite_i_1
    while (LITE_IDX[lite_i]<=light_index and lite_i<N_bins_LUT_LITE)
    {
        lite_i++;
    }
    if (lite_i==0) lite_i_1=0;  else    lite_i_1=lite_i-1;
    double weight_lite=(LITE_IDX[lite_i]-light_index)/(LITE_IDX[lite_i]-LITE_IDX[lite_i-1]);
    An[spp][time_period][light_index]=weight_lite*An_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i_1]
                                    +(1-weight_lite)*An_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i];
    
    Anb[spp][time_period][light_index]=weight_lite*Anb_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i_1]
                                        +(1-weight_lite)*Anb_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i];
    
    E[spp][time_period][light_index]=weight_lite*E_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i_1]
                                    +(1-weight_lite)*E_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i];
    
    Eb[spp][time_period][light_index]=weight_lite*Eb_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i_1]
                                    +(1-weight_lite)*Eb_LUT[data->MERRA2_timestamp-1][spp][time_period][lite_i];
    
    tf_air[spp][time_period]=tf_LUT_air[data->MERRA2_timestamp-1][spp][time_period];
    tf_soil[spp][time_period]=tf_LUT_soil[data->MERRA2_timestamp-1][spp][time_period];

//#else  //COUPLE_MERRA2_LUT
    }
        else
        {
    for (size_t mon=time_period*24;mon<time_period*24+24;mon++)
    {
        tmp=data->global_tmp[mon][globY_][globX_];
        hum=data->global_hum[mon][globY_][globX_];
        swd=data->global_swd[mon][globY_][globX_];
        windspeed=data->global_windspeed[mon][globY_][globX_];
        Ts=data->global_soiltmp[mon][globY_][globX_];
        
        if(data->mechanism_year<1850)
            CO2=280.0/390.0*data->global_CO2[mon][globY_][globX_];
        else if(data->mechanism_year>=1850 and data->mechanism_year<1950)
            CO2=(280.0+0.314*(data->mechanism_year-1850))/390.0*data->global_CO2[mon][globY_][globX_];
        else if(data->mechanism_year>=1950 and data->mechanism_year<2001)
            CO2=(311+1.290*(data->mechanism_year-1950))/390.0*data->global_CO2[mon][globY_][globX_];
        else
            CO2=data->global_CO2[mon][globY_][globX_];

            
            //As growth temperature defined in Lombardozzi et all 2015 and Atkin et al 2008 is the preceding 10 days running mean of air temperature, here for simplicity, use mean temperature of current month
            Tg=0;
            for (size_t mon1=time_period*24;mon1<time_period*24+24;mon1++)
            {
                Tg+=data->global_tmp[mon1][globY_][globX_];
            }
            Tg/=24.0;
            
            ////Currently, ambient CO2 concentration is 350 umol
            //farquhar(Vcmax25/1e6,CA,data->global_tmp[mon][globY_][globX_],Ts,data->global_hum[mon][globY_][globX_],data->global_swd[mon][globY_][globX_],shade,pt,farquhar_results);
            
            Farquhar_couple(pt,spp,data,tmp,Ts,hum,swd,Tg,CO2,windspeed,Pa,shade,Vcmax25,farquhar_results);
            
            tf_air[spp][time_period]+=farquhar_results[0];
            tf_soil[spp][time_period]+=farquhar_results[5];
            An[spp][time_period][light_index]+=farquhar_results[1];
            E[spp][time_period][light_index]+=farquhar_results[2];
            Anb[spp][time_period][light_index]+=farquhar_results[3];
            Eb[spp][time_period][light_index]+=farquhar_results[4];

            }
            
            tf_air[spp][time_period]/=24.0;
            tf_soil[spp][time_period]/=24.0;
            An[spp][time_period][light_index]*=3600.0*360.0;
            E[spp][time_period][light_index]*=3600.0*540.0;
            Anb[spp][time_period][light_index]*=3600.0*360.0;
            Eb[spp][time_period][light_index]*=3600.0*540.0;
            //printf("Finished cal actual lat %d lon %d mechyear %d timestamp %d co2 %f\n",globY_,globX_,data->mechanism_year,data->MERRA2_timestamp-1,time_period,spp,CO2);

    
//#endif //COUPLE_MERRA2_LUT
            } //COUPLE_MERRA2_LUT
            //printf("lai %d lon %d mon %d lite %d An %f Anb %f E %f Eb %f\n",globY_,globX_,time_period,light_index,An[spp][time_period][light_index],Anb[spp][time_period][light_index],E[spp][time_period][light_index],Eb[spp][time_period][light_index]);
            //exit(0);
            
    return 1;
}

bool SiteData::farquhar (double Vmax,double CA, double ta, double ts,double ea, double q, double shade, int C4, double outputs[5]) {
    //printf("Cal farquhar V %f Ca %f ta %f ea %f q %f shade %f C4 %f\n",Vmax,CA,ta,ea,q,shade,C4);
    double shade_thresh, shade_thresh2;
    /* Scalers for turning boundary resitance */
    double gh_adj;
    double an;
    double g_adj;
    double capgam, vm, kc, ko, ds, v;
    double ci, tl, gg_adj, gg2, maxg, hta, hdriv, hl, lite;
    double f1, f2, f3, ginc, lasta, laste, lasttf, lastg, lastfunc, maxfunc;
    int flagv, flagg, gmax, j;
    double g;
    double M;
    double B=0.01;
    double rn;
    double tf, a, e, ab, eb;
    
    
    double GB=3.0,KAPPA=0.5,LAM=45000.0,CP=1280.1,GH=0.03,DO1=0.01,ALPHA3=0.08,ALPHA4=0.06;
    int PRECISION=10,BND=0;
    
    if (C4)
    {
        M = 4.0;
    }
    else
    {
        M = 8.0;
    }
    
    /* Meteorological conditions */
    rn = q * shade;
    
    
    /* Conversion of Solar influx from W/(m^2s) to mole_of_quanta/(m^2s) PAR, *
     * empirical relationship from McCree is lite = rn * 0.0000023            *
     * Prentice's relationship lite = rn / 540000.0 is ~25% - 30% lower       *
     * lite = rn * 0.0000023; converts watts short wave to Einsteins PAR      *
     * 1/540000 ~ 0.0000018                                                   */

    lite = rn * 0.0000023; /* assuming par 400-700nm, 0.5 sw is par */
    
    hta = 2541400.0 * exp(-5415.0 / (ta + 273.2)); //mol per mol humidity of air if it were saturated
    hdriv = 5415.0 * hta / ((ta + 273.2) * (ta + 273.2)); //Note in tropical this had an extra *hta in denomintaor - I belive this to be incorrect
    
    ginc = 1.0;
    gmax = B + 10.0 * ginc;
    g = B;
    
    shade_thresh = exp(-1.0 * KAPPA * 1.);
    shade_thresh2 = exp(-1.0 * KAPPA * 4.);

    v = Vmax;
    
    if(C4) {
        gg2 = 0.04;// gg_adj / v;
    } else {
        gg2 = 0.02;
    }
    
    flagv = 0;
    for (j=0; j<PRECISION; j++) {
        if (flagv == 0)
            flagg = 0;
        else
            flagg = 1;
        while ((flagg == 0) && (g < gmax)) {
            if (BND) {
                g_adj = GB / (g + GB);
                gh_adj = GB / 1.275 * 0.029 / GH;
            } else {
                g_adj = 1.;
                gh_adj = 1.;
            }
            tl = ta + (rn * KAPPA - LAM * g * g_adj * (hta - ea)) / (CP * GH * gh_adj + LAM * g * g_adj * hdriv);
            if (tl < ta)
                tl = ta;
            hl = 2541400.0 * exp(-5415.0 / (tl + 273.2));
            
            /* CO2 compenstaion point */
            capgam = 0.209 / (9000.0 * exp(-5000.0 * (1.0 / 288.2 - 1.0 / (tl + 273.2)))); //Equation 3 in Foley
            
            /* Temperature dependent Michaelis-Menten coefficients */
            ko = 0.25 * exp(1400.0 * (1.0 / 288.2 - 1.0 / (tl + 273.2)));
            kc = 0.00015 * exp(6000.0 * (1.0 / 288.2 - 1.0 / (tl + 273.2)));
            
            /* Saturation defficit */
            ds = hl - ea;
            
            /* Temperature dependence of Vmax */
            vm = v * exp(3000.0 * (1.0 / 288.2 - 1.0 / (tl + 273.2))); //A2 ED or 12 Foley
            
            /* Cold and warm shutdown- 092998 - GCH FROM SWP */

            if (C4)
            /* 4/25/00 coorected t func from foley to limit c4s */
                vm /= (1.0 + exp(0.4 * (10.0 - ta))) * (1.0 + exp(0.4 * (ta - 50.0)));
            else
                vm /= (1.0 + exp(0.4 * (5.0 - ta))) * (1.0 + exp(0.4 * (ta - 45.0)));
            
            /* Substomatal moler fraction of CO2 */
            if (BND) {
                ci = (g - B) * (1 + ds / DO1) * (capgam * (g + GB) - CA * GB) + M * CA * GB * g / 1.6;
                ci = ci / ((g - B) * (1 + ds / DO1) * g + M * GB * g / 1.6);
            } else {
                ci = (g - B) * (CA - capgam) * (1.0 + ds / DO1) / M;//13 Foley
                ci = CA - 1.6 * ci / g; //15 foley - note no boundary layer CO2 conductance
            }
            if (C4) {
                f1 = ALPHA4 * lite * KAPPA * 1.00;
                f2 = vm;
                f3 = 18000.0 * vm * ci;
                
                if (f1 < f2)
                    an = f1;
                else
                    an = f2;
                
                if (f3 < an)
                    an = f3;
                
            } else {
                /* Limiting rates of gross photosynthesis for C3 plants */
                f1 = ALPHA3 * lite * 1.0 * KAPPA * (ci - capgam) / (ci + 1.0 * capgam); //A1 in ED (main bit)
                f2 = vm * (ci - capgam) / (ci + kc * (1.0 + 0.209 / ko));
                
                if (f1 < f2)
                    an = f1;
                else
                    an = f2;
            }
            
            /* Dark "day" respiration */

            an -= gg2 * vm; /* this must equal gg * conts * tf */
            
            
            if (BND)
                f1 = an - g * g_adj * (CA - ci);
            else
                f1 = an - g * g_adj * (CA - ci) / 1.6;
            
            if (g == B) {
                maxfunc = f1;
                maxg = B;
                e = g * g_adj * ds;

                a = -gg2 * vm;
                tf = exp(3000.0 * (1.0 / 288.2 - 1.0 / (ta + 273.2)));
                
                if (C4)
                /* 4/25/00 coorected t func from foley to limit c4s */
                    tf /= (1.0 + exp(0.4 * (10.0 - ta))) * (1.0 + exp(0.4 * (ta - 50.0)));
                else
                    tf /= (1.0 + exp(0.4 * (5.0 - ta))) * (1.0 + exp(0.4 * (ta - 45.0)));
                
                /* tl changed to ta and denominator added  by gch 09/29/99 */
                eb = e;
                ab = a;
                flagv = 1;
                lastg = g;
                lastfunc = f1;
                laste = e;
                lasta = an;
                lasttf = tf;
                if (lite < 0.00000001)
                    flagg = 1;
                g += ginc / exp(j * log(10.0));
            } else {
                if (((f1 < 0.0) && (lastfunc >= 0.0)) || ((f1 > 0.0) && (lastfunc <= 0.0))) {
                    maxg = lastg;
                    maxfunc = lastfunc;
                    e = laste;
                    a = lasta;
                    tf = lasttf;
                    flagg = 1;
                    flagv = 0;
                } else {
                    lastg = g;
                    lastfunc = f1;
                    laste = g * g_adj * ds;
                    lasta = an;
                    lasttf = exp(3000.0 * (1.0 / 288.2 - 1.0 / (ta + 273.2)));
                    
                    if (C4)
                        lasttf /= (1.0 + exp(0.4 * (10.0 - ta))) * (1.0 + exp(0.4 * (ta - 50.0)));
                    else
                    /* tl changed to ta and denominator added, gch 09/29/99 */
                        lasttf /= (1.0 + exp(0.4 * (5.0 - ta))) * (1.0 + exp(0.4 * (ta - 45.0)));
                    
                    g += ginc / exp(j * log(10.0));
                }
            }
        }
        g = maxg + ginc / exp((j + 1.0) * log(10.0));
    }
    //Add separate rate for root respiration and saved in the 6th element in outputs array
    double tf_soil=0;
    tf_soil = exp(3000.0 * (1.0 / 288.2 - 1.0 / (ts + 273.2)));
    if (C4)
        tf_soil /= (1.0 + exp(0.4 * (10.0 - ts))) * (1.0 + exp(0.4 * (ts - 50.0)));
    else
    /* tl changed to ta and denominator added, gch 09/29/99 */
        tf_soil /= (1.0 + exp(0.4 * (5.0 - ts))) * (1.0 + exp(0.4 * (ts - 45.0)));
    
    outputs[0] = tf; outputs[1] = a; outputs[2] = e; outputs[3] = ab; outputs[4] = eb;outputs[5] = tf_soil;
    //printf("capgam %f kc %f ko %f vm %f \n",capgam,kc,ko,vm);
    //printf("Tl %f Ta %f Ag %f An %f light %f Ci %f Ds %f\n",tl,ta,a,an,lite,ci,ds);
    return 1;
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
   for (size_t i=0;i<N_LAI;i++)
   {
       cs->site_lai_profile[i]=0;
   }
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
   for (size_t i=0;i<N_LAI;i++)
   {
       cs->site_aa_lai_profile[i]=0;
   }
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
    
    cs->site_avg_fopen=0.0;
    cs->site_avg_leafAn_pot=0.0;
    cs->site_avg_leafE_pot=0.0;
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
      for (size_t i=0;i<N_LAI;i++)
      {
          cs->site_lai_profile[i]+=cs->lai_profile[lu][i]*cs->area_fraction[lu];
      }
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
      for (size_t i=0;i<N_LAI;i++)
      {
         cs->site_aa_lai_profile[i]+=cs->aa_lai_profile[i]*cs->area_fraction[lu];
      }
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
       
       cs->site_avg_leafAn_pot+=cs->lu_avg_leafAn_pot[lu]*cs->area_fraction[lu];
       cs->site_avg_leafE_pot+=cs->lu_avg_leafE_pot[lu]*cs->area_fraction[lu];
       cs->site_avg_fopen+=cs->lu_avg_fopen[lu]*cs->area_fraction[lu];
       
      cs->site_npp2               += cs->npp2[lu] * cs->area_fraction[lu];
      cs->site_aa_gpp             += cs->aa_gpp[lu] * cs->area_fraction[lu];
      cs->site_aa_npp2            += cs->aa_npp2[lu] * cs->area_fraction[lu];
      if (data->time_period == 0) {
         cs->nep2[lu] = cs->total_c[lu] - cs->last_total_c[lu];
         cs->last_total_c[lu] = cs->total_c[lu];
      }
       cs->nep3[lu]=cs->total_c[lu] - cs->last_total_c_last_month[lu];
       cs->last_total_c_last_month[lu]=cs->total_c[lu];
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
    cs->site_nep3 = cs->site_total_c- cs->last_site_total_c_last_month;
    cs->last_site_total_c_last_month=cs->site_total_c;
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
    for (size_t i=0;i<N_LAI;i++)
    {
        cs->aa_lai_profile[i]=0;
    }
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
   cs->lu_avg_fopen[lu]          = 0.0;
   cs->lu_avg_leafAn_pot[lu]     = 0.0;
   cs->lu_avg_leafE_pot[lu]       = 0.0;
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
    for (size_t i=0;i<N_LAI;i++)
    {
        cs->lai_profile[lu][i]=0;
    }
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
          for (size_t i=0;i<N_LAI;i++)
          {
              cs->aa_lai_profile[i]+=cp->aa_lai_profile[i]*frac;
          }
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
          
          cs->lu_avg_leafAn_pot[lu]+=cp->avg_leafAn_pot*frac;
          cs->lu_avg_leafE_pot[lu]+=cp->avg_leafE_pot*frac;
          cs->lu_avg_fopen[lu]+=cp->avg_fopen*frac;
          
         cs->npp2[lu]               += cp->npp2 * frac;
         cs->aa_gpp[lu]             += cp->aa_gpp * frac;
         cs->aa_npp2[lu]            += cp->aa_npp2 * frac;
         /* water */
         cs->water[lu]              += cp->water * frac;
         cs->theta[lu]              += cp->theta * frac;
         /* TODO: why are these different? -justin */
         cs->total_water_uptake[lu] += cp->total_water_uptake / (cs->area_fraction[lu]*data->area);
         cs->perc[lu]               += cp->perc * frac;
          //Modified add check
          if (cp->perc>1e18)
          {
              printf("Wrong in perc 1: lat-%f lon-%f perc-%f lu-%d frac-%f\n",cs->sdata->lat_,cs->sdata->lon_,cp->perc,lu,frac);
              //exit(1);
          }
         cs->soil_evap[lu]          += cp->soil_evap * frac;
         /* height */
         /* insulate against empty patches */
         if (cp->tallest != NULL)
            cs->lu_avg_height[lu]   += cp->tallest->hite * frac;

         cs->basal_area[lu]         += cp->basal_area * frac;
         cs->lai[lu]                += cp->lai * frac;
          for (size_t i=0;i<N_LAI;i++)
          {
              cs->lai_profile[lu][i]+=cp->lai_profile[i]*frac;
          }
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

