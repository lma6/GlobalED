#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "edmodels.h"
#include "site.h"
#include "patch.h"
#include "read_site_data.h"
#include "cohort.h"


#define R_gas 8.314  //!< \b R idealgasconstant
#define maxiter 200 //!< \b maxiter maximum number of iterations
#define epsilon 0.97   //!< epsilon emissivity See Campbell and Norman, 1998, page 163 (CHECK)
#define sbc 5.6697e-8  //!< stefan-Boltzmann constant Wm-2 k-4. Actually varies somewhat with temperature
#define scatt 0.15     //!< leaf reflectance + transmittance
#define f_abs 0.15      //!< spectral correction
#define O 205.0     //!<  Oxygen partial pressure gas units are mbar
#define Q10 2.0        //!< Q10 factor

//test_larch
#define parameterizationCase 4 //1-Bonan 2011; 2-Kattge 2007 without thermal accmlimation; 3-Kattge 2007 with acclimation;
                               //4-Kattge 2007 with acclimation but Tg is 40
                               //5-Kattge 2007 with PFT dependent

//Reference of this rountine is from Bonan et al 2011 including temperature dependency of Vcmax, Jcmax and Rd.
//Parameterization from Kattge 2007
void SiteData::Initilize(int pt,int spp,double Tg,UserData* data)
{
    C4=pt;
    errTolerance = 0.001;
    eqlTolerance = 1.0e-6;
    
    PhotoFluxDensity=0,  //!< Photosynthetic Flux Density (umol photons m-2 s-1
    R_abs=0, //!< Absorbed incident radiation (watts m-2)
    Tair=0,  //!< Air temperature at 2m, (C)
    CO2=0,   //!< CO2 concentration (umol mol-1 air)
    RH=0,   //!<  Relative Humidity (%, i.e., 80)
    wind=0, //!<  Windspeed at 2 meters (km s-1)
    width=0, //!< Leaf width (m)
    Press=0;  //!<  Air pressure (kPa)
    
    AssimilationNet=0,    //!< Net photosynthesis (umol CO2 m-2 s-1)
    AssimilationGross=0, //!< Gross photosynthesis (umol CO2 m-2 s-1) (Adjusted for respiration)
    Transpiration=0,     //!< Transpiration mol H2O m-2 s-1
    Tleaf=0,  //!< Leaf temperature C
    Ci=0,     //!< Internal CO2 concentration umol mol-1
    StomatalConductance=0,     //!< Stomatal conductance umol H2O m-2 s-1
    BoundaryLayerConductance=0,    //!< Boundary layer conductance umol H2O m-2 s-1
    DarkRespiration=0,    //!< Plant respiration    umol m-2 s-1
    VPD=0,    //!< Vapor Pressure Density, kPa */
    Ci_Ca=0;  //!< Ratio of internal to external CO2, unitless
    
    Theta=0.7;
    
    //Kattge et al 2007 used here for parameterization with temperature acclimation
    if (C4==0)
    {
#if parameterizationCase==1
        EaVc=65330;  //Need to find value, now use same with C3
        Eaj=43540;   //Delta_Ja
        Ear=46390;   //66400
        Eap=53100;
        
        Sv=485;
        Sj=495;
        Sr=490;
        Sp=490;
        
        Hv=149250;
        Hj=152040;
        Hr=150650;
        Hp=150650;
        
        rJ2V=1.97;
        Q10R=1; //from Atkin et al 2008
#endif
        
#if parameterizationCase==2
        EaVc=71513;  //Kattge 2007 Table.3
        Eaj=49884;   //Kattge 2007 Table.3
        Ear=66400;   //Need to check. from photo...master paramater.xls
        Eap=53100;   //As no acclimation in Kattge 2007, same as Bonan 2011
        
        Sv=649.12;
        Sj=646.22;
        Sr=490;     //As no acclimation in Kattge 2007, same as Bonan 2011, acclimation is achieved using Q10R factor
        Sp=490;     //As no acclimation in Kattge 2007, same as Bonan 2011
        
        Hv=200000;  //Kattge 2007 Table.3
        Hj=200000;  //Kattge 2007 Table.3
        Hr=150650;  //As no acclimation in Kattge 2007, same as Bonan 2011
        Hp=150650;  //As no acclimation in Kattge 2007, same as Bonan 2011
        
        rJ2V=1.97; //Kattge 2007 Table.3
        
        Q10R=1; //from Atkin et al 2008
#endif
        
#if parameterizationCase==3
        //printf("Using thermal acclimation based parameterization\n");
        //As Kattge only considered temperature range from 11 to 35 degrees Celsius
        if (Tg>35)  Tg=35;
        if (Tg<11)  Tg=11;
        
        
        EaVc=71513;  //Kattge 2007 Table.3
        Eaj=49884;   //Kattge 2007 Table.3
        Ear=66400;   //Need to check. from photo...master paramater.xls
        Eap=53100;   //As no acclimation in Kattge 2007, same as Bonan 2011
        
        Sv=668.39-1.07*Tg;
        Sj=659.70-0.75*Tg;
        Sr=490;     //As no acclimation in Kattge 2007, same as Bonan 2011, acclimation is achieved using Q10R factor
        Sp=490;     //As no acclimation in Kattge 2007, same as Bonan 2011
        
        Hv=200000;  //Kattge 2007 Table.3
        Hj=200000;  //Kattge 2007 Table.3
        Hr=150650;  //As no acclimation in Kattge 2007, same as Bonan 2011
        Hp=150650;  //As no acclimation in Kattge 2007, same as Bonan 2011
        
        rJ2V=2.59-0.035*Tg; //Kattge 2007 Table.3
        
        Q10R=pow(10,-0.00794*(Tg-25.0)); //from Atkin et al 2008, Eq 8, C value from log-log Rmt-N plots
        
//        if (!strcmp(data->title[spp],"evergreen"))
//            rJ2V=1.07;
#endif
        
#if parameterizationCase==4
        //test_larch
        //Tg = 40.0;
        Tg = 30.0;
        
        EaVc=71513;  //Kattge 2007 Table.3
        Eaj=49884;   //Kattge 2007 Table.3
        Ear=66400;   //Need to check. from photo...master paramater.xls
        Eap=53100;   //As no acclimation in Kattge 2007, same as Bonan 2011
        
        Sv=668.39-1.07*Tg;
        Sj=659.70-0.75*Tg;
        Sr=490;     //As no acclimation in Kattge 2007, same as Bonan 2011, acclimation is achieved using Q10R factor
        Sp=490;     //As no acclimation in Kattge 2007, same as Bonan 2011
        
        Hv=200000;  //Kattge 2007 Table.3
        Hj=200000;  //Kattge 2007 Table.3
        Hr=150650;  //As no acclimation in Kattge 2007, same as Bonan 2011
        Hp=150650;  //As no acclimation in Kattge 2007, same as Bonan 2011
        
        rJ2V=2.59-0.035*Tg; //Kattge 2007 Table.3
        
        Q10R=pow(10,-0.00794*(Tg-25.0)); //from Atkin et al 2008
        
        //test_larch
        //if (!strcmp(data->title[spp],"evergreen"))
//         if (!strcmp(data->title[spp],"evergreen") || !strcmp(data->title[spp],"larch"))
//            rJ2V=1.07;
#endif
        
//test_larch
#if parameterizationCase==5
        if(!strcmp(data->title[spp],"pine"))  //test_mor, should be evergreen_short
           Tg = 30.0;
        if(!strcmp(data->title[spp],"late_succ_conifer"))  //test_mor, should be evergreen_long
            Tg = 10.0;
        if(!strcmp(data->title[spp],"evergreen"))
            Tg = 10.0;
        if(!strcmp(data->title[spp],"cold_decid"))
            Tg = 30.0;
        if(!strcmp(data->title[spp],"c3_grass"))
            Tg = 30.0;
        if(!strcmp(data->title[spp],"c4_grass"))
            Tg = 30.0;
        if(!strcmp(data->title[spp],"early_succ"))
            Tg = 30.0;
        if(!strcmp(data->title[spp],"mid_succ"))
            Tg = 30.0;
        if(!strcmp(data->title[spp],"late_succ"))
            Tg = 30.0;
           
        EaVc=71513;  //Kattge 2007 Table.3
        Eaj=49884;   //Kattge 2007 Table.3
        Ear=66400;   //Need to check. from photo...master paramater.xls
        Eap=53100;   //As no acclimation in Kattge 2007, same as Bonan 2011
           
        Sv=668.39-1.07*Tg;
        Sj=659.70-0.75*Tg;
        Sr=490;     //As no acclimation in Kattge 2007, same as Bonan 2011, acclimation is achieved using Q10R factor
        Sp=490;     //As no acclimation in Kattge 2007, same as Bonan 2011
           
        Hv=200000;  //Kattge 2007 Table.3
        Hj=200000;  //Kattge 2007 Table.3
        Hr=150650;  //As no acclimation in Kattge 2007, same as Bonan 2011
        Hp=150650;  //As no acclimation in Kattge 2007, same as Bonan 2011
           
        rJ2V=2.59-0.035*Tg; //Kattge 2007 Table.3
           
        Q10R=pow(10,-0.00794*(Tg-25.0)); //from Atkin et al 2008
#endif
        
    }
    else
    {
        EaVc=55900;
        EaVp=75100;
        Ear=39800;
        Eaj=32800;
        
        Sj=702.6;
        Hj=220000;
    }

    iter_total=0,
    Jm25=0,
    Vpm25=0,
    TPU25=0,
    Rd25=0,
    g0=0.01, // test_larch 0.02
    g1=10;
    stomaRatio=0.5;
    LfWidth=0.05;
    LfAngFact=1;
}

//ea- Specfific hummidity mol mol-1
//shortwaveRad: shortwave radiation w/m2
//Ca: ambient CO2 centration
//Vm25 Maxmimum Rubisco campacity at 25 ceilus
//wind: Windspeed at 2.5 m, m s-1
//Press Atmospheric pressure (kpa m-2)
void SiteData::Farquhar_couple(int pt, int spp,UserData* data,double Ta, double Ts,double ea, double swd, double Tg,double Ca, double speed, double Pa, double shade, double Vm25, double outputs[6])
{
    Initilize(pt,spp,Tg,data);
    
    if(C4==0)
    {
        Vcm25=Vm25;
        Jm25=Vcm25*rJ2V;      //Medlyn 2002 Fig.3
        TPU25=Jm25*0.06;
        Rd25=Vcm25*0.015;
        
        //test_larch
        //Atkin et al 2005 suggests PFTs in cold site has high dark respiration than those in warm sites -- Lei
//        if(spp==5)
//            Rd25=Vcm25*0.045;
        
    }else
    {
//        Vcm25=Vm25;
//        Jm25=Vcm25*6.0;
//        Vpm25=Vcm25*1.4;
//        Rd25=Vcm25*0.015;
        
        Vcm25=Vm25;
        Jm25=Vcm25*5;
        Vpm25=Vcm25*1.4;
        Rd25=Vcm25*0.01;
    }
    
    
    PhotoFluxDensity = swd*shade*0.5*4.55;
    Tair=Ta;
    RH=ea/(0.611*exp(17.502*Tair/(240.97+Tair))/1e2);  //Convert specfific humidity (mol/mol) to relative humidify (0-100%)
    RH=fmin(1.0, RH);
    wind=speed;
    Press=Pa;
    CO2=Ca;
    
    double PAR = (PhotoFluxDensity/4.55); //PAR is watts m-2
    double NIR = PAR; // If total solar radiation unavailable, assume NIR the same energy as PAR waveband
    R_abs = (1-scatt)*PAR + 0.15*NIR + 2*(epsilon*sbc*pow(Tair+273,4)); // times 2 for projected area basis
    // shortwave radiation (PAR (=0.85) + NIR (=0.15) solar radiation absorptivity of leaves: =~ 0.5
    //transfer variables to local scope

    GasEx();   // Gas exchange calculations here
    
    double tf=exp(3000.0 * (1.0 / 288.2 - 1.0 / (Tair + 273.2)));
    if (C4==1)
        tf /= (1.0 + exp(0.4 * (10.0 - Tair))) * (1.0 + exp(0.4 * (Tair - 50.0)));
    else
        tf /= (1.0 + exp(0.4 * (5.0 - Tair))) * (1.0 + exp(0.4 * (Tair - 45.0)));
    

    outputs[5]=exp(3000.0 * (1.0 / 288.2 - 1.0 / (Ts + 273.2)));
    if (C4==1)
        outputs[5] /= (1.0 + exp(0.4 * (10.0 - Ts))) * (1.0 + exp(0.4 * (Ts - 50.0)));
    else
        outputs[5] /= (1.0 + exp(0.4 * (5.0 - Ts))) * (1.0 + exp(0.4 * (Ts - 45.0)));
    
    outputs[0]=tf;
    outputs[1]=AssimilationNet;
    outputs[2]=Transpiration;
    outputs[3]=DarkRespiration;
    outputs[4]=g0*BoundaryLayerConductance/(g0+BoundaryLayerConductance)*(Es(Tleaf)-Es(Tair)*RH)/Press;//Transpiration when stomatal closed, may need to be tested.
    
    //test_mor
//    for (int i=0;i<5;i++)
//    {
//        if (outputs[i]<0 or outputs[i]>5e3)  outputs[i]=0;
//    }
    outputs[1]/=1e6;
    outputs[3]/=-1e6;
    
    
//    double Ds=(Es(Tleaf)-RH*Es(Tair))/Press;
//    double Jmax = Jm25*exp(((Tleaf-25)*Eaj)/(R_gas*(Tleaf+273)*298))*
//    (1+exp((Sj*298-Hj)/(R_gas*298)))/
//    (1+exp((Sj*(Tleaf+273)-Hj)/(R_gas*(Tleaf+273)))); // de Pury 1997
//    double Vcmax = Vcm25*exp(((Tleaf-25)*EaVc)/(R_gas*(Tleaf+273)*298))*
//    (1+exp((Sv*298-Hv)/(R_gas*298)))/
//    (1+exp((Sv*(Tleaf+273)-Hv)/(R_gas*(Tleaf+273)))); // Used peaked response, DHF
//    //TPU = TPU25*exp(Eap*(Tleaf-25)/(298*R_gas*(Tleaf+273)));  //orginal one, does account thermal breakdown
//    double TPU = TPU25*exp(((Tleaf-25)*Eap)/(R_gas*(Tleaf+273)*298))*
//    (1+exp((Sp*298-Hp)/(R_gas*298)))/
//    (1+exp((Sp*(Tleaf+273)-Hp)/(R_gas*(Tleaf+273))));
//
//
//    const long Lambda = 44000; //latent heat of vaporization of water J mol-1 - not used in this implementation
//    const double Cp = 29.3; // thermodynamic psychrometer constant and specific hear of air, J mol-1 C-1
//    const double psc = 6.66e-4; //psycrometric constant units are C-1
//
//    double HeatConductance,  //heat conductance J m-2 s-1
//    VaporConductance, //vapor conductance ratio of stomatal and heat conductance mol m-2 s-1
//    RadiativeConductance, //radiative conductance J m-2 s-1
//    RadiativeAndHeatConductance, //radiative+heat conductance
//    psc1,  // apparent psychrometer constant Campbell and Norman, page 232 after eq 14.11
//    Ea,   //ambient vapor pressure kPa
//    thermal_air; // emitted thermal radiation Watts  m-2
//
//
//
//    HeatConductance = BoundaryLayerConductance*(0.135/0.147);  // heat conductance, HeatConductance = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Boundary Layer Conductance to Heat
//    // Since BoundaryLayerConductance is .147*sqrt(u/d) this scales to 0.135*sqrt(u/d) - HeatConductance on page 109 of Campbell and Norman, 1998
//    // Wind was accounted for in BoundaryLayerConductance already  as BoundaryLayerConductance (turbulent vapor transfer) was calculated from CalcTurbulentVaporConductance() in GasEx.
//    // units are J m-2 s-1
//    VaporConductance = StomatalConductance*BoundaryLayerConductance/(StomatalConductance+BoundaryLayerConductance);      //vapor conductance, StomatalConductance is stomatal
//    double ds=(Es(Tleaf)-Ea)/Press;
//
//    double thermal_leaf=epsilon*sbc*pow(Tleaf+273,4)*2;
//
//    double Res = R_abs - thermal_leaf - Cp*HeatConductance*(Tleaf - Tair) - Lambda*VaporConductance*1.0*(Es(Tleaf)-Ea)/Press; // Residual function: f(Ti), KT Paw (1987)
//
//    printf("swd %f R_abs %f thermLeaf %f Rn %f H %f lamdaE %f Res %f Tl %f Ta %f gsw %f gbw %f gH %f ds %f\n",swd,R_abs,thermal_leaf,R_abs - thermal_leaf,Cp*HeatConductance*(Tleaf - Tair),Lambda*VaporConductance*1.0*(Es(Tleaf)-Ea)/Press,Res,Tleaf,Tair,StomatalConductance,BoundaryLayerConductance,HeatConductance,ds);
    
    
    
//    printf("Tl %f Ta %f Vcmax %f Jmax %f TPU %f Ag %f An %f light %f Ci %f Ds %f gb %f VPD %f\n",Tleaf,Tair,Vcmax,Jmax,TPU,AssimilationGross,AssimilationNet,PhotoFluxDensity,Ci,Ds*Press,BoundaryLayerConductance,this->VPD);
   
    Initilize(pt,spp,Tg,data);
}

void SiteData::GasEx(void)
{
    double Tleaf_old;  //previous leaf temperture (for iteration)
    int iter=1;
    iter_total=0;
    Tleaf = Tair; Tleaf_old = -1000;
    Ci = 0.7*CO2;
    BoundaryLayerConductance = CalcTurbulentVaporConductance();
    StomatalConductance = CalcStomatalConductance();
    while ((fabs(Tleaf_old -Tleaf)>0.01) && (iter < maxiter))
    {
        Tleaf_old=Tleaf;
        Ci=SearchCi(Ci);
        StomatalConductance=CalcStomatalConductance();
        EnergyBalance();
        iter2 =++iter; //iter=iter+1, iter2=iter;
    }
    
}


void SiteData::PhotosynthesisC3(double Ci)
{
//    for (Tleaf=0;Tleaf<50;Tleaf++)
//    {
//        PhotoFluxDensity=1350;
//        Ci=350;
    //parameters for C3 Photosythesis;
    const double curvature=0.999 ; //!< \b Curvature -factor of Av and Aj colimitation
    
    const int        Kc25 = 404.4;//!< \b Kc25, MM Constant of rubisco for CO2 of C3 plants (de Pury and Farquar, 1997) (umol m-2 s-1)
    const int        Ko25 = 278.4;//!< \b Ko25, MM Constant of rubiscuo for O2 from above reference (umol m-2 s-1)
    const long      Eac = 79430;//!< \b Eac, Energy Activation kJ mol-1
    
    const long       Eao = 36380;//!< \b Eao, activation energy values
    //** \endcode
    // These variables hold temporary calculations
    double alpha, Kc, Ko, gamma, Ia,Jmax, Vcmax, TPU, J, Av, Aj, Ap, Ac, Km, Ca, Cc, P;
    //gamma = 36.9 + 1.88*(Tleaf-25)+0.036*Square(Tleaf-25);  // CO2 compensation point in the absence of mitochondirial respiration, in ubar}
    gamma=42.75*exp(37830*(Tleaf-25)/298/R_gas/(Tleaf+273));
    
    //* Light response function parameters */
    Ia = PhotoFluxDensity*(1-scatt);    //* absorbed irradiance */
    alpha = (1-f_abs)/2; // *!apparent quantum efficiency, params adjusted to get value 0.3 for average C3 leaf
    
    AssimilationNet = 0;
    
    //* other input parameters and constants */
    P  = Press/100; //Press is kPa. Used to convert mole fraction to partial pressure
    Ca = CO2*P; //* conversion to partial pressure */
    Kc = Kc25*exp(Eac*(Tleaf-25)/(298*R_gas*(Tleaf+273)));
    Ko = Ko25*exp(Eao*(Tleaf-25)/(298*R_gas*(Tleaf+273)));
    Km = Kc*(1+O/Ko); //* effective M-M constant for Kc in the presence of O2 */

    DarkRespiration=Rd25*exp(((Tleaf-25)*Ear)/(R_gas*(Tleaf+273)*298))*
    (1+exp((Sr*298-Hr)/(R_gas*298)))/
    (1+exp((Sr*(Tleaf+273)-Hr)/(R_gas*(Tleaf+273))));
    DarkRespiration*=Q10R;
    
    Jmax = Jm25*exp(((Tleaf-25)*Eaj)/(R_gas*(Tleaf+273)*298))*
    (1+exp((Sj*298-Hj)/(R_gas*298)))/
    (1+exp((Sj*(Tleaf+273)-Hj)/(R_gas*(Tleaf+273)))); // de Pury 1997
    Vcmax = Vcm25*exp(((Tleaf-25)*EaVc)/(R_gas*(Tleaf+273)*298))*
    (1+exp((Sv*298-Hv)/(R_gas*298)))/
    (1+exp((Sv*(Tleaf+273)-Hv)/(R_gas*(Tleaf+273)))); // Used peaked response, DHF
    //TPU = TPU25*exp(Eap*(Tleaf-25)/(298*R_gas*(Tleaf+273)));  //orginal one, does account thermal breakdown
    TPU = TPU25*exp(((Tleaf-25)*Eap)/(R_gas*(Tleaf+273)*298))*
    (1+exp((Sp*298-Hp)/(R_gas*298)))/
    (1+exp((Sp*(Tleaf+273)-Hp)/(R_gas*(Tleaf+273))));
    Cc = Ci; // assume infinite gi
    
    StomatalConductance = CalcStomatalConductance(); // Initial value
    BoundaryLayerConductance=  CalcTurbulentVaporConductance();
    Av = (Vcmax*(Cc-gamma))/(Cc+Km);
    J =  (((alpha*Ia + Jmax) - sqrt(Square(alpha*Ia+Jmax) - 4*alpha*Ia*(Jmax)*Theta)) / (2*Theta)) ;
    Aj = J*(Cc-gamma)/(4*(Cc+2*gamma));
    Ap = 3*TPU;
    Ac = ((Av+Aj) - sqrt(Square(Av+Aj)-4*curvature*Av*Aj))/(2*curvature); // curvatureaccount for colimitation between Av and Aj */
    
    if (Cc > gamma)
        AssimilationNet = fmin(Ac, Ap) -DarkRespiration;
        //AssimilationNet = Ac -DarkRespiration;  //Turn off TPU limitation
    else
    {
        AssimilationNet = Av-DarkRespiration;
    }
    
    
    AssimilationGross = fmax(AssimilationNet+DarkRespiration,0.0);
    StomatalConductance = CalcStomatalConductance(); // Update StomatalConductance using new value of AssimilationNet
    
//    printf("gb %f\n",BoundaryLayerConductance);
//    printf("Tl %f Av %f Aj %f Ap %f Ac %f Ag %f An %f Vcmax %f Jmax %f light %f Ci %f Cc %f gamma %f Km %f Kc %f Ko %f\n",Tleaf,Av,Aj,Ap,Ac,AssimilationGross,AssimilationNet,Vcmax,Jmax,PhotoFluxDensity,Ci,Cc,gamma,Km,Kc,Ko);
//    exit(0);
}

void SiteData::PhotosynthesisC4(double Ci)
{
    const double    curvature=0.995; //!< \b curvature factor of Av and Aj colimitation
    
    
    const int       Kc25 = 650,    //!< \b Kc25, Michaelis constant of rubisco for CO2 of C4 plants (2.5 times that of tobacco), ubar, Von Caemmerer 2000
    Ko25 = 450,                //!< \b Ko25, Michaelis constant of rubisco for O2 (2.5 times C3), mbar
    Kp25 = 57;                 //*!< \b Kp25, Michaelis constant for PEP caboxylase for CO2 - was 60 in Kim's paper */
    const long       Eao = 36000;  //*!< \b EAO, activation energy for Ko */
    const int        Vpr25 = 80;    //*!<   \b Vpr25, PEP regeneration limited Vp at 25C, value adopted from vC book */
    const double    gbs = 0.003;   //*!< \b gbs, bundle sheath conductance to CO2, umol m-2 s-1 gbs x Cm is the inward diffusion of CO2 into the bundle sheath  */
    const double    x = 0.4;       //*!< \b x Partitioning factor of J, yield maximal J at this value */
    const double    alpha = 0.001; //*!< \b alpha, fraction of PSII activity in the bundle sheath cell, very low for NADP-ME types  */ 0.001
    const double    gi = 5.0;      //*!< \b gi, conductance to CO2 from intercelluar to mesophyle, mol m-2 s-1, assumed  was 1, changed to 5 as per Soo 6/2012*/
    const double    beta = 0.99;   //*!< \b beta, smoothing factor */
    const double    gamma1 = 0.193; //*!< \b gamma1, half the reciprocal of rubisco specificity, to account for O2 dependence of CO2 comp point, note that this become the same as that in C3 model when multiplied by [O2] */
    
    double Kp, Kc, Ko, Km;         //!<\b Kp, \b Kc, \b Ko, \b Km, Calculated Michaelis params as a function of temperature
    double Ia, I2;                 // secondary calculated light variables
    double Vpmax, Jmax, Vcmax, Eac, Om, Rm, J, Ac1, Ac2, Ac, Aj1,
    Aj2, Aj, Vp1, Vp2, Vp, P,  Ca, Cm, Vpr,
    Os, GammaStar, Gamma, a1, b1, c1; //secondary calculated variables
    
    //* Light response function parameters */
    Ia = PhotoFluxDensity*(1-scatt);    //* absorbed irradiance */
    I2 = Ia*(1-f_abs)/2;    //* useful light absorbed by PSII */
    //* other input parameters and constants */
    P  = Press/100;
    Ca = CO2*P; //* conversion to partial pressure Atmospheric partial pressure of CO2, kPa*/
    Om = O;   //* mesophyle O2 partial pressure */
    Eac=EaVc;
    
    Kp = Kp25*pow(Q10,(Tleaf-25.0)/10.0);
    Vpr = Vpr25*pow(Q10,(Tleaf-25.0)/10.0);
    Kc = Kc25*exp(Eac*(Tleaf-25)/(298*R_gas*(Tleaf+273))); //Kc adjusted for temperature
    Ko = Ko25*exp(Eao*(Tleaf-25)/(298*R_gas*(Tleaf+273)));
    Km = Kc*(1+Om/Ko); //* effective M-M constant for Kc in the presence of O2 */
    DarkRespiration = Rd25*exp(Ear*(Tleaf-25)/(298*R_gas*(Tleaf+273)));
    // The following are Arrhenius Equations for parameter temperature dependencies
    // Vpm25 (PEPC activity rate) , Vcm25  (Rubisco Capacity rate) and Jm25 (Whole chain electron transport rate) are the rates at 25C for Vp, Vc and Jm
    Vpmax = Vpm25*exp(EaVp*(Tleaf-25)/(298*R_gas*(Tleaf+273)));
    Vcmax = Vcm25*exp(EaVc*(Tleaf-25)/(298*R_gas*(Tleaf+273)));
    Jmax = Jm25*exp((((Tleaf+273)-298)*Eaj)/(R_gas*(Tleaf+273)*298))*(1+exp((Sj*298-Hj)/(R_gas*298)))
    /(1+exp((Sj*(Tleaf+273)-Hj)/(R_gas*(Tleaf+273.0))));
    Rm = 0.5*DarkRespiration;
    
    Cm=Ci; //* mesophyle CO2 partial pressure, ubar, one may use the same value as Ci assuming infinite mesohpyle conductance */
    double gs_last=0;
    
    StomatalConductance = CalcStomatalConductance();
    Vp1 = (Cm*Vpmax)/(Cm+Kp); //* PEP carboxylation rate, that is the rate of C4 acid generation  Eq 1 in Kim 2007*/
    Vp2 = Vpr;
    Vp = fmax(fmin(Vp1, Vp2),0);
    //* Enzyme limited A (Rubisco or PEP carboxylation */
    Ac1 = (Vp+gbs*Cm-Rm);
    Ac2 = (Vcmax-DarkRespiration);
    //* Quadratic expression to solve for Ac */
    a1 = 1-(alpha/0.047)*(Kc/Ko);
    b1 = -(Ac1 + Ac2 + gbs*Km + (alpha/0.047)*(gamma1*Vcmax + DarkRespiration*Kc/Ko));
    c1 = Ac1*Ac2-(Vcmax*gbs*gamma1*Om+DarkRespiration*gbs*Km);
    Ac = QuadSolnLower(a1,b1,c1);
    Ac = fmin(Ac1,Ac2);
    //* Light and electron transport limited  A mediated by J */
    J=minh(I2,Jmax,Theta);  //* rate of electron transport */
    Aj1 = (x*J/2-Rm+gbs*Cm);  // Eq 4 in Kim, 2007
    Aj2 = (1-x)*J/3-DarkRespiration;       //Eq 4 in Kim, 2007
    Aj = fmin(Aj1,Aj2);      //Eq 4 in Kim, 2007
    AssimilationNet = ((Ac+Aj) - sqrt(Square(Ac+Aj)-4*beta*Ac*Aj))/(2*beta); //* smooting the transition between Ac and Aj */
    AssimilationNet=minh(Ac,Aj, curvature);
    gs_last=StomatalConductance;
    Os = alpha*AssimilationNet/(0.047*gbs)+Om; //* Bundle sheath O2 partial pressure, mbar */
    GammaStar = gamma1*Os;
    Gamma = (DarkRespiration*Km + Vcmax*GammaStar)/(Vcmax-DarkRespiration);
    AssimilationGross = fmax(0, AssimilationNet + DarkRespiration);
    
//    printf("Ci-%f Tleaf-%f Vcm25-%f Vpm25-%f Jm25-%f Vcmax-%f Vpmax-%f Jmax-%f Rm-%f Vp1-%f Vp2-%f Ac1-%f Ac2-%f Aj1-%f Aj2-%f An-%f\n",Ci,Tleaf,Vcm25,Vpm25,Jm25,Vcmax,Vpmax,Jmax,Rm,Vp1,Vp2,Ac1,Ac2,Aj1,Aj2,AssimilationNet);
    
}

//Change Newton iteration solving Tleaf from orignal version saved in EnergyBalance2, which does not
//change newTi during iteration.
void SiteData::EnergyBalance()
{
    const long Lambda = 44000; //latent heat of vaporization of water J mol-1 - not used in this implementation
    const double Cp = 29.3; // thermodynamic psychrometer constant and specific hear of air, J mol-1 C-1
    const double psc = 6.66e-4; //psycrometric constant units are C-1
    //psc=Cp/Lambda = 29.3/44000 See Campbell and Norman, pg 232, after eq 14.11
    
    //The following are secondary variables used in the energy balance
    double HeatConductance,  //heat conductance J m-2 s-1
    VaporConductance, //vapor conductance ratio of stomatal and heat conductance mol m-2 s-1
    RadiativeConductance, //radiative conductance J m-2 s-1
    RadiativeAndHeatConductance, //radiative+heat conductance
    psc1,  // apparent psychrometer constant Campbell and Norman, page 232 after eq 14.11
    Ea,   //ambient vapor pressure kPa
    thermal_air; // emitted thermal radiation Watts  m-2
    double lastTi, newTi;
    int    iter;
    
    HeatConductance = BoundaryLayerConductance*(0.135/0.147);  // heat conductance, HeatConductance = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Boundary Layer Conductance to Heat
    // Since BoundaryLayerConductance is .147*sqrt(u/d) this scales to 0.135*sqrt(u/d) - HeatConductance on page 109 of Campbell and Norman, 1998
    // Wind was accounted for in BoundaryLayerConductance already  as BoundaryLayerConductance (turbulent vapor transfer) was calculated from CalcTurbulentVaporConductance() in GasEx.
    // units are J m-2 s-1
    VaporConductance = StomatalConductance*BoundaryLayerConductance/(StomatalConductance+BoundaryLayerConductance);      //vapor conductance, StomatalConductance is stomatal conductance and is given as gvs in Campbell and Norman.
    // note units are moles m-2 s-1.
    RadiativeConductance = (4*epsilon*sbc*pow(273+Tair,3)/Cp)*2; // radiative conductance, *2 account for both sides
    RadiativeAndHeatConductance = HeatConductance + RadiativeConductance;
    thermal_air = epsilon*sbc*pow(Tair+273,4)*2; //Multiply by 2 for both surfaces
    psc1 = psc*RadiativeAndHeatConductance/VaporConductance;
    this->VPD = Es(Tair)*(1-RH); // vapor pressure deficit Es is saturation vapor pressure at air temperature
    // iterative version
    newTi=Tleaf*0.7;
    iter=0;
    lastTi=Tleaf;
    double Res, dRes; //temporary variables
    double thermal_leaf;
    Ea = Es(Tair)*RH; // ambient vapor pressure
    while ((fabs(lastTi-newTi)>0.001) && (iter <maxiter))
    {
        lastTi=newTi;
        //double Tleaf2= Tair + (R_abs- thermal_air-Lambda*VaporConductance*this->VPD/Press)/(Cp*RadiativeAndHeatConductance+Lambda*Slope(Tair)*VaporConductance); // eqn 14.6a
        thermal_leaf=epsilon*sbc*pow(lastTi+273,4)*2;
        Res = R_abs - thermal_leaf - 2*Cp*HeatConductance*(lastTi - Tair) - Lambda*VaporConductance*(Es(lastTi)-Ea)/Press; // Residual function: f(Ti), KT Paw (1987)
        dRes= -4*epsilon*sbc*pow(273+lastTi,3)*2-2*Cp*HeatConductance-Lambda*VaporConductance*Slope(lastTi); // derivative of residual: f'(Ti)
        newTi = lastTi - Res/dRes; // newton-rhapson iteration
        iter++;
        //printf("iter %d newTi %f lastTi %f Res %f dRes %f\n",iter,newTi,lastTi,Res,dRes);
    }
    Tleaf=newTi;
    //printf("test energy Ta %f Tl %f s %f\n",Tair,Tleaf,Slope(Tair));
    //exit(0);
    Transpiration =VaporConductance*(Es(Tleaf)-Ea)/Press; //Don't need Lambda - cancels out see eq 14.10 in Campbell and Norman, 1998
    // umol m-2 s-1. note 1000 converts from moles to umol since units of VaporConductance are moles.
    //printf("ml marks flag Tr %f %f %f %f %f %f %f\n",Transpiration,VaporConductance,StomatalConductance,BoundaryLayerConductance,Es(Tleaf),Ea,Tleaf);
    //printf("Check Ta %f Tl %f gsw %f gbw %f gH %f Es %f R_abs %f thermal_leaf %f LamdaE %f Res %f dRes %f\n",Tair,Tleaf,StomatalConductance,BoundaryLayerConductance,HeatConductance,Es(Tair),R_abs,thermal_leaf,Lambda*VaporConductance*1.0*(Es(Tleaf)-Ea)/Press,Res,dRes);
   //exit(0);
}

//void SiteData::EnergyBalance2()
//{
//    const long Lambda = 44000; //latent heat of vaporization of water J mol-1 - not used in this implementation
//    const double Cp = 29.3; // thermodynamic psychrometer constant and specific hear of air, J mol-1 C-1
//    const double psc = 6.66e-4; //psycrometric constant units are C-1
//    //psc=Cp/Lambda = 29.3/44000 See Campbell and Norman, pg 232, after eq 14.11
//
//    //The following are secondary variables used in the energy balance
//    double HeatConductance,  //heat conductance J m-2 s-1
//    VaporConductance, //vapor conductance ratio of stomatal and heat conductance mol m-2 s-1
//    RadiativeConductance, //radiative conductance J m-2 s-1
//    RadiativeAndHeatConductance, //radiative+heat conductance
//    psc1,  // apparent psychrometer constant Campbell and Norman, page 232 after eq 14.11
//    Ea,   //ambient vapor pressure kPa
//    thermal_air; // emitted thermal radiation Watts  m-2
//    double lastTi, newTi;
//    int    iter;
//
//    HeatConductance = BoundaryLayerConductance*(0.135/0.147);  // heat conductance, HeatConductance = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Boundary Layer Conductance to Heat
//    // Since BoundaryLayerConductance is .147*sqrt(u/d) this scales to 0.135*sqrt(u/d) - HeatConductance on page 109 of Campbell and Norman, 1998
//    // Wind was accounted for in BoundaryLayerConductance already  as BoundaryLayerConductance (turbulent vapor transfer) was calculated from CalcTurbulentVaporConductance() in GasEx.
//    // units are J m-2 s-1
//    VaporConductance = StomatalConductance*BoundaryLayerConductance/(StomatalConductance+BoundaryLayerConductance);      //vapor conductance, StomatalConductance is stomatal conductance and is given as gvs in Campbell and Norman.
//    // note units are moles m-2 s-1.
//    RadiativeConductance = (4*epsilon*sbc*pow(273+Tair,3)/Cp)*2; // radiative conductance, *2 account for both sides
//    RadiativeAndHeatConductance = HeatConductance + RadiativeConductance;
//    thermal_air = epsilon*sbc*pow(Tair+273,4)*2; //Multiply by 2 for both surfaces
//    psc1 = psc*RadiativeAndHeatConductance/VaporConductance;
//    this->VPD = Es(Tair)*(1-RH); // vapor pressure deficit Es is saturation vapor pressure at air temperature
//    // iterative version
//    newTi=-10;
//    iter=0;
//    lastTi=Tleaf;
//    double Res, dRes; //temporary variables
//    double thermal_leaf;
//    Ea = Es(Tair)*RH; // ambient vapor pressure
//    while ((fabs(lastTi-newTi)>0.001) && (iter <maxiter))
//    {
//        lastTi=newTi;
//        Tleaf= Tair + (R_abs- thermal_air-Lambda*VaporConductance*this->VPD/Press)/(Cp*RadiativeAndHeatConductance+Lambda*Slope(Tair)*VaporConductance); // eqn 14.6a
//        thermal_leaf=epsilon*sbc*pow(Tleaf+273,4)*2;
//        Res = R_abs - thermal_leaf - Cp*HeatConductance*(Tleaf - Tair) - Lambda*VaporConductance*0.5*(Es(Tleaf)-Ea)/Press; // Residual function: f(Ti), KT Paw (1987)
//        dRes= -4*epsilon*sbc*pow(273+Tleaf,3)*2-Cp*HeatConductance*Tleaf-Lambda*VaporConductance*Slope(Tleaf); // derivative of residual: f'(Ti)
//        newTi = Tleaf + Res/dRes; // newton-rhapson iteration
//        iter++;
//        printf("iter %d newTi %f Res %f dRes %f\n",iter,newTi,Res,dRes);
//    }
//    Tleaf=newTi;
//    printf("test energy Ta %f Tl %f s %f\n",Tair,Tleaf,Slope(Tair));
//    Transpiration =VaporConductance*(Es(Tleaf)-Ea)/Press; //Don't need Lambda - cancels out see eq 14.10 in Campbell and Norman, 1998
//    // umol m-2 s-1. note 1000 converts from moles to umol since units of VaporConductance are moles.
//    //printf("ml marks flag Tr %f %f %f %f %f %f %f\n",Transpiration,VaporConductance,StomatalConductance,BoundaryLayerConductance,Es(Tleaf),Ea,Tleaf);
//}



double SiteData::CalcStomatalConductance()
{
    //** \code
    double Ds, //! \b Ds, VPD at leaf surface
    aa,    //! \b aa, a value in quadratic equation
    bb,    //! \b bb, b value in quadratic equation
    cc,    //! \b cc, calcuation variable (x) in quadratic equation
    hs,    //! \b hs, solution for relative humidity
    Cs,    //! \b Cs, estimate of mole fraction of CO2 at the leaf surface
    Gamma, //! \b Gamma, CO2 compensation point in the absence of mitochondirial respiration, in ubar
    StomatalConductance;    //! \b StomatalConductance, temporary variable to hold stomatal conductance
    Gamma = 10.0;
    //** \endcode
    
    double P=Press/100;
    Cs = (CO2 - (1.37*AssimilationNet/BoundaryLayerConductance))*P; // surface CO2 in mole fraction
    if (Cs == Gamma) Cs = Gamma + 1;
    if (Cs <= Gamma) Cs = Gamma + 1;
    // Quadratic equation to obtain hs by combining StomatalConductance with diffusion equation
    aa = g1*AssimilationNet/Cs;
    bb = g0+BoundaryLayerConductance-(g1*AssimilationNet/Cs);
    cc = (-RH*BoundaryLayerConductance)-g0;
    hs = QuadSolnUpper(aa,bb,cc);
    if (hs > 1) hs = 1;
    if (hs<0) hs = 0;
    //Ds = (1-hs)*Es(Tleaf); // VPD at leaf surface
    //StomatalConductance = (g0+g1*(AssimilationNet*hs/Cs));
    //if (StomatalConductance < g0) StomatalConductance=g0; //Limit StomatalConductance to mesophyll conductance
    
    
    Ds=(Es(Tleaf)-RH*Es(Tair))/Press;
    
    StomatalConductance=(0.01+8*AssimilationNet/((Cs-Gamma)*(1+Ds/0.01)));
    if (StomatalConductance < g0) StomatalConductance=g0; //Limit StomatalConductance to mesophyll conductance
    
    //printf("Ds %f An %f Cs %f CO2 %f Gamma %f gsc %f\n",Ds,AssimilationNet,Cs,CO2,Gamma,StomatalConductance);
    
    return StomatalConductance;  // moles H2O m-2 s-1
}

double SiteData::CalcTurbulentVaporConductance(void)
{
    double ratio; /*!< temporary holding variable for stomatal ratio calculations*/
    double Char_Dim; /*!<  characteristic dimension of leaf */
    ratio = Square(stomaRatio+1)/(Square(stomaRatio)+1);
    Char_Dim = LfWidth*0.72; // characteristic dimension of a leaf, leaf width in m
    // wind is in m per second
    //return (1.4*0.147*sqrt(fmax(0.1,wind)/Char_Dim))*ratio;
    return (1.4*0.147*sqrt(fmax(0.1,wind)/Char_Dim));
    // multiply by 1.4 for outdoor condition, Campbell and Norman (1998), p109, gva
    // multiply by ratio to get the effective blc (per projected area basis), licor 6400 manual p 1-9
}

double SiteData::Es(double Temperature)
{
    double result;
    // a=0.611 kPa, b=17.502 C and c=240.97 C
    //Units of Es are kPa
    result=(0.611*exp(17.502*Temperature/(240.97+Temperature)));
    return result;
}

double SiteData::Slope(double Temperature)
{
    double VPSlope;
    // units of b and c are  degrees C
    const double b= 17.502; const double c= 240.97;
    VPSlope=(Es(Temperature)*(b*c)/Square(c+Temperature)/Press);
    return VPSlope;
}

double SiteData::SearchCi(double CO2i)
{
    bool isCiConverged = true;
    int iter;
    double fprime, Ci1, Ci2, Ci_low, Ci_hi, Ci_m;
    double temp;
    Ci1 = CO2i;
    Ci2 = CO2i + 1.0;
    Ci_m = (Ci1+Ci2)/2.0;
    
    iter_Ci = 0;
    iter = 0;

    
    do
    {
        iter++;
        //Secant search method
        if (fabs(Ci1-Ci2) <= errTolerance) {break;}
        if (iter >= maxiter)
        {
            isCiConverged = false;
            break;
        }
        fprime = (EvalCi(Ci2)-EvalCi(Ci1))/(Ci2-Ci1);  // f'(Ci)
        if (fprime != 0.0)
        {
            Ci_m = fmax(errTolerance, Ci1-EvalCi(Ci1)/fprime);
        }
        else
            Ci_m = Ci1;
        Ci1 = Ci2;
        Ci2 = Ci_m;
        temp=EvalCi(Ci_m);
        double temp2=maxiter;
    } while ((fabs(EvalCi(Ci_m)) >= errTolerance) || (iter < maxiter));
    
    
    
    
    // C4 photosynthesis fails to converge at low soil water potentials using secant search, 6/8/05 SK
    // Bisectional type search is slower but more secure
    //Bisectional search
    if (iter > maxiter)
    {
        Ci_low = 0.0;
        Ci_hi = 2.0*CO2;
        isCiConverged = false;
        
        while (fabs(Ci_hi-Ci_low) <= errTolerance || iter > (maxiter*2))
        {
            Ci_m = (Ci_low + Ci_hi)/2;
            if (fabs(EvalCi(Ci_low)*EvalCi(Ci_m)) <= eqlTolerance) break;
            else if (EvalCi(Ci_low)*EvalCi(Ci_m) < 0.0) {Ci_hi = fmax(Ci_m, errTolerance);}
            else if (EvalCi(Ci_m)*EvalCi(Ci_hi) < 0.0)  {Ci_low = fmax(Ci_m, errTolerance);}
            else {isCiConverged = false; break;}
        }
        
    }
    
    CO2i = Ci_m;
    Ci_Ca = CO2i/CO2;
    iter_Ci = iter_Ci + iter;
    iter_total = iter_total + iter;
    return CO2i;
    
}
double SiteData::EvalCi(double Ci)
{
    double newCi;
    
    if (C4== 0) PhotosynthesisC3(Ci);
    if (C4== 1) PhotosynthesisC4(Ci);
    if (fabs(StomatalConductance) > eqlTolerance)
    {
        newCi = fmax(1.0,CO2 - AssimilationNet*(1.6/StomatalConductance+1.37/BoundaryLayerConductance)*(Press/100.0));
    }
    else
        newCi = fmax(1.0,CO2 - AssimilationNet*(1.6/eqlTolerance+1.37/BoundaryLayerConductance)*(Press/100.0));
    return (newCi-Ci);
}

//These two functions solve the quadratic equation.
double SiteData::QuadSolnUpper (double a, double b, double c )
{
    
    /** solves the uppper part of the quadratic equation ax2+bx2=c
     
     @param[in] a
     @param[in] b
     @param[in] c
     
     \return lower portion of x        */
    if (a==0) return 0;
    else if ((b*b - 4*a*c) < 0) return -b/a;   //imaginary roots
    else  return (-b+sqrt(b*b-4*a*c))/(2*a);
}

double SiteData::QuadSolnLower (double a, double b, double c )
{
    /** solves the lower part of the quadratic equation ax2+bx=c
     
     @param[in] a
     @param[in] b
     @param[in] c
     \return lower portion of x
     */
    if (a==0) return 0;
    else if ((b*b - 4*a*c) < 0) return -b/a;   //imaginary roots
    else  return (-b-sqrt(b*b-4*a*c))/(2*a);
}

//*! hyperbolic min
double SiteData::minh(double fn1,double fn2,double theta2)
{
    
    /**
     @param [in] fn1 first value to be compared for min
     @param [in] fn2 second value to be compared for min
     @param [in] theta2  curvature factor
     
     \return hyperbolic minimum
     */
    double x, res;
    
    x = ((fn1+fn2)*(fn1+fn2)-4*theta2*fn1*fn2);
    if (x<0)
    {
        res = fmin(fn1,fn2);
        return res;
    }
    if (theta2==0.0)
    {
        res= fn1*fn2/(fn1+fn2);
        return res;
    }
    else
    {
        res = ((fn1+ fn2) - sqrt(x))/(2*theta2); // hyperbolic minimum
        return res;
    }
}

double SiteData::Square(double a)
{
    return a * a;
    
} /*!< Squares number */

double SiteData::Min(double a, double b, double c)
{
    return (fmin(fmin(a,b),c));
    
} /*!< Finds minimum of three numbers */
