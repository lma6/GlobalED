#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "edmodels.h"
#include <cmath>
#include "site.h"
#include "patch.h"
#include "read_site_data.h"

#ifdef ED

# define epsilonE
#define ch2o 4.218e+3  //specific heat of liquid water (J deg-1 kg-1)
#define cice 2.106e+3 //specific heat of ice (J deg-1 kg-1)

//// Some variables need to be declared later, later
//
// rhow - density of liquid water (all types) (kg m-3)
// wpud - liquid content of puddles per soil area (kg m-2)
// wipud - ice content of puddles per soil area (kg m-2)
// hsoi(N_SIOL_LAYERS) - soil thickness for each layers
// hydraul(N_SIOL_LAYERS) - saturated hydraulic conductivity (m/s)
// tsoil(N_SIOL_LAYERS) - soil temperature for each layers (Cellius degree)
// tsno(N_SIOL_LAYERS) - snow temperature for each layers (Cellius degree)
// wsoi(N_SIOL_LAYERS) - fraction of soil pore space containing liquid water
// wisoi(N_SIOL_LAYERS) - fraction of soil pore space containing ice
// hvasug - latent heaf of vaporization and sublimation for soil surface (J Kg-1)
// hvasui - latent heaf of vaporization and sublimation for snow surface (J Kg-1)
// gadjust - h2o flux due to adjustments in function wadjust (kg h2o m-2 s-1)
// qglif(4) - fraction of soil evaporation from soil liquid, soil ice, puddle liquid and puddle ice respectively
// porosflo(N_SIOL_LAYERS) - porosity after reduction by ice content
// csoi(N_SIOL_LAYERS) - specific heat of soil, no pore spaces (J kg-1 deg-1)
// rhosoi(N_SIOL_LAYERS) - soil density (without pores, not bulk) (kg m-3)


// This function set soil characteristics and soil evaporation fraction from 4 parts and as well as latent heat
void site:ini_soil(site *currents, UserData* data)
{
    double fsand, fclay, fsilt;
    double powliq, powice;
    double dry_soil_conc;
    int isoil;
    
    // zwpmax - assumed maximum value of zwpud&
    // zwpud - fraction of soil surface covered by puddle&
    // zwsoi - volumetric water content of top soil layer &
    //
    double zwpmax=0.5,zwpud=0,zwsoi=0,rwork1=0,rwork2=0,zvap=0,zsub=0;
    
    //Extract percentage of three components
    
    
    for (isoil=0;isoil<N_SIOL_LAYERS;isoil++)
    {
        fsand=currents->sdata->fsand;
        fclay=currents->sdata->fclay;
        fsilt=1.0-fsand-fclay;
        
        powliq=currents->sdata->poros[isoil]*currents->soil_water[isoil]*(1-currents->soil_ice[isoil]);
        powice=currents->sdata->poros[isoil]*currents->soil_ice[isoil];
        
        //dry_soil_conc - dry-soil conductivity, zcondry in IBIS2
        dry_soil_conc=fsand*0.3+fsilt*0.265+fclay*0.250;
        
        consoil(isoil)=dry_soil_conc*pow(0.56*100.0,powliq)*pow(2.24*100.0,powice);
    }
    
    // zwpud - fraction of surface area covered by puddle (range: 0 - zwpmax)
    // zwpmax - maximum value of zwpud (currently assumed to 0.5, same value in IBIS2)
    // 1-zwpud - fraction of surface area covered by soil (range: (1-zwpmax)-1.0 )
    // zwsoi: volumetric water content of top soil layer (range: 0 - 1.0)
    // qglif(0): fraction of soil evaporation (fvapg) from soil liquid
    // qglif(1): fraction of soil evaporation (fvapg) from ice
    // qglif(2): fraction of soil evaporation (fvapg) from puddle liquid
    // qglif(3): fraction of soil evaporation (fvapg) from puddle ice
    
    zwpud=max(0.0,min(zwpmax,zwpmax*(wpud+wipud)/wpudmax));
    zwsoi = (1.0-wisoi(1)*wsoi(0)+wisoi(0));
    
    if (zwsoi>epsilonE)
    {
        rwork1=1.0/zwsoi;
        if (zwpud>epsilonE)  // If there is water in puddles
        {
            rwork2=1.0/(wpud+wipud);
            qglif(0)=(1.0-zwpud)*(1-wisoi(0))*wsoi(0)*rwork1;
            qglif(1)=(1-zwpud)*wisoi(0)*rwork1;
            qglif(2)=zwpud*wpud*rwork2;
            qglif(3)=zwpud*wipud*rwork2;
        }else
        {
            qglif(0)=(1-wisoi(0))*wsoi(0)*rwork1;
            qglif(1)=wisoi(0)*rwork1;
            qglif(2)=0.0;  //As no water in puddle liquid and lice
            qglif(3)=0.0;
        }
    }
    else  //for a 100% dry soil surface, assign all soil evaporation to the puddles.
          //Note that for small puddle sizes, this would lead to negative puddle depths.
          //However, for a 100% dry soil with small puddles, evaporation is likely to be very small or less than zero (condensation),
        //so negative puddle depths are not likely to occur; --note from IBIS2.0
    {
        if (zwpud>epsilonE)   //If puddles still have water
        {
            rwork2=1.0/(wpud+wipud);
            qglif(0)=0.0;
            qglif(1)=0.0;
            qglif(2)=zwpud*wpud*rwork2;
            qglif(3)=zwpud*wipud*rwork2;  // wpud - puddle water in liquid?  wipud - ice in puddle?
        }
        else
        {
            if (tsoil(0)>tmelt)  //Above freezing, tmelt is to be defined
            {
                qglif(0)=0.0;
                qglif(1)=0.0;
                qglif(2)=1.0;  //means all soil evaporation is from liquid water in puddle
                qglif(3)=0.0;
            }
            else  //Below freezing
            {
                qglif(0)=0.0;
                qglif(1)=0.0;
                qglif(2)=0.0;
                qglif(3)=1.0; //means all soil evaporation is from ice in puddle
            }
        }
    }
    
    // Set letent heat values
    zvap = hvapf(tsoil(0),ta);   // hvapf is to be defined as well as ta; later
    zsub = hsubf(tsoil(0),ta);   // same as above  later
    
    hvasug=(qglif(0)+qglif(1))*zvap+(qglif(2)+qglif(3))*zsub; // latent heat from soil surface , later
    hvasu=hsubf(tsno(0),ta);
    
}

bool soilControl(site *new_site, UserData* data)
{
    int isoil,k; //loop indices
    // zfrez - factor decreasing runoff fraction for tsoil < tmelt
    // zrunf - fraction of rain that doesn't stay in puddle (runoff fraction)
    // wipre - storing variable
    // zdpud - used to compute transfer from puddle to infiltration
    // cx - average specific heat for soil, water and ice
    // zwsoi - what?
    // owsoi[N_SIOL_LAYERS] - old value for wsoi
    // otsoi[N_SIOL_LAYERS] - old value for tsoi
    // c0pud[N_SIOL_LAYERS] - layer heat capacity due to puddles (=0 except for top)
    // c1pud[N_SIOL_LAYERS] - updated av. specifilayer heat capacity due to puddle
    // wflo[N_SIOL_LAYERS+1] - drainage at the bottom, returned by soilHydrology
    // fwtop - evaporation rate from soil (for soilHydrology)
    // fhtop - heat flux through soil surface (for soilThermal)
    // fwpud - portion of puddle that infiltrates in soil (rate)
    // fsqueez - excess amount of water (soilHydrology)
    // dh - correction if wate at tsoi < tmelt or ice at temp > tmelt
    // dw - what ?
    // zporos - what ?
    double zfrez, zrunf, rwork, wipre, zdpud, cx, chav, zwsoi;
    double owsoi[N_SIOL_LAYERS], otsoi[N_SIOL_LAYERS], c0pud[N_SIOL_LAYERS], c1pud[N_SIOL_LAYERS], wflo[N_SIOL_LAYERS+1];
    double fwtop, fhtop, fwpud, fsqueez, dh, dw, zporos;
    
    // Do we need to the same thing? later
    //call const (c0pud, npoi * nsoilay, 0.0)
    //call const (c1pud, npoi * nsoilay, 0.0)
    
    // For soil, set soil infiltration rate fwtop (for soilHydrology) and upper heat flux fhtop (for soil Thermal)
    //
    // also step puddle model wpud, wipud
    //
    // procedure is:
    //
    // 1. Immediately transfer any excess puddle liquid to runoff.
    // 2. Apportion raing (raining?) between puddle liquid (wpud) or runoff (grunof)
    // 3. Apportion evaporation/condensation (fvapg) between infiltration rate (fwtop), soil ice (wisol[]), puddle liquid (wpud) or puddle ice (wipud)
    // 4. Transfer some puddle liquid to fwtop
    // 5. Compute upper heat flux fhtop: include fwtop*ch2o*tsoi(0) to be consistent with whflo in soilThermal, and accounts for change rain temperature from tsoi(0) and runoff temperature from tsoi to max(tsoi(0),tmelt)
    
    
    //Step 1, immediately transfer any excess puddle liquid to runoff, the following runoff formulation could
    // give rise to very small amounts of negative runoff --note from IBIS2
    grunof = min(wpud,max(0.0,wpud+wipud-wpudmax))/dtime; // Transfer any excess liquid to runoff, dtime is to be defined, later
    wpud = wpud -grunof*dtime; //dtime is to be defined, later
    
    // Step 2. Apportion sfc-level rain between puddle liquid and runoff
    
    zfrez = max(0.0, min(1.0, (tsoi(0)-tmelt+0.5)*0.5));  // is 0 for tsoil(0)<tmelt; otherwise, will be 1
    zrunf = zfrez*max(0.0,min(1.0,(wpud+wipud)/wpudmax));  // is 0 for below freezing. is (wpud(i) + wipud(i)/wpudmax when puddle water does not reach maximum. Otherwise, it is 1
    
    //Add part of rainfall to puddles
    wpud += (1.0-zrunf)(raing*dtime);  // raing is precipitation rate? later
    
    //Add rest of rainfall to runoff
    grunof += zrunf*raing;
    
    rwork = fvapg*dtime;
    
    if (fvapg>=0)  //If soil evaporation is positive, reduce water from 4 water stores
    {
        fwtop -= qglif(0)*fvapg;  //Why here use fvapg instead of rwork, later
        wpud -=qglif(2)*rwork;
        wipud -= qglif(3)*rwork;
        
        wipre = wisoi(0);
        wisoi = max(0.0, wipre-qglif(i,1)*rwork/(rhow*poros(0)*hsoi(0)));  // I think this line estimate whether evaporation amount from ice will exceed avaliable ice content, if exceed, all ice are removed for fvapg and set wisoi as 0
        
        if (1.0-wisoi(0)>epsilonE)
            wsoi(0) = wsoi(0)*(1-wipre)/(1-wisoi(0)); //Update wsoi as wisoi has changed
    }
    else
    {
        // Condensation: give all to puddles (to avoid wsoi, wisol>1), why?
        fwtop = 0.0;
        wpud = wpud - (qglif(0)+qglif(2))*rwork;
        wipud = wipud - (qglif(1)+qglif(3))*rwork;
    }
    
    
    // Step 2, transfer some puddle liquid to infultration; can lead to small amounts of negetive wpud (in soilHydrology) due to round-off error
    
    zdpud = rhow*dtime * pow(max(0.0,1-wisoi(0)),2.0) * hydraul(0);
    
    fwpud = max(0.0, min(wpud,zdpud))/dtime;
    
    c0pud(0) = ch2o * wpud + cice * wipud;
    
    // Step 3 compute upper soil heat flux
    // traing - rainfall temperature at soil level (cellius degree), later
    // soil heat are detetmined by ground heaf flux, heat from rainfall and loss heat due to runoff
    fhtop = heatg + raing * ch2o * (traing - tsoi(0)) - grunof * ch2o * max(tmelt - tsoi(0),0.0);
    
    gadjust = 0.0;
    
    
    // Reduce soil moisture due to transpiration (upsoiu & upsoil) from turvap. need to do that before other time
    // stepping below since specific heat of this transport is neglected
    //
    // first set porosflo, reeduce porosity due to ice content, used as the effective porosity for uptake here and liquid hydraulics
    // later in soilHydrology. To avoid divide-by-zeros, use small epsilon limit; This will always cancel with epsilon or 0 in
    // numerators.
    //
    // Also increment soil temperature to balance transpired water differential between temperature of soil and leaf, Physically should
    // apply this to the tree, but would be awkard in turvap.
    // Also, ave old soil moisture and temperature otsoi so implicit soil2Hrdrylogy can aposterioir deduce fluxes.
    //        --- Notes from IBIS2.0
    
    for (isoil=0;isoil<N_SIOL_LAYERS;isoil++)
    {
        porosflo(isoil) = poros(isoil) * max(epsilonE,(1.0 - wisoi(isoil))); // porosity for liquid water, remove some for ice
        
        porosflo(isoil) = max(porosflo(isoil), epsilonE);
        
        // psoiu(isoil), upsoil(isoil) are transpiration rate from upper and lower canopies at each soil layer
        // reduce wosil due to canopy transpiration
        wsoi(isoil) = wsoi(isoil) - dtime *(upsoiu(isoil)+upsoil(isoil))/(rhow*porosflo(isoil)*hsoi(isoil));
        
        // Calculate average specific heat over three components in soil (i.e. water, ice and non-pore) using porosity[kg m-3]*density[kg m-3]*volumetic specific heat[J kg-1 deg-1], 1st term is heat in puddles, 2nd term is heat in non-pores, 3rd heat in liquid water, 4th term is heat in ice
        cx = c0pud(isoil)
           + ((1-poros(isoil)) * csoi(isoil) * rhosol(isoil)
           + poros(isoil) * (1.0-wisoi(isoil)) * wsoi(isoil) * ch2o * rhow
           + poros(isoil) * wisoi(isoil) * cice*rhow) * hsoi(isoil);
        
        //Update soil temperature of each layer by considering heat loss for transpiration at upper and lower canopies
        tsoi(isoil) = tsoi(isoil)-dtime * ch2o * (upsoiu(isoil) * (tu - tsoi(isoil)) + upsoil(isoil) * (tl - tsoi(isoil)))/cx; //why devide by cx rather than times cx
        
        owsoi(isoil) = wsoi;
        otsoi(isoil) = tsoi;
    }
    soilHydrology(owsoi[N_SIOL_LAYERS], fwtop, fwpud, fsqueez, wflo);
    return 1;
}


bool soilHydrology(double owsoi, double fwtop, double fwpud, double fsqueez, double wflo);
{
    return 1;
    
}

bool soilThermal()
{
    
    return 1;
}

bool wadjust()
{
    
    return 1;
}
