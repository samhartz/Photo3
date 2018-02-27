from dictionaries import *
from math import exp, pi, sqrt, log

def gaCalc(h):
    KVC = 0.41 # von Karman constant
    UZ = 2.  # windspeed (m/s)
    ZW = 2. # height at which windspeed is measured (m)
    """calculate atmospheric conductance (mm/s) from plant height (m)"""
    return 1000.*KVC**2.*UZ/log((ZW-0.64*h)/(0.13*h))**2.
def VPD(ta, qa):
    """Vapor pressure deficit (Pa)"""
    return esat(ta) - (qa*P_ATM)/.622
def drh(ta, qa):
    """Specific humidity deficit (kg/kg)"""
    return (esat(ta)*.622)/P_ATM - qa 
def psi_atm(ta, qa):
    """Atmospheric water potential (MPa)"""
    return R*ta/VW*log(qa*P_ATM/.622/esat(ta))/1000000.
def qi(tl, psi_l):
    """Specific humidity internal to leaf (kg/kg)"""
    return .622*esat(tl)/P_ATM*exp(psi_l*1000000.*VW/R/tl)
def hgtNew(phi, shf, hgt):
    """Boundary layer height (m)"""
    if phi > 0.:
       return ((1. + 2.*B)*shf*1000.)/(RHO_A*CP_A*hgt*GAMMA_THETA)*dt + hgt
    else:
       return ho
def theta_aNew(phi, theta_a, shf, hgt_old, hgt_new):
    """Change in temperature"""
    if phi > 0.:
        return theta_a + dt*shf/(RHO_A*CP_A*hgt_old) + (theta(hgt_old) - theta_a)/hgt_old*(hgt_new - hgt_old)
    else:
        return theta_c(tio, ho)#K
def qaNew(t): 
    """Specific Humidity (kg/kg)"""
    if phi(tch(t)) > 0.:
        return  qa[t] + (RHO_W*Ev[t]/10.**6*dt)/(RHO_A*hgt[t]) + (q(hgt[t]) - qa[t])/hgt[t]*(hgt[t + 1] - hgt[t])
    else:
        return qio
def shfNew(species, tl, ta, lai):
    """Sensible heat flux (W/m^2), per unit ground area"""
    return CP_A*RHO_A*ga[species]*(tl-ta)/1000.*lai

#Soil and plant water fluxes
def leak(sType, s):
    """Leakage (um/s) """
    return .11574*Ks[sType]*s**(2.*b[sType] + 3.)                                   
def psi_s(sType, s):
    """Soil Potential (MPa)"""
    return psi_ss[sType]*(s**-b[sType])  
def RAI(species, s):
    """Root area index (-)"""
    return RAIW[species]*s**-A_ROOT
def gsr(species, sType, s, zr):
    """Soil-Root Conductance, per unit ground area (um/(s-MPa))"""
    return (leak(sType, s)*sqrt(RAI(species, s))*1000000.)/(float(pi)*g*RHO_W*zr)
def gp(species, psi_l):
    """Plant conductance, per unit leaf area (um/(s-MPa))"""
    if psi_l<-10:
        return 0.
    else:
        return gpmax[species]*exp(-(-psi_l/2.)**2.)
def gsrp(species, sType, s, gp, lai, zr):
    """Soil-Root-Plant Conductance, per unit ground area (um/(s-MPa))"""
    return (lai*gsr(species, sType, s, zr)*gp)/(gsr(species, sType, s, zr) + lai*gp)
def gsrfp(species, sType, s, gp, lai, zr):
    """Soil-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
    return (lai*gsr(species, sType, s, zr)*gp/F_CAP)/(gsr(species, sType, s, zr) +  lai*gp/F_CAP) 
def gw(species, vw):
    """Xylem-storage conductance, per unit leaf area (um/(MPa-s))"""
    return gwmax[species]*(vw/vwt[species])**4. 
def evap(sType, s): 
    """Soil evaporation rate, per unit ground area (mm/day)"""
    if s > sh[sType]:
        return EVMAX*(s - sh[sType])/(1. - sh[sType])
    else:
        return 0.
def sNew(sType, species, ev, s, zr, **kwargs): 
    """Soil moisture"""
    if 'vw' in kwargs: ## CapOn==True
        vw = kwargs['vw']
        gp = kwargs['gp']
        psi_l = kwargs['psi_l']
        return (dt/(n[sType]*zr*10.**6)*(-ev +lai*gw(species, vw)*((psi_wf(species, vw))-ev*(1.-F_CAP)/(lai*gp)-psi_l)\
         - (evap(sType, s)*1000.)/(24.*60*60)- leak(sType, s))) + s
    else:
        return (dt/(n[sType]*zr*10.**6)*(-ev - (evap(sType, s)*1000.)/(24.*60*60)- leak(sType, s))) + s

#DEFINITIONS FOR HYDRAULIC CAPACITANCE
def psi_wf(species, vw): 
    """Water potential of stored water (MPa)"""
    return (1./cap[species])*vw/vwt[species] - (1./cap[species])
def vwf(species, vw, ev, gp, psi_l, lai):
    """Stored water volume, per unit leaf area (m3/m2)"""
    return min(vw - gw(species, vw)*((psi_wf(species, vw)) - (ev*(1. - F_CAP))/(lai*gp) - psi_l)*dt/10.**6, vwt[species])
def psi_xf(ev, gp, psi_l):
    """Water potential at connection node x (MPa)"""
    return ev*(1. - F_CAP)/(lai*gp) + psi_l
def qwf(vw, ev, gp, psi_l, lai):
    """Stored water flux, per unit ground area"""
    return (vw - vwf(species, vw, ev, gp, psi_l, lai))*lai*10.**6/dt
def qsf(vw, ev, gp, psi_l, lai):
    """Soil water flux, per unit ground area"""
    return ev - qwf(vw, ev, gp, psi_l, lai)

#Definitions for photosynthesis    
def v_cmax(species, tl, ared):
    """Maximum carboxylation rate (umol/(m^2s))"""
    return ared*Vcmax0[species]*exp(HaV[species]/(R*TO)*(1. - TO/tl))/(1. + exp((SvC[species]*tl - HdV[species])/(R*tl)))
def gamma(species, tl):
    """CO2 compensation point (umol/mol)"""
    return gamma_o[species]*(1. + GAMMA_1*(tl - TO) + GAMMA_2*(tl - TO)**2.);
def jmax(species, tl):
    """Max. e- transport rate (umol/(m^2s))"""
    return Jmax0[species]*exp(HaJ[species]/(R*TO)*(1. - TO/tl))/(1. + exp((SvQ[species]*tl - HdJ[species])/(R*tl))) 
def j(species, phi, tl):
    """Electron transport rate (umol/(m^2s))"""
    return min((phi*10.**6)/(EP*NA)*KAPPA_2*.5, jmax(species, tl)) 
def jpar(species, phi, tl):
    """Electron transport rate (umol/(m^2s), based off of PAR, not total solar radiatoion)"""
    return min(phi*KAPPA_2, jmax(species, tl)) 
def k_o(species, tl):
    """Michaelis-menten coefficient for O2"""
    return Ko0[species]*exp(Hko[species]/(R*TO)*(1. - TO/tl))
def k_c(species, tl):
    """Michaelis-menten coefficient for CO2"""
    return Kc0[species]*exp(Hkc[species]/(R*TO)*(1. - TO/tl))
def a_c(species, ci, tl, ared):
    """Rubisco-limited photosynthetic rate (umol/(m^2s^1))"""
    return v_cmax(species, tl, ared)*(ci - gamma(species, tl))/(ci + k_c(species, tl)*(1. + (oi[species]*1000.)/k_o(species, tl)))
def a_q(species, phi, ci, tl):
    """Light-limited photosynthetic rate (umol/(m^2s^1))"""
    return (j(species, phi, tl)*(ci - gamma(species, tl)))/(4.*(ci + 2.*gamma(species, tl)))
def a_phiciTl(species, phi, ci, tl, ared):
    """Net photosynthetic demand for CO2 (umol/(m^2s^1))"""
    return max(min(a_c(species, ci, tl, ared), a_q(species, phi, ci, tl)),0)
def a_psilc02(species, psi_l):  
    """Vulnerability curve for water potential (-)"""
    if psi_l < psi_laoO[species] :
        return 0.
    elif psi_laoO[species] <= psi_l <= psi_la1O[species] :
        return (psi_l - psi_laoO[species] )/(psi_la1O[species]  - psi_laoO[species])
    else: 
        return 1.
def an(species, phi, psi_l, tl, ci, ared, **kwargs): 
    """Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
    if 'z' in kwargs: # CAM species
        z = kwargs['z']
        M = kwargs['M']
        return a_sc(species, phi, psi_l, tl, ci, z, M, ared) + a_sv(species, phi, tl, psi_l, z, M) 
    else: # C3 and C4 species
        return a_psilc02(species, psi_l)*a_phiciTl(species, phi, ci, tl, ared)
def a_sc(species, phi, psi_l, tl, ci, z, m, ared):
    """Flux from stomata to Calvin cycle (umol/(m^2s))"""
    return max(0, a_psilc02(species, psi_l)*(a_phiciTl(species, phi, ci, tl, ared) - r_dc(species, phi, tl))*(1. - fC(species, z, m)))
def r_d(species, tl):
    """Dark respiration flux (umol/(m^2s))"""
    return Rd0[species]*exp(HkR[species]/(R*TO)*(1. - TO/tl))
def r_dv(species, phi, tl):
    """Flux of dark respiration to vacuole (umol/(m^2s))"""
    return r_d(species, tl)*exp(-phi)
def r_dc(species, phi, tl):
    """Flux of dark respiration to calvin cycle (umol/(m^2s))"""
    return r_d(species, tl)*(1. - exp(-phi))
def fO(z):
    """Circadian order function (-)"""
    return exp(-(z/MU)**CIRC_3)
def fM(species, z, m, tl):
    """Malic acid storage function"""
    return fO(z)*(Mmax[species]*((TH - tl)/(TH - TW)*(1. - ALPHA_2) + ALPHA_2) - m)/(ALPHA_2*Mmax[species]*((TH - tl)/\
        (TH - TW)*(1. - ALPHA_2) + ALPHA_2) + (Mmax[species]*((TH - tl)/(TH - TW)*(1. - ALPHA_2) + ALPHA_2) - m))
def fC(species, z, m):
    """Carbon circadian control function"""
    return (1. - fO(z))*m/(ALPHA_1*Mmax[species] + m)
def a_sv(species, phi, tl, psi_l, z, m):
    """Flux from stomata to vacuole (umol/(m^2s))"""
    # if phi>0 and M < .005:
    #     return 0.
    # else:
    #     if Mmax[species]*((Th - Tl)/(Th - Tw)*(1. - alpha_2) + alpha_2) > M and (1. - k*(Tl - Topt1)**2.) >0:
    #         return (Ammax[species]*(1. - k*(Tl - Topt1)**2.) - Rdv(phi, Tl))*fM(z, M, Tl)*Apsilc02(psi_l)
    #     else:
    #         return 0.

    if Mmax[species]*((TH - tl)/(TH - TW)*(1. - ALPHA_2) + ALPHA_2) > m and (1. - k*(tl - Topt1)**2.) >0:
        return (Ammax[species]*(1. - k*(tl - Topt1)**2.) - r_dv(species, phi, tl))*fM(species, z, m, tl)*a_psilc02(species, psi_l)
    else:
        return 0.
def a_vc(species, phi, cc, tl, z, m):
    """Flux from vacuole to calvin cycle (umol/(m^2s))"""
    return (a_phiciTl(species, phi, cc, tl, 1.) - r_dc(species, phi, tl))*fC(species, z, m)
def MT(species, z, m, tl, phi): 
    """Malic acid equilibrium value"""
    if phi>0.:
        return Mmax[species]*(CIRC_1*((TH - tl)/(TH - TW) + 1.)*(BETA*(z - MU))**3. - BETA*(TH - tl)/(TH - TW)*(z - MU) + \
            CIRC_2*(TH - tl)/(TH - TW) -(1- fO(z))*(1-m/(m+ALPHA_1*Mmax[species]))) 
    else:
        return Mmax[species]*(CIRC_1*((TH - tl)/(TH - TW) + 1.)*(BETA*(z - MU))**3. - BETA*(TH - tl)/(TH - TW)*(z - MU) + \
            CIRC_2*(TH - tl)/(TH - TW)+ (1-fO(z))) 
def zNew(species, phi, m, z, tl):
    return max(0, dt*(m - MT(species, z, m, tl, phi))/(Mmax[species]*60.*tr) + z)
def fD(vpd):
    """Stomatal response to vapor pressure deficit (-)"""
    return 3/13./sqrt(vpd/1000.)
def gsN(species, phi, ta, psi_l, qa, tl, ci, ared, **kwargs):
    """Stomatal conductance to CO2, per unit leaf area (mol/m2/s)"""
    if an(species, phi, psi_l, tl, ci, ared, **kwargs) < 0.:
        return 0.
    else:
        return a1new[species]*an(species, phi, psi_l, tl, ci, ared, **kwargs)/ca*fD(VPD(ta, qa))
def gwN(species, phi, ta, psi_l, qa, tl, ci, ared, **kwargs): 
    """Stomatal conductance to water, per unit leaf area (mol/m2/sec)"""
    #return gsN(phi, Ta, psi_l, qa, Tl, ci, t)*(1.6*(1. + gmgsratio[pType[species]]))/(1.6 + gmgsratio[pType[species]]) + (gcut[species]*po/(1000.*R*Ta))
    return gsN(species, phi, ta, psi_l, qa, tl, ci, ared, **kwargs)*1.6 + (gcut[species]*P_ATM/(1000.*R*ta))
def evf(species, phi, ta, psi_l, qa, tl, ci, lai, ared, **kwargs):
    """Transpiration, per unit ground area (um/sec)"""
    return lai*(1./(gwN(species, phi, ta, psi_l, qa, tl, ci, ared, **kwargs)*R*ta/P_ATM*1000000.)+1./(ga[species]*1000.))**(-1.)\
    *RHO_A/RHO_W*(qi(tl, psi_l)-qa)
def csNew(species, an):
    """CO2 concentration at leaf surface (ppm)"""
    return ca - an/ga[species]
def ciNew(species, cs, ta, qa):
    """CO2 concentration in mesophyll cytosol (ppm)""" 
    return cs*(1.-1./(a1new[species]*fD(VPD(ta, qa))))  
def cmNew(species, cs, ta, qa):
    """CO2 concentration in mesophyll (ppm)"""
    return ciNew(species, cs, ta, qa)  
#     return CiNew(i) - ((1.+VPD(Ta[i], qa[i])/Dxnew[species])*Cs[i])/(a1new[species]*gmgsratio[pType[species]])   
def ccNew(species, cs, ta, qa, z, m):
    """CO2 concentration in mesophyll cytosol resulting from malic acid decarboxylation (ppm)"""
    return cmNew(species, cs, ta, qa) + fC(species, z, m)*c0
def mNew(species, phi, psi_l, cc, tl, z, m): 
    """Malic acid concentration"""
    return max(((dt/ VCM)*(a_sv(species, phi, tl, psi_l, z, m) - a_vc(species, phi, cc, tl, z, m) + r_dv(species, phi, tl))) + m, 0.)
def gsvNew(species, phi, ta, psi_l, qa, tl, ci, ared, **kwargs):
    """Stomatal conductance to water vapor, per unit leaf area (mm/s)"""
    return (R*ta)/P_ATM*1000.*gwN(species, phi, ta, psi_l, qa, tl, ci, ared, **kwargs)

#DEFINITIONS FOR C4
def v_p(psi_l, ci):
    """CO2 concentrating flux (umol/m2/s)"""
    return min((ci*VPMAX0)/(ci + KP), VPR)
def cbsNew(an, cm, psi_l):
    """CO2 concentration in bundle sheath cell (ppm)"""
    return (v_p(psi_l, cm) - an)/GBS + cm
def taSimple(t):
    return TaNight + (TaDay - TaNight)*phiPar(tch(t))/phi_max
def qaSimple(t):
    return qaNight + (qaDay - qaNight)*phiPar(tch(t))/phi_max
def taStep(t):
    if phiPar(tch(t))>0:
        return TaDay
    else:
        return TaNight
def qaStep(t):
    if phiPar(tch(t))>0:
        return qaDay
    else:
        return qaNight
def qaV(t, vpdiff):
    if phiPar(tch(t))>0:
       return qaDay
    else:
       return (esat(TaNight)-esat(TaDay) +vpdiff + (qaDay*P_ATM)/.622)*0.622/P_ATM
def esat(ta):
    """Saturated vapor pressure (Pa)"""
    return A_SAT*exp((B_SAT*(ta - 273.))/(C_SAT + ta - 273.))
def deltaS(ta):
    return (B_SAT*(C_SAT+ta-273.)-B_SAT*(ta-273.))/(C_SAT+ta-273.)**2.*esat(ta)

#Atmospheric boundary layer
def theta(h):
    """Potential temp at top of boundary layer (K), h input in m"""
    return GAMMA_THETA*h/1000. + 293.6
def q(h):
    """Relative Humidity at Top of Boundary layer (kg/kg), h input in m"""
    return -.00285*h/1000 + .01166
def p(h):
    """Pressure Change with Altitude"""
    return P_ATM*(1. - 2.25577*10**-5*h)**5.25588
def theta_c(ta, h):
    """Potential Temperature Calculation from Atmopheric Temperature"""
    return ta*(P_ATM/p(h))**(R_A/CP_A)
def Temp(theta, h):
    """Temperature Calculation from Potential Temperature (K)"""
    return theta*(p(h)/P_ATM)**(R_A/CP_A)
def phiPar(t):
    phi_max = 500. # Maximum solar radiation (W/m^2)
    delta = 12. # day length (h)
    to = 6. # (h)
    return max((4.*phi_max)/delta**2.*(-(t%24.)**2. + (delta + 2.*to)*(t%24.) - to*(to + delta)), 0.)
def phiStep(t):
    phi_max = 500. # Maximum solar radiation (W/m^2)
    t1=6.;
    t2=18.;
    if 0. <= (t%24.) < t1:
        return 0.
    elif t1 <= (t%24.) <= t2:
        return phi_max
    else: 
        return 0.  

#Time Conversions
def steps(duration, timeStep):
    """Change Duration of Simulation to to number of timesteps according to timestep value"""
    return(duration*24*60)/timeStep
def tcm(t):
    """TimeStep conversion, step counter to minutes"""
    return timestepM*t
def tch(t):
    """TimeStep conversion, step counter to hours"""
    return tcm(t+1.)/60.
def td(days):
    """Number of timesteps (days)"""
    return (days*60.*24.)/timestepM
def fBal(p, sType, species, phi, ta, qa, tl, c1, s, lai, gp, ared, zr, kwargs):
    psi_l, tl =p

    if 'vw' in kwargs: # capacitance is true
        vw = kwargs['vw']
        return (evf(species, phi, ta, psi_l, qa, tl, c1, lai, ared, **kwargs)\
               -(gsrfp(species, sType, s, gp, lai, zr)*(psi_s(sType, s) - psi_l) + lai*gw(species, vw)*(psi_wf(species, vw) - psi_l))/\
               (1. + (gsrfp(species, sType, s, gp, lai, zr)*(1. - F_CAP))/(lai*gp) + (gw(species, vw)*(1. - F_CAP))/gp),\
                phi - shfNew(species, tl, ta, lai) -LAMBDA_W*RHO_W*evf(species, phi, ta, psi_l, qa, tl, c1, lai, ared, **kwargs)/1000000.)
    else:
        if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
            return (phi*lai - shfNew(species, tl, ta, lai) -LAMBDA_W*RHO_W*evf(species, phi, ta, psi_l, qa, tl, c1, lai, ared, **kwargs)/1000000., \
                evf(species, phi, ta, psi_l, qa, tl, c1, lai, ared, **kwargs) - gsrp(species, sType, s, gp, lai, zr)*(psi_s(sType, s) - psi_l))
        else:
            return (phi - shfNew(species, tl, ta, lai) -LAMBDA_W*RHO_W*evf(species, phi, ta, psi_l, qa, tl, c1, lai, ared, **kwargs)/1000000., \
                evf(species, phi, ta, psi_l, qa, tl, c1, lai, ared, **kwargs) - gsrp(species, sType, s, gp, lai, zr)*(psi_s(sType, s) - psi_l))