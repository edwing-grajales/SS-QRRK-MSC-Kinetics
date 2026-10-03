#python 3
"""
Created on Tue 05 14:43:02 2020

@author: ejgrajales

v2B3 (2026): v2B2 + the robust integration of F_E and of the Dean/Gilbert
Delta factors ("line B"), and nothing else. Tagged "# v2B3 CHANGE".
  Why: scipy.quad on [0, 20000] kcal/mol misses the integrand at low T,
  because all of it lies within a few RT (~0.4 kcal/mol at 200 K) of E0 or of
  E = 0. F_E then comes out << 1 (impossible) and Delta_N far too small, so
  beta_c -> 1 and the Dean factor is wrong by orders of magnitude. The new
  integrals use panels of RT/4 placed where the integrand is, work in
  logarithms (no overflow for molecules with >= 73 modes) and stop where the
  integrand is < e^-60 of its maximum.
  Everything else (constants, QRRK sums, 200-term loops) is exactly as in v2B2.
  Set ROBUST_INTEGRALS = False to get the v2B2 numbers back.

v2B2 (2026): minimal fixes only; numerically identical to v2B1. Every change
is tagged "# v2B2 FIX":
  1. addsheet() no longer writes into the input workbook. In pandas >= 2.0,
     pd.ExcelWriter(file) empties the file as soon as it is created, and the
     next line (writer.book = book) raises AttributeError. The input is left
     at 0 bytes, and the next run fails with "ValueError: Excel file format
     cannot be determined". The fit tables now go to a separate file,
     Fit_<name>.xlsx.
  2. writer.save() was removed in pandas 2.0; writer.close() saves.
  3-6. NumPy >= 1.23/2.x compatibility (np.asscalar removed; float() and
     scalar assignment of 1-element arrays).
  7. Typo in the overflow check of the F(E) denominator loop (Keq_t -> Keqd_t);
     it never triggered, so the results do not change.
  8. lmfit >= 1.2 rejects scalar data in Model.fit (Redlich-Kwong step).
"""

import os
import sys
import math
from lmfit import Model
from scipy import stats
from scipy.integrate import quad
from scipy.special import gamma# Factorial.
from scipy.special import gammaln                     # v2B3 CHANGE
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt

start_time = time.time()

# File names. This should be specify by the user
file = 'H-Assisted_SM.xlsx'   # Name of the main Excel file (example: H-assisted propen-2-ol, Combust. Flame 2021)

# Importing basic data
data = pd.read_excel(file,'Rates')
rxn = pd.read_excel(file,'Energy')
SS = pd.read_excel(file,'SS-QRRK')

k1 = np.float64(np.array(data['k1'].values.tolist()))        # HPL Rate constants for A + B --> AB
km1 = np.float64(np.array(data['km1'].values.tolist()))      # HPL Rate constants for AB --> A + B
k2 = np.float64(np.array(data['k2'].values.tolist()))        # HPL Rate constants for AB --> P
energym1 = np.array(rxn['energym1'].values.tolist())         # Endo or exothermicity of reaction AB --> A + B
energy2 = np.array(rxn['energy2'].values.tolist())           # Endo or exothermicity of reaction AB --> P
T = np.float64(np.array(data['T'].values.tolist()))          # Temperatures

Fr = np.float64(np.array(SS['Fr'].values.tolist()))          # Frequencies Intermediate AB [cm-1]
Fr = Fr[~np.isnan(Fr)]                                       # Frequencies without NaN elements. Due to different elements with Temperature.
E_Z = np.float64(np.array(SS['E_Z'].values.tolist()))        # Zero point Energy of AB [kcal mol-1 K-1]
E_Z = E_Z[~np.isnan(E_Z)].item()                            # v2B2 FIX 3
alfa_d = np.float64(np.array(SS['alfa_d'].values.tolist()))  # Energy transferred  in the deactivation process [cm^(-1)]
alfa_d = alfa_d[~np.isnan(alfa_d)]   
Td = np.float64(np.array(SS['Td'].values.tolist()))          # Temperature-dependent parameter for the deactivation process
Td = Td[~np.isnan(Td)]  
PaA = np.float64(np.array(SS['PaA'].values.tolist()))        # Diameter[angstrom]-e/k[K]-MW. Lennard-Jones Parameters for AB. Taken from nn-pronane. Prausnitz.
PaA = PaA[~np.isnan(PaA)]  
PaB = np.float64(np.array(SS['PaB'].values.tolist()))        # Diameter[angstrom]-e/K[K]-MW-T_c[K]-P_c[bar]. Lennard-Jones Parameters for Bath Gas. Taken from Nitrogen. Prausnitz.
PaB = PaB[~np.isnan(PaB)]  
P = np.float64(np.array(SS['P'].values.tolist()))            # Test Pressures [atm]
P = P[~np.isnan(P)]
# MaxI = int(SS['MaxI'])                      # Maximum iteration steps. 100 is a good initial value.
# MaxI = MaxI[~np.isnan(MaxI)]
MaxI = 200
ROBUST_INTEGRALS = True    # v2B3 CHANGE (line B). False -> v2B2 results (scipy.quad on [0, 20000])

# Determining the type of reaction
energy= [energym1, energy2]
for energy in energy:
    if energy == 0:
        print('Reaction is endothermic')
    else:
        if energy == 1:
            print('Reaction is exothermic')
        else:
            print('Error specifying the type of reaction (endothermic, enrgy=0 and exothermic=1)')
            print('Program stopped')
            sys.exit()

# Useful values
R = 0.00198588 # [=] kcal/(mol*K)
R_RK = 83.14462#Universal gas constant[cm^3 bar mol^(-1) K^(-1)]
NA = 6.022140857e23#Avogadro number.
kB = 1.38064852e-23#Boltzmann constant [m^2 kg s^-2 K^-1]
kB_2 = 3.29982916e-27#Boltzmann constant [kcal K^-(1)]

#Function definition depending of type of reaction
def fit(x,lnA,n,E,T0,**arb):
    ener = arb['heat']                                            # I could add a 'print' in the following line to check this value                        
    if ener == 0:
        return lnA+n*np.log((x)/300)-E*(x+T0)/(R*(x**2+T0**2))    # Endothermic equation fitting
    else:                                                         # cle energy == 1.0:
        return lnA+n*np.log((x+T0)/300)-E*(x+T0)/(R*(x**2+T0**2)) # Exothermic equation fitting

# Fitting using model
def optim (k,ener):
    arb = dict(heat=ener)
    Amodel = Model(fit)
    result = Amodel.fit(np.log(k),x=T,lnA=20.0,n=1.0,E=30.0,T0=90.0,**arb)
    A = result.params
    P = np.array([np.exp(A["lnA"]._val),A["n"]._val,A["E"]._val,A["T0"]._val])
    ter_1 = P[2]*np.divide(T**4 + 2*P[3]*T**3 - P[3]**2*T**2,(T**2 + P[3]**2)**2)#T**4 was given problem, the digit was to big for int34, then I changed to int64 with np.float64
    ter_2 = P[1]*R*T
    if ener == 0:
        Ea = ter_1 + ter_2              #Endothermic
    else:
        Ea = ter_1 + ter_2*T/(T + P[3]) #Exothermic
    A_inf = k*np.exp(Ea/(R*T))   
    return np.transpose(np.concatenate(([A_inf],[Ea]))), result

# Adding fitting data to the file
def addsheet(result,sheetname):
    names = ['T', 'A_inf', 'Ea']
    results = np.column_stack((T, result))
    df = pd.DataFrame(data=results, columns=names)
    # v2B2 FIX 1: never write into the input workbook `file`.
    out = 'Fit_' + sheetname + '.xlsx'
    if os.path.abspath(out) == os.path.abspath(file):
        raise RuntimeError('Refusing to overwrite the input workbook')
    df.to_excel(out, index=False, sheet_name=sheetname)

# Fitting plots
def plots(k,result,lab):
    fig,axs = plt.subplots(nrows=2,ncols=1,sharex=True,gridspec_kw={'height_ratios':[1,2]})
    axs[0].scatter(T,result.residual,color='royalblue',marker='o',label='Residuals')
    axs[0].axhline(y=0.0, color='royalblue',linestyle='-')
    axs[0].set_ylabel('Residuals')
    axs[0].legend(loc=4)
    axs[0].title.set_text('Fitting of '+lab+' for SS-QRRK Parameters')
    axs[1].plot(T,np.log(k),'ro',label='Theoretical values')
    axs[1].plot(T,result.best_fit,'m-',label='Fitting')
    axs[1].legend()
    axs[1].set_xlabel('Temperature (K)')
    axs[1].set_ylabel(lab)
    axs[1].text((T[0]+T[-1])/2,(np.log(k[0])+np.log(k[-1]))/2, 'reduced chi-square=' + str(result.redchi)[0:6])
    fig.savefig(lab+'.png')
    
# Solving the system
resultm1,datam1 = optim(km1,energym1[0]) 
result2,data2 = optim(k2,energy2[0])
plots(km1,datam1,'ln(k\u2098\u2081)')
plots(k2,data2,'ln(k\u2082)')
addsheet(resultm1,'km1') #Optional
addsheet(result2,'k2')   #Optional


## SS-QRRK Algorithm for Chemical Activated Reactions

# Collision efficiency computation
Te = T; E0 = resultm1[:,1]; E02 = result2[:,1]                       # Varaiables in format SS-QRRK
Am1 = resultm1[:,0]; A2 = result2[:,0]                               # Varaiables in format SS-QRRK
vgmean = stats.gmean(Fr)/349.75                                      # Geometric mean of the frequencies [kcal mol^(-1)]
proFrcal = np.prod(Fr/349.75)                                        # Productory of the frequencies in [kcal mol^(-1)]^n, with "n" the number of data. Used to computed density of states 
alfa = alfa_d/349.75*(np.array(Te)/300)**Td                          # Total energy transferred in the deactivation process [kcal mol^(-1)]
B = (len(Fr)-1)**2*sum((Fr/349.75)**2)/(len(Fr)*(sum(Fr/349.75))**2) # Empirical parameter to obtain F_E

#Numerical Integration
def density(E):
    w = []# Parameter function of the threshold and zero-point vibrational energy needed to finally get the value F_E.
    if E > E_Z:
        w.append(10**(-1.506*(E/E_Z)**0.25))
    else:
        w.append((5*E/E_Z+2.73*(E/E_Z)**0.5+3.51)**-1) 
    a = 1.0-B*w[0]# Correction factor                                   # v2B2 FIX 4: scalar, so quad() receives a float
    # d_states = (E+a*E_Z)**(len(Fr)-1)/(gamma(len(Fr)-1+1)*np.prod(Fr)/349.75**len(Fr))# (Number of states)/ [kcal mol⁻¹]. 
    ter_n = (E+a*E_Z)**(len(Fr)-1)
    ter_d1 = gamma(len(Fr)-1+1)
    ter_d2 = np.prod(Fr)/349.75**len(Fr)
    ter_d = ter_d1*ter_d2
    d_states = ter_n/ter_d
    return d_states

def ddeltaN(E,T):
    return density(E)*np.exp(-E/(R*T))

def ddelta_one(E,T):
    return density(E)*np.exp(-E/(R*T))
    
def ddelta_two(E,T,Er,FET): # Just to distinguish the notation Er=E0 and FET=F_E
    return density(E)*np.exp(-E/(R*T))*np.exp((E-Er)/(FET*R*T))
    

# -----------------------------------------------------------------------------
# v2B3 CHANGE (line B): robust evaluation of the F_E, Delta_N, Delta_1 and
# Delta_2 integrals (see the header). Set ROBUST_INTEGRALS = False above to get
# the v2B2 numbers back.
# -----------------------------------------------------------------------------
GL_X, GL_W = np.polynomial.legendre.leggauss(20)      # 20-point Gauss-Legendre rule

def ln_density(E):
    """ln rho(E): the same Whitten-Rabinovitch formula as density(), but in
    logarithms (no overflow for large molecules) and vectorised over E."""
    x = np.asarray(E, dtype=float)/E_Z
    with np.errstate(divide='ignore', invalid='ignore'):
        w = np.where(x > 1.0, 10.0**(-1.506*np.abs(x)**0.25),
                     1.0/(5.0*x + 2.73*np.sqrt(np.abs(x)) + 3.51))
    a = 1.0 - B*w
    s = len(Fr)
    return (s-1)*np.log(E + a*E_Z) - gammaln(s) - np.sum(np.log(Fr/349.75))

def gl_integrate(f, edges):
    """Composite Gauss-Legendre quadrature of f over consecutive panels."""
    a, b = edges[:-1, None], edges[1:, None]
    half, mid = 0.5*(b - a), 0.5*(a + b)
    return float(np.sum(half*GL_W*f(mid + half*GL_X)))

def panels(lo, hi, h, extra=()):
    """Panel edges no wider than h (= RT/4), refined near E = 0 and with break
    points at E0 and at E_Z (where a(E) changes formula)."""
    n = max(int(math.ceil((hi - lo)/h)), 1)
    e = [np.linspace(lo, hi, n + 1)]
    if lo == 0.0:
        e.append(h*2.0**-np.arange(1, 30))
    e.append(np.array([p for p in extra if lo < p < hi]))
    return np.unique(np.concatenate(e))

def energy_integrals(T, E0, alpha):
    """Return F_E, ln(Delta_N), Delta_1, Delta_2, Delta_Gilbert, Delta_Dean."""
    RT = R*T
    s = len(Fr)
    g = lambda E: ln_density(E) - E/RT                    # ln[rho(E) exp(-E/RT)]
    h = RT/4.0                                            # panels of RT/4: the natural energy scale
    E_hi = E0 + RT*(s + 20.0*math.sqrt(s) + 100.0)        # well past the Boltzmann peak
    for _ in range(10):                                   # extend until the tail is < e^-60
        gmax = float(np.max(g(np.linspace(0.0, E_hi, 2000))))
        if g(E_hi) - gmax < -60.0:
            break
        E_hi *= 2.0
    f = lambda E: np.exp(g(E) - gmax)                     # scaled integrand, maximum = 1
    below = panels(0.0, E0, h, (E0, E_Z))
    I_below = gl_integrate(f, below)                      # int_0^E0
    I_above = gl_integrate(f, panels(E0, E_hi, h, (E0, E_Z)))   # int_E0^inf
    dN = I_below + I_above
    FE = math.exp(math.log(I_above) - (float(g(E0)) - gmax))/RT
    fert = FE*RT
    c = fert/(alpha + fert)
    D1 = I_below/dN
    D2 = gl_integrate(lambda E: f(E)*np.exp(-(E0 - E)/fert), below)/dN
    # Delta_Dean = D1 - c*D2, written without subtracting two numbers close to 1
    DD = gl_integrate(lambda E: f(E)*((1.0 - c) - c*np.expm1(-(E0 - E)/fert)), below)/dN
    x, y = E0/fert, E0/(alpha + fert)
    DG = -math.expm1(-x) - math.exp(-x)*y                 # = 1 - exp(-x)(1 + y)
    return FE, gmax + math.log(dN), D1, D2, DG, DD


deltaN,delta_one,delta_two,DGilbert,DDean,Bc_SS,Bc_Gilbert,Bc_Dean = [],[],[],[],[],[],[],[]
# FE_I, FE_II = [], [], []
FE_test = []

for i in range (len(Te)):
    if ROBUST_INTEGRALS:                                                                                                   # v2B3 CHANGE (line B)
        FE_i, lnDN_i, D1_i, D2_i, DG_i, DD_i = energy_integrals(Te[i], E0[i], alfa[i])
        FE_test.append(FE_i); deltaN.append(lnDN_i)                                                                        # note: deltaN now holds ln(Delta_N)
        delta_one.append(D1_i); delta_two.append(D2_i); DGilbert.append(DG_i); DDean.append(DD_i)
    else:                                                                                                                  # v2B2 integrals, unchanged
        FE_test.append((quad(ddelta_one,E0[i],20000,args=(Te[i]))[0])/(R*Te[i]*density(E0[i])*np.exp(-E0[i]/(R*Te[i]))))  # Normalized number of states above threshol energy numerical.     
        deltaN.append(quad(ddeltaN, 0, 20000, args=(Te[i]))[0])                                                               # Delta_N of Dean et al.
        delta_one.append((quad(ddelta_one, 0, E0[i], args=(Te[i]))[0])/deltaN[i])                                              # Delta_one of Dean et al.
        delta_two.append((quad(ddelta_two, 0, E0[i], args=(Te[i], E0[i], FE_test[i]))[0])/deltaN[i])                           # Delta_two of Dean et al.
        DGilbert.append(1-np.exp(-E0[i]/(FE_test[i]*R*Te[i]))*(1+E0[i]/(alfa[i]+FE_test[i]*R*Te[i])))                         # Complementaryn factor from Gilbert.
        DDean.append(delta_one[i]-delta_two[i]*FE_test[i]*R*Te[i]/(alfa[i]+FE_test[i]*R*Te[i]))                               # Complementary factor from Carstensen.
   
    Bc_SS.append((alfa[i]/(alfa[i]+R*FE_test[i]*Te[i]))**2)                                                               # Collisonal efficiency of SS-QRRK.   
    if FE_test[i] >= 1e6:                                                                                                 # A. Y. Chang, J. W. Bozzelli 2 and A. M. Dean. Zeitschrift für Physikalische Chemie, 214, 11, 153321568 (2000).
        j = i-1 
        Bc_Gilbert.append(Bc_Gilbert[j])
        Bc_Dean.append(Bc_Dean[j])                                                                                        # Complemented collisional efficiency corrected from Dean.   
    else:
        Bc_Gilbert.append(Bc_SS[i]/DGilbert[i])                                                                           # Complemented collisional efficiency from Gilbert.
        Bc_Dean.append(Bc_SS [i]/DDean[ i])                                                                               # Complemented collisional efficiency corrected from Dean.

# Lennard-Jones collision rate Constant r_c=k_c*[M]
mA = PaA[2]/NA/1e3; mB = PaB[2]/NA/1e3                                                                                    # Mass of A particle [kg] and mass of B particle (bath gas) [kg]
k_HS = NA*math.pi*((2**(1/6)*1e-8*(PaA[0]+PaB[0]))/2)**2*1e2*np.sqrt(8*np.array(Te)*kB/(mA*mB/(mA+mB)*math.pi))           # Hard-sphere collision rate constant [cm^3*mol^(-1)*s^(-1)].

ekAB = math.sqrt(PaA[1]*PaB[1])
Ome22 = []
for i in range (len(Te)):
    if Te[i]/ekAB >=3 and Te[i]/ekAB <=300:
        Ome22.append((0.697+0.5185*np.log10(Te[i]/ekAB))**-1)      # Unitless
    else:
        Ome22.append((0.636+0.5670*np.log10(Te[i]/ekAB))**-1)      # Unitless

kc_SS = np.transpose(np.transpose(Bc_SS)*np.array(Ome22)*np.array(k_HS))             # SS-QRRK deactivation rate constant
kc_Gilbert = np.transpose(np.transpose(Bc_Gilbert)*np.array(Ome22)*np.array(k_HS))   # Gilbert deactivation rate constant
kc_Dean = np.transpose(np.transpose(Bc_Dean)*np.array(Ome22)*np.array(k_HS))         # Dean deactivation rate constant

# Ideal gas and initial value for EK
M0,rc_SS = [],[]
for i in range(len(Te)):
    M0.append(np.array(P)/(R*4.1294e4*Te[i]))
    rc_SS.append(kc_SS[i]*np.array(P)/(R*4.1294e4*Te[i]))# [s-1]
M = np.matrix(M0); V = 1/M

# Redlich-Kwong EOS
aa = 0.42747*R_RK**2*PaB[3]**2.5/PaB[4]; bb = 0.08664*R_RK*PaB[3]/PaB[4]; pp = np.array(P)*1.01325;pp = pp  
vol1 = np.zeros((len(Te),len(pp))); test = np.zeros((len(Te),len(pp)))
for q in range(len(pp)):
    vol,vol0,vol00  = [],[],[]
    for k in range(len(Te)):
        def f(x,v):
            return R_RK*x/(v-bb)-aa/(x**0.5*v*(v+bb))    
        Amodel = Model(f)
        fresult = Amodel.fit(np.array([pp[q]]), x=Te[k], v=V[k,q])                        # v2B2 FIX 8: lmfit >= 1.2 needs array data (same fix as unimolecular v4)
        vol.append(fresult.best_values)
        vol0 = pd.DataFrame(vol)
        vol00 = np.array(vol0)
        vol1[k][q] = vol00[k].item()                                                    # v2B2 FIX 5: np.asscalar removed in NumPy 1.23
        test[k][q] = ((R_RK*Te[k]/(vol00[k]-bb)-aa/(Te[k]**0.5*vol00[k]*(vol00[k]+bb)) - pp[q])**2).item()  # v2B2 FIX 6
M_RK = 1/vol1

# Chemical activated rate constants

n = np.array(E0)/vgmean; nm1 = np.array(E0)/vgmean                 # Initial value of n, which is common n all calculations
mm1 = np.array(E0)/vgmean; m2 = np.array(E02)/vgmean               # m quantum number for the reaction -1 and 2, respectively

# denominator of factor F(E)
Fd,Keqd,k2_m1d = [],[],[]
for i in range (len(Te)):
    k2_m1dt,Keqd_t,Fd_t,Keq_t = np.float64(0.0),np.float64(0.0),0.0,0.0
    for j in range(MaxI):
        k2_m1dt = Am1[i]*gamma(nm1[i]+1)*gamma(nm1[i]-mm1[i]+len(Fr)-1+1)/gamma(nm1[i]-mm1[i]+1)/gamma(nm1[i]+len(Fr)-1+1)       # The value of "+1" always stands for the gamma function that is of the way (n+1). This normally goes to 'inf'. "RuntimeWarning: overflow encountered in double_scalars" expected.
        if k2_m1dt == np.inf:
            nm1[i] -= 1
            k2_m1dt = Am1[i]*gamma(nm1[i]+1)*gamma(nm1[i]-mm1[i]+len(Fr)-1+1)/gamma(nm1[i]-mm1[i]+1)/gamma(nm1[i]+len(Fr)-1+1)        
        Keqd_t = np.exp(-n[i]*vgmean/(kB_2*Te[i]*NA))*(1-np.exp(-vgmean/(kB_2*Te[i]*NA)))**(len(Fr))*gamma(n[i]+len(Fr)-1+1)/gamma(n[i]+1)/gamma(len(Fr)-1+1)
        if Keqd_t == 0.0 or Keqd_t == np.inf or Keqd_t == -np.inf:   # v2B2 FIX 7 (typo)
            n[i] -= 1
            Keqd_t = np.exp(-n[i]*vgmean/(kB_2*Te[i]*NA))*(1-np.exp(-vgmean/(kB_2*Te[i]*NA)))**(len(Fr))*gamma(n[i]+len(Fr)-1+1)/gamma(n[i]+1)/gamma(len(Fr)-1+1)
        Fd_t = Fd_t+k2_m1dt*Keqd_t
        n[i] += 1;  nm1[i] += 1
    Fd.append(Fd_t)     
    Keqd.append(Keqd_t)
    k2_m1d.append(k2_m1dt)

# Actual rate constants
n = np.array(E0)/vgmean; nm1 = np.array(E0)/vgmean; n2 = np.array(E0)/vgmean
k2_m1,k2_2,Keq,F,kstab_SS,kstab_G,kstab_D,kp_SS,kp_G,kp_D = [],[],[],[],[],[],[],[],[],[]
km1_G,km1_D = [],[]
for i in range (len(Te)):
    k2_m1t,k2_2t,Keq_t,F_t,kstab0_SS_t,kstab0_G_t,kstab0_D_t,kp0_SS_t,kp0_G_t,kp0_D_t = np.float64(0.0),np.float64(0.0),np.float64(0.0),0.0,0.0,0.0,0.0,0.0,0.0,0.0
    for j in range (MaxI):
        k2_m1t = Am1[i]*gamma(nm1[i]+1)*gamma(nm1[i]-mm1[i]+len(Fr)-1+1)/gamma(nm1[i]-mm1[i]+1)/gamma(nm1[i]+len(Fr)-1+1)       # The value of "+1" always stands for the gamma function that is of the way (n+1). This normally goes to 'inf'. "RuntimeWarning: overflow encountered in double_scalars" expected.
        if k2_m1t == np.inf:
            nm1[i] -= 1
            k2_m1t = Am1[i]*gamma(nm1[i]+1)*gamma(nm1[i]-mm1[i]+len(Fr)-1+1)/gamma(nm1[i]-mm1[i]+1)/gamma(nm1[i]+len(Fr)-1+1)
        k2_2t = A2[i]*gamma(n2[i]+1)*gamma(n2[i]-m2[i]+len(Fr)-1+1)/gamma(n2[i]-m2[i]+1)/gamma(n2[i]+len(Fr)-1+1)           # The value of "+1" always stands for the gamma function that is of the way (n+1). This normally goes to 'inf'. "RuntimeWarning: overflow encountered in double_scalars" expected.
        if k2_2t == np.inf:
            n2[i] -= 1
            k2_2t = A2[i]*gamma(n2[i]+1)*gamma(n2[i]-m2[i]+len(Fr)-1+1)/gamma(n2[i]-m2[i]+1)/gamma(n2[i]+len(Fr)-1+1)   
        Keq_t = np.exp(-n[i]*vgmean/(kB_2*Te[i]*NA))*(1-np.exp(-vgmean/(kB_2*Te[i]*NA)))**(len(Fr))*gamma(n[i]+len(Fr)-1+1)/gamma(n[i]+1)/gamma(len(Fr)-1+1)
        if Keq_t == 0.0 or Keq_t == np.inf or Keq_t == -np.inf:
            n[i] -= 1
            Keq_t = np.exp(-n[i]*vgmean/(kB_2*Te[i]*NA))*(1-np.exp(-vgmean/(kB_2*Te[i]*NA)))**(len(Fr))*gamma(n[i]+len(Fr)-1+1)/gamma(n[i]+1)/gamma(len(Fr)-1+1)
        F_t = k2_m1t*Keq_t/Fd[i]        
        kstab0_SS_t = kstab0_SS_t+kc_SS[i]*M_RK[i,:]*F_t/(kc_SS[i]*M_RK[i,:]+k2_m1t+k2_2t)
        kstab0_G_t = kstab0_G_t+kc_Gilbert[i]*M_RK[i,:]*F_t/(kc_Gilbert[i]*M_RK[i,:]+k2_m1t+k2_2t)
        kstab0_D_t = kstab0_D_t+kc_Dean[i]*M_RK[i,:]*F_t/(kc_Dean[i]*M_RK[i,:]+k2_m1t+k2_2t)        
        kp0_SS_t = kp0_SS_t+k2_2t*F_t/(kc_SS[i]*M_RK[i,:]+k2_m1t+k2_2t)
        kp0_G_t = kp0_G_t+k2_2t*F_t/(kc_Gilbert[i]*M_RK[i,:]+k2_m1t+k2_2t)
        kp0_D_t = kp0_D_t+k2_2t*F_t/(kc_Dean[i]*M_RK[i,:]+k2_m1t+k2_2t)
        n[i] += 1;  nm1[i] += 1; n2[i] += 1
    k2_m1.append(k2_m1t);k2_2.append(k2_2t) 
    Keq.append(Keq_t); F.append(F_t)
    kstab_SS.append(k1[i]*kstab0_SS_t);kstab_G.append(k1[i]*kstab0_G_t);kstab_D.append(k1[i]*kstab0_D_t)
    kp_SS.append(k1[i]*kp0_SS_t);kp_G.append(k1[i]*kp0_G_t);kp_D.append(k1[i]*kp0_D_t)
    km1_G.append(k1[i]*(1-kstab0_G_t-kp0_G_t));km1_D.append(k1[i]*(1-kstab0_D_t-kp0_D_t))
    
names1 = list(P); names1.insert(0,'T (K)'); names1.insert(len(P)+1,'Bc')    
names2 = list(P); names2.insert(0,'T (K)');

kstab_SS = np.matrix(kstab_SS); kp_SS = np.matrix(kp_SS) 
kstab_SS = np.column_stack((np.transpose(np.matrix(Te)),kstab_SS,np.array(Bc_SS))); kp_SS = np.column_stack((np.transpose(np.matrix(Te)),kp_SS))
kstab_G = np.matrix(kstab_G); kp_G = np.matrix(kp_G)
kstab_G = np.column_stack((np.transpose(np.matrix(Te)),kstab_G,np.array(Bc_Gilbert))); kp_G = np.column_stack((np.transpose(np.matrix(Te)),kp_G))
kstab_D = np.matrix(kstab_D); kp_D = np.matrix(kp_D)
km1_G = np.matrix(km1_G); km1_D = np.matrix(km1_D)
kstab_D = np.column_stack((np.transpose(np.matrix(Te)),kstab_D,np.array(Bc_Dean))); kp_D = np.column_stack((np.transpose(np.matrix(Te)),kp_D))
km1_G = np.column_stack((np.transpose(np.matrix(Te)),km1_G));km1_D = np.column_stack((np.transpose(np.matrix(Te)),km1_D))
df_kstab_SS = pd.DataFrame(data=kstab_SS, columns=names1); df_kp_SS = pd.DataFrame(data=kp_SS, columns=names2)
df_kstab_G = pd.DataFrame(data=kstab_G, columns=names1); df_kp_G = pd.DataFrame(data=kp_G, columns=names2)
df_kstab_D = pd.DataFrame(data=kstab_D, columns=names1); df_kp_D = pd.DataFrame(data=kp_D, columns=names2)
df_km1_G = pd.DataFrame(data=km1_G, columns=names2);df_km1_D = pd.DataFrame(data=km1_D, columns=names2)  

writer = pd.ExcelWriter('P-Dependent ' + file, engine='xlsxwriter')
df_kstab_SS.to_excel(writer,index=False,sheet_name='SS-QRRK_MSC-k_stab'); df_kp_SS.to_excel(writer,index=False,sheet_name='SS-QRRK_MSC-k_p') 
df_kstab_G.to_excel(writer,index=False,sheet_name='SS-QRRK_MSC_G-k_stab'); df_kp_G.to_excel(writer,index=False,sheet_name='SS-QRRK_MSC_G-k_p') 
df_kstab_D.to_excel(writer,index=False,sheet_name='SS-QRRK_MSC_D-k_stab'); df_kp_D.to_excel(writer,index=False,sheet_name='SS-QRRK_MSC_D-k_p')
df_km1_G.to_excel(writer,index=False,sheet_name='SS-QRRK_MSC_G-k_m1'); df_km1_D.to_excel(writer,index=False,sheet_name='SS-QRRK_MSC_D-k_m1')

writer.close()                                                                                                            # v2B2 FIX 2: writer.save() removed in pandas 2.0

print("--- %s seconds ---" % (time.time() - start_time))





