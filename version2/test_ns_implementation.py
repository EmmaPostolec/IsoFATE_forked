
from constants import *
from orbit_params import *
from isofate_coupler import isocalc
from atmodeller_coupler import *
from atmodeller import *

# M1 star
R_star = 0.5*Rs # [m]
M_star = 0.5*Ms
T_star = 3600 # [K]
L = Luminosity(R_star, T_star) # [W]
Mp = 5*Me
P = 12/s2day
f_atm = 0.01
a = SemiMajor(M_star, P) # [m]
Fp = Insolation(L, a)  # [W/m2]
T = EqTemp(Fp, A = 0) # planetary eq temp [K]
F0 = Fp*1e-3 # use for M star 
flux_model = 'power law'
stellar_type = 'M1'
t_sat = 5e8 # XUV saturation time [yr]
d = a # orbital distance [m]
time = 5e10 # [yr]
mechanism = 'XUV' # if using fixed phi, be sure to change Rp = r_core below and rad_evol = False
Johnson = False
RR = True
rad_evol = True
Rp_override = False
n_TO_final = 'null'
f_atm_final = 'null'
n_steps = int(1e5)
species = 'H/He/D/O/C/N/S'
f_pred = False
thermal = True
H2O_reservoir = False
wmf = 0
M_atm = Mp*f_atm # initial atmospheric mass [kg]

mu_avg = (N_H*mu_H + N_He*mu_He + N_D*mu_D + N_O*mu_O + N_C*mu_C +N_N*mu_N + N_S*mu_S)/(N_H+N_He+N_D+N_O+N_C+N_N+N_S) # average molecular weight [kg/mol]


if species == 'H/He/D/O/C/N/S':
    N_H = (12.6/13.6)*(1 - OtoH_protosolar - CtoH_protosolar - DtoH_solar)*M_atm/mu_avg # initial H number [atoms]
    N_O = (12.6/13.6)*OtoH_protosolar*M_atm/mu_avg # initial O number [atoms]
    N_C = (12.6/13.6)*CtoH_protosolar*M_atm/mu_avg # initial O number [atoms]
    N_D = (12.6/13.6)*DtoH_solar*M_atm/mu_avg # initial D number [atoms]
    N_N = (12.6/13.6)*NtoH_protosolar*M_atm/mu_avg
    N_S = (12.6/13.6)*StoH_protosolar*M_atm/mu_avg




test = isocalc(f_atm, Mp, Mstar, F0, Fp, T, d, time = 5e9, mechanism = 'XUV', rad_evol = True,
N_H = 0, N_He = 0, N_D = 0, N_O = 0, N_C = 0, N_N=0, N_S=0, melt_fraction_override = False,
mu = mu_solar, eps = 0.15, activity = 'medium', flux_model = 'power law', stellar_type = 'M1',
Rp_override = False, t_sat = 5e8,
n_steps = int(1e5), t0 = 1e6, rho_rcb = 1.0, RR = True, thermal = True, 
beta = -1.23, n_atmodeller = int(1e2), save_molecules = False, mantle_iron_dict = False)

print(test)