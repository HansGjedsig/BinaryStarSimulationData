import numpy as np
import random
import pandas as pd
import time
from BinaryStarFunctions import*
from LoadbarFunc import progress_bar_while
t_start_all = time.time()


sun_mass = 1.989*10**30 #kg
sun_radius = 696342*1000 # sol radius [m]

#NUMBER OS SIMS EACH LOOP
res = 10**3 #max 10**3

#NUMBER OF constricted systems wanted
iterate = 5*10**5

keep=True #True if only constricted systems are wanted
dataToExcel=False #true if data export to excel
dataToCSV = True #true if data export to CSV
dataToTxt = True #true if systemdata is wanted in a txt

name = 'PSR J2108+4516' #folder of same name is needed


#pre SN variable ranges

M_NS_min, M_NS_max= 1.4, 1.4
w_min, w_max= 0, 230
M1_min, M1_max= 1.55, 8
M2_min, M2_max,M2_type,flatdist= 17.5, 23, 'MSS', False
ang_min, ang_max= 50.3, 58.3
p_orb_min, p_orb_max= 190, 320


#post SN constraints (set to +- inf if no constraints needed)

v_sys_min, v_sys_max= -np.inf, np.inf
misAng_min, misAng_max= -np.inf, np.inf
e_min, e_max= 0.09-(0.09*0.03), 0.09+(0.09*0.03)
p_orbf_min, p_orbf_max= 269-(269*0.03), 269+(269*0.03)
alpha_ej_min, alpha_ej_max = -np.inf, np.inf



##################################################################################

post_SN_variable_logic = [(M_NS_max-M_NS_min)!=0,
                           (w_max-w_min)!=0,
                           (M1_max-M1_min)!=0,
                           (M2_max-M2_min)!=0,
                           (p_orb_max-p_orb_min)!=0]

post_SN_variable_all = [['Mass of newborn NS',M_NS_min, M_NS_max,'M$_{\odot}$'],
                        ['Kick velocity',w_min, w_max,'km/s'],
                        ['Mass of star going SN',M1_min, M1_max,'M$_{\odot}$'],
                        ['Mass of companion star',M2_min, M2_max,'M$_{\odot}$'],
                        ['Pre-SN orbital period',p_orb_min, p_orb_max,'d']]

constraints_logic = [np.isfinite(v_sys_min) & np.isfinite(v_sys_max),
                     np.isfinite(misAng_min) & np.isfinite(misAng_max),
                     np.isfinite(e_min) & np.isfinite(e_max),
                     np.isfinite(p_orbf_min) & np.isfinite(p_orbf_max),
                     np.isfinite(alpha_ej_min) & np.isfinite(alpha_ej_max)]
constraints_all = [['Systemic recoil velocity',v_sys_min,v_sys_max, 'km/s'],
                   ['Misalignment angle',misAng_min,misAng_max, 'Degree'],
                   ['Eccentricity',e_min,e_max,' '],
                   ['Orbital period',p_orbf_min,p_orbf_max,'d'],
                   ['$\\alpha_{ej}$',alpha_ej_min,alpha_ej_max,' ']]

w_sys = []
theta_sys = []
phi_sys = []
M1_sys = []
M2_sys = []
R2_sys = []
M_NS_sys = []
p_orb_sys = []
p_orb_f_sys = []
e_sys = []
a_i_sys = []
a_f_sys = []
v_sys_sys = []
v_rel_sys = []
thetacrit_sys = []
misAng_sys = []
Viable_sys = []
alpha_ej_sys =[]
N_sims = 0
N_viable = 0
t2=0
sinang_min, sinang_max = 1-np.cos(ang_min*np.pi/180), 1-np.cos(ang_max*np.pi/180)
while N_viable < iterate:
    t_start = time.time()
    
    #Creates random kick angles using the inverse transformation method
    theta = np.arcsin(np.sqrt(np.random.rand(res)))*2
    phi = np.random.rand(len(theta))*2*np.pi
    
    
    ang = np.arccos(1-(sinang_min+np.random.rand(res)*(sinang_max-sinang_min)))*180/np.pi 
    M2 = M2_max-(ang-ang_min)/(ang_max-ang_min)*(M2_max-M2_min)
    if flatdist:
        M2 = M2_max-np.random.rand(len(theta))*(M2_max-M2_min)
    
    M2 = M2*sun_mass
    R2 = RadiusMainSequence(M2)
    if M2_type == 'NS':
        R2 = np.ones(len(M2))*15000
    
    M_NS = (M_NS_min+np.random.rand(res)*(M_NS_max-M_NS_min))*sun_mass
    
    M1 = (M1_min+np.random.rand(res)*(M1_max-M1_min))*sun_mass
    
    p_orb = (p_orb_min+np.random.rand(res)*(p_orb_max-p_orb_min))*24*60*60
    
    w = (w_min+np.random.rand(res)*(w_max-w_min))*1000
    
    M = M1 + M2
    DeltaM = M1 - M_NS
    M_f = M2 + M_NS
    
    a_i = SemimajorAxKepler3(M, p_orb)
    rocheLobe_i = RocheLobeRadius(M2, M1, a_i)
    
    
    v_rel = RelativeVelocity(M, a_i)
    thetacrit = ThetaCritical(M, DeltaM, w, v_rel)
    
    
    misAng = MisalignmentAngle(v_rel, w, theta, phi)
    v_sys = SystemicVelocity(M_NS, M1, M2, w, a_i, theta, phi)
    
    a_f = PostSNSemimajorAx(M, DeltaM, v_rel, a_i, w, theta)
    e = PostSNeccentricity(M_NS, M2, a_i, a_f, v_rel, w, theta, phi)
    rocheLobe_f = RocheLobeRadius(M2, M_NS, a_f*(1-e))
    
    p_orb_f = OrbitalPeriodKepler3(a_f, M_f)
    
    alpha_ej = (w/1000)/211*(0.1/(DeltaM/sun_mass))*((M_NS/sun_mass)/1.5)
    
    #LOGIC TESTS
    close = R2<rocheLobe_i
    
    bounded = theta>thetacrit
    
    unmerged = rocheLobe_f>R2
    
    p_orb_check = (p_orb_f>(p_orbf_min*24*60*60)) & (p_orb_f<(p_orbf_max*24*60*60))
    
    e_check = (e>e_min) & (e<e_max)
    
    misAng_check = (misAng*180/np.pi>misAng_min) & (misAng*180/np.pi<misAng_max)
    
    v_sys_check = (v_sys/1000>v_sys_min) & (v_sys/1000<v_sys_max)
    
    checks = p_orb_check & misAng_check & v_sys_check & e_check
    
    viable = bounded & unmerged & close
    
    viableChecks = viable & checks
    
    if keep:
        viable = viableChecks
    
    Viable_sys.append(checks[viable])
    R2_sys.append(R2[viable]/sun_radius)
    M_NS_sys.append(M_NS[viable]/sun_mass)
    M2_sys.append(M2[viable]/sun_mass)
    w_sys.append(w[viable]/1000)
    theta_sys.append(theta[viable]*180/np.pi)
    phi_sys.append(phi[viable]*180/np.pi)
    M1_sys.append(M1[viable]/sun_mass)
    p_orb_sys.append(p_orb[viable]/(24*60*60))
    p_orb_f_sys.append(p_orb_f[viable]/(24*60*60))
    e_sys.append(e[viable])
    a_i_sys.append(a_i[viable]/sun_radius)
    a_f_sys.append(a_f[viable]/sun_radius)
    v_sys_sys.append(v_sys[viable]/1000)
    v_rel_sys.append(v_rel[viable]/1000)
    thetacrit_sys.append(thetacrit[viable]*180/np.pi)
    misAng_sys.append(misAng[viable]*180/np.pi)
    alpha_ej_sys.append(alpha_ej[viable])
    N_sims = N_sims + res
    
    N_viable = N_viable + len(w[viableChecks])
    t2=progress_bar_while(N_viable, iterate,t_start,t2)

R2_sys = np.concatenate(R2_sys[:])
M_NS_sys = np.concatenate(M_NS_sys[:])
M2_sys = np.concatenate(M2_sys[:])
w_sys = np.concatenate(w_sys[:])
theta_sys = np.concatenate(theta_sys[:])
phi_sys = np.concatenate(phi_sys[:])
M1_sys = np.concatenate(M1_sys[:])
p_orb_sys = np.concatenate(p_orb_sys[:])
p_orb_f_sys = np.concatenate(p_orb_f_sys[:])
e_sys = np.concatenate(e_sys[:])
a_i_sys = np.concatenate(a_i_sys[:])
a_f_sys = np.concatenate(a_f_sys[:])
v_sys_sys = np.concatenate(v_sys_sys[:])
v_rel_sys = np.concatenate(v_rel_sys[:])
thetacrit_sys = np.concatenate(thetacrit_sys[:])
misAng_sys = np.concatenate(misAng_sys[:])
Viable_sys = np.concatenate(Viable_sys[:])
alpha_ej_sys = np.concatenate(alpha_ej_sys[:])

data_sys = {'Viable':Viable_sys,
            'kick [km/s]':w_sys,
            'theta [deg]': theta_sys,
            'phi [deg]': phi_sys,
            'M1 [s_m]': M1_sys,
            'M2 [s_m]': M2_sys,
            'M_NS [s_m]':M_NS_sys,
            'R2 [s_r]': R2_sys,
            'p_orb_i [d]': p_orb_sys,
            'p_orb_f [d]': p_orb_f_sys,
            'e_f':e_sys,
            'Initial seperation [s_r]': a_i_sys,
            'Final semi major axis [s_r]': a_f_sys,
            'v_sys [km/s]': v_sys_sys,
            'Initial relative velocity [km/s]': v_rel_sys,
            'Critical angle [deg]': thetacrit_sys,
            'misalignment [deg]':misAng_sys,
            'alpha_ej':alpha_ej_sys}

df = pd.DataFrame(data_sys)


if dataToExcel:
    print('\n'+'Exporting viable systems data to excel spreadsheet')
    df.to_excel(name+'/'+name+'.xlsx', index=False)

if dataToCSV:
    print('\n'+'Exporting systems data to CSV')
    df.to_csv(name+'/'+name+'.csv', index=False)

t_slut = time.time()-t_start_all
th = int(t_slut/(60*60))
tm = int((t_slut/(60*60)-th)*60)
ts =  int(((t_slut/(60*60)-th)*60-tm)*60)

if dataToTxt:
    with open(name+'/'+name+'.txt', 'w') as f:
        f.write(name)
        f.write('\n'+f"Total simulation time= {th:02d}h {tm:02d}m {ts:02d}s")
        f.write('\n'+'Number of simulations= '+str("{:e}".format(N_sims)))
        f.write('\n'+'Number of saved systems= '+str("{:e}".format(N_viable)))
        f.write('\n'+'\n'+'Simulation ranges')
        f.write('\n'+'M_NS_min, M_NS_max= '+str(M_NS_min)+', '+str(M_NS_max))
        f.write('\n'+'w_min, w_max= '+str(w_min)+', '+str(w_max))
        f.write('\n'+'M1_min, M1_max= '+str(M1_min)+', '+str(M1_max))
        f.write('\n'+'M2_min, M2_max= '+str(M2_min)+', '+str(M2_max))
        f.write('\n'+'ang_min, ang_max= '+str(ang_min)+', '+str(ang_max))
        f.write('\n'+'p_orb_min, p_orb_max= '+str(p_orb_min)+', '+str(p_orb_max))
        f.write('\n'+'\n'+'constraints')
        if constraints_logic[0]:
            f.write('\n'+'v_sys_min, v_sys_max= '+str(v_sys_min)+', '+str(v_sys_max))
        else:
            f.write('\n'+'v_sys_min, v_sys_max= -np.inf, np.inf')
        if constraints_logic[1]:
            f.write('\n'+'misAng_min, misAng_max= '+str(misAng_min)+', '+str(misAng_max))
        else:
            f.write('\n'+'misAng_min, misAng_max= -np.inf, np.inf')
        if constraints_logic[2]:
            f.write('\n'+'e_min, e_max= '+str(e_min)+', '+str(e_max))
        else:
            f.write('\n'+'e_min, e_max= -np.inf, np.inf')
        if constraints_logic[3]:
            f.write('\n'+'p_orbf_min, p_orbf_max= '+str(p_orbf_min)+', '+str(p_orbf_max))
        else:
            f.write('\n'+'p_orbf_min, p_orbf_max= -np.inf, np.inf')
        if constraints_logic[4]:
            f.write('\n'+'alpha_ej_min, alpha_ej_max= '+str(alpha_ej_min)+', '+str(alpha_ej_max))
        else:
            f.write('\n'+'alpha_ej_min, alpha_ej_max= -np.inf, np.inf')
        
        
        


print('\n'+'Done! YaY')