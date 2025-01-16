#! /usr/bin/env python3.8
import sys,matplotlib
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm
sys.path.append('../../utils/')
matplotlib.use('Agg')
from numpy import *
from pylab import *
import getopt
import os, re
import pdb
import math as math
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d, interp2d, RectBivariateSpline
from time import time as pytime

#./cloud_layers.py ../run/carma/hjt_west_1800 xt 1800_west layer_plots/ > log xt &

# Read in inputs from command line
opts, args = getopt.getopt(sys.argv[1:], "h")
file = args[0]
name = args[1]
path = args[2]

# Check if folder where the plots are generated exists; if 
# not, make one. 
if not os.path.isdir(path):
	os.makedirs(path)

start_time = pytime()  # Record start time

# Read in given output file and plot results. 
print('Read in output file')
infile=open(file,'r')
line = infile.readline().split()

longitude = np.loadtxt('../run/carma/longitudes_WASP-52b_clear_1xsolar.txt')
ilong = len(longitude)

nz,ngroup,nelem,nbin,ngas,nstep,iskip = map(int,line)
r,ms,dr,rl,ru = zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup)),zeros((nbin,ngroup))
p,t,z = zeros(nz),zeros(nz),zeros(nz)

ntime=int(nstep/iskip)
print(nstep,iskip,ntime)

time, distance, rotation = zeros(ntime), zeros(ntime), zeros(ntime)

tio2, het_fe, het_mg2sio4, pure_fe, het_cr, cr, mns, na2s,zns,kcl,al2o3 = zeros((nz,nbin,ilong,ntime)),zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime))
mass_tio2, mass_het_fe, mass_het_mg2sio4, mass_pure_fe, mass_het_cr, mass_cr, mass_mns, mass_na2s, mass_zns, mass_kcl, mass_al2o3 = zeros((nz,nbin,ilong,ntime)),zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime)), zeros((nz,nbin,ilong,ntime))

t_tio2, t_het_fe, t_het_mg2sio4, t_pure_fe, t_het_cr, t_cr, t_mns, t_na2s,t_zns,t_kcl,t_al2o3 = zeros((nz,nbin,ilong)),zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong))
t_mass_tio2, t_mass_het_fe, t_mass_het_mg2sio4, t_mass_pure_fe, t_mass_het_cr, t_mass_cr, t_mass_mns, t_mass_na2s, t_mass_zns, t_mass_kcl, t_mass_al2o3 = zeros((nz,nbin,ilong)),zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong)), zeros((nz,nbin,ilong))

p_ti, p_fe, p_mg, p_cr, p_mn, p_na,p_zn,p_k,p_al = zeros((nz,nbin)),zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))
pm_ti, pm_fe, pm_mg, pm_cr, pm_mn, pm_na,pm_zn,pm_k,pm_al = zeros((nz,nbin)),zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin)), zeros((nz,nbin))
p_pure_cr,p_pure_fe = zeros((nz,nbin)),zeros((nz,nbin))
pm_pure_cr,pm_pure_fe = zeros((nz,nbin)),zeros((nz,nbin))

mmr_ti, mmr_fe, mmr_mg, mmr_cr, mmr_zn, mmr_al, mmr_kc = zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime))
svp_ti, svp_fe, svp_mg, svp_cr, svp_zn, svp_al, svp_kc = zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime)), zeros((nz,ilong,ntime))

t_mmr_ti, t_mmr_fe, t_mmr_mg, t_mmr_cr, t_mmr_zn, t_mmr_al, t_mmr_kc = zeros((nz,ilong)), zeros((nz,ilong)), zeros((nz,ilong)), zeros((nz,ilong)), zeros((nz,ilong)), zeros((nz,ilong)), zeros((nz,ilong))
t_svp_ti, t_svp_fe, t_svp_mg, t_svp_cr, t_svp_zn, t_svp_al, t_svp_kc = zeros((nz,ilong)), zeros((nz,ilong)), zeros((nz,ilong)), zeros((nz,ilong)), zeros((nz,ilong)), zeros((nz,ilong)), zeros((nz,ilong))

M_KCL = 1.98
M_ZNS = 4.1
M_NA2S = 1.856
M_MNS = 3.3
M_CR = 7.19
M_FE = 7.874
M_MG2SIO4 = 3.21
M_TIO2 = 4.23
M_AL2O3 = 3.95

for i in range(ngroup):
	for j in range(nbin):
		line = infile.readline().split()
		r[j,i] = float(line[2])
		ms[j,i] = float(line[3])
		dr[j,i] = float(line[4])
		rl[j,i] = float(line[5])
		ru[j,i] = float(line[6])
for i in range(nz):
	line = infile.readline().split()
	z[i] = float(line[1])
	p[i] = float(line[3])
	t[i] = float(line[4])
			
#np.savetxt('new_ti_rad.txt', r[:,0] )
#np.savetxt('new_al_rad.txt', r[:,1] )
#np.savetxt('new_fe_rad.txt', r[:,2] )
#np.savetxt('new_mg_rad.txt', r[:,3] )
#np.savetxt('new_cr_rad.txt', r[:,4] )
#np.savetxt('new_mn_rad.txt', r[:,5] )
#np.savetxt('new_na_rad.txt', r[:,6] )
#np.savetxt('new_pfe_rad.txt', r[:,7] )
#np.savetxt('new_pcr_rad.txt', r[:,8] )
#np.savetxt('new_k_rad.txt', r[:,9] )
#np.savetxt('new_zn_rad.txt', r[:,10] )

plt.rcParams['figure.figsize'] = (12,6)
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['font.size'] = 20

levels_num = np.array([1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8])/1e7
levels_mass = np.array([1e-16,1e-14,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2])/1e7

count = np.zeros((ilong))

for i in range(ntime-1):
	if i == 0:
		time[i] = float(infile.readline())
		long_index = 0
	else:
		td  = infile.readline().split()
		time[i] = float(td[0])
		distance[i] = float(td[1])
		rotation[i] = float(td[2])
		current_longitude = distance[i] - (len(longitude)*rotation[i])
		if int(round(current_longitude,0)) == ilong:
			long_index = 0
		else:
			long_index = int(round(current_longitude,0))
	for j in range(nbin):
		for k in range(nz):
			print(i,j,k)
			line = infile.readline().split()
			
			tio2[k,j,long_index,i] = float(line[2])/(math.log(ru[j,0]) - math.log(rl[j,0]))
			mass_tio2[k,j,long_index,i] = tio2[k,j,long_index,i]*(M_TIO2*((4./3.)*np.pi*(r[j,0]*1e-4)**3.))
			
			al2o3[k,j,long_index,i] = float(line[3])/(math.log(ru[j,1]) - math.log(rl[j,1]))
			mass_al2o3[k,j,long_index,i] = al2o3[k,j,long_index,i]*(M_AL2O3*((4./3.)*np.pi*(r[j,1]*1e-4)**3.))
			
			het_fe[k,j,long_index,i] = float(line[4])/(math.log(ru[j,2]) - math.log(rl[j,2]))
			mass_het_fe[k,j,long_index,i] = het_fe[k,j,long_index,i]*(M_FE*((4./3.)*np.pi*(r[j,2]*1e-4)**3.))
			
			het_mg2sio4[k,j,long_index,i] = float(line[5])/(math.log(ru[j,3]) - math.log(rl[j,3]))
			mass_het_mg2sio4[k,j,long_index,i] = het_mg2sio4[k,j,long_index,i]*(M_MG2SIO4*((4./3.)*np.pi*(r[j,3]*1e-4)**3.))
			
			het_cr[k,j,long_index,i] = float(line[6])/(math.log(ru[j,4]) - math.log(rl[j,4]))
			mass_het_cr[k,j,long_index,i] = het_cr[k,j,long_index,i]*(M_CR*((4./3.)*np.pi*(r[j,4]*1e-4)**3.))
			
			mns[k,j,long_index,i] = float(line[7])/(math.log(ru[j,5]) - math.log(rl[j,5]))
			mass_mns[k,j,long_index,i] = mns[k,j,long_index,i]*(M_MNS*((4./3.)*np.pi*(r[j,5]*1e-4)**3.))
			
			na2s[k,j,long_index,i] = float(line[8])/(math.log(ru[j,6]) - math.log(rl[j,6]))
			mass_na2s[k,j,long_index,i] = na2s[k,j,long_index,i]*(M_NA2S*((4./3.)*np.pi*(r[j,6]*1e-4)**3.))
			
			pure_fe[k,j,long_index,i] = float(line[9])/(math.log(ru[j,7]) - math.log(rl[j,7]))
			mass_pure_fe[k,j,long_index,i] = pure_fe[k,j,long_index,i]*(M_FE*((4./3.)*np.pi*(r[j,7]*1e-4)**3.))
			
			cr[k,j,long_index,i] = float(line[10])/(math.log(ru[j,8]) - math.log(rl[j,8]))
			mass_cr[k,j,long_index,i] = cr[k,j,long_index,i]*(M_CR*((4./3.)*np.pi*(r[j,8]*1e-4)**3.))
			
			kcl[k,j,long_index,i] = float(line[11])/(math.log(ru[j,9]) - math.log(rl[j,9]))
			mass_kcl[k,j,long_index,i] = kcl[k,j,long_index,i]*(M_KCL*((4./3.)*np.pi*(r[j,9]*1e-4)**3.))
			
			zns[k,j,long_index,i] = float(line[12])/(math.log(ru[j,10]) - math.log(rl[j,10]))
			mass_zns[k,j,long_index,i] = zns[k,j,long_index,i]*(M_ZNS*((4./3.)*np.pi*(r[j,10]*1e-4)**3.))
			
			mmr_ti[k,long_index,i] = (float(line[13])/1.e6)*p[k]/1.e6 #partial pressure in bar
			svp_ti[k,long_index,i] = (float(line[14])/1.e6)*p[k]/1.e6
			
			mmr_fe[k,long_index,i] = (float(line[15])/1.e6)*p[k]/1.e6 #partial pressure in bar
			svp_fe[k,long_index,i] = (float(line[16])/1.e6)*p[k]/1.e6
			
			mmr_mg[k,long_index,i] = (float(line[17])/1.e6)*p[k]/1.e6 #partial pressure in bar
			svp_mg[k,long_index,i] = (float(line[18])/1.e6)*p[k]/1.e6
			
			mmr_cr[k,long_index,i] = (float(line[19])/1.e6)*p[k]/1.e6 #partial pressure in bar
			svp_cr[k,long_index,i] = (float(line[20])/1.e6)*p[k]/1.e6

			#21
			#22

			#23
			#24

			#25
			#26
			
			mmr_kc[k,long_index,i] = (float(line[27])/1.e6)*p[k]/1.e6 #partial pressure in bar
			svp_kc[k,long_index,i] = (float(line[28])/1.e6)*p[k]/1.e6
			
			mmr_al[k,long_index,i] = (float(line[29])/1.e6)*p[k]/1.e6 #partial pressure in bar
			svp_al[k,long_index,i] = (float(line[30])/1.e6)*p[k]/1.e6
			
			
			
	if i > 1 :
		count[long_index] += 1
			
		t_tio2[:,:,long_index] += tio2[:,:,long_index,i]
		t_mass_tio2[:,:,long_index] += mass_tio2[:,:,long_index,i]
			
		t_het_fe[:,:,long_index] += het_fe[:,:,long_index,i]
		t_mass_het_fe[:,:,long_index] += mass_het_fe[:,:,long_index,i]
			
		t_het_mg2sio4[:,:,long_index] += het_mg2sio4[:,:,long_index,i]
		t_mass_het_mg2sio4[:,:,long_index] += mass_het_mg2sio4[:,:,long_index,i]
			
		t_pure_fe[:,:,long_index] += pure_fe[:,:,long_index,i]
		t_mass_pure_fe[:,:,long_index] += mass_pure_fe[:,:,long_index,i]
			
		t_het_cr[:,:,long_index] += het_cr[:,:,long_index,i]
		t_mass_het_cr[:,:,long_index] += mass_het_cr[:,:,long_index,i]
			
		t_cr[:,:,long_index] += cr[:,:,long_index,i]
		t_mass_cr[:,:,long_index] += mass_cr[:,:,long_index,i]
			
		t_mns[:,:,long_index] += mns[:,:,long_index,i]
		t_mass_mns[:,:,long_index] += mass_mns[:,:,long_index,i]
			
		t_na2s[:,:,long_index] += na2s[:,:,long_index,i]
		t_mass_na2s[:,:,long_index] += mass_na2s[:,:,long_index,i]
			
		t_zns[:,:,long_index] += zns[:,:,long_index,i]
		t_mass_zns[:,:,long_index] += mass_zns[:,:,long_index,i]
			
		t_kcl[:,:,long_index] += kcl[:,:,long_index,i]
		t_mass_kcl[:,:,long_index] += mass_kcl[:,:,long_index,i]
			
		t_al2o3[:,:,long_index] += al2o3[:,:,long_index,i]
		t_mass_al2o3[:,:,long_index] += mass_al2o3[:,:,long_index,i]
		
		t_mmr_ti[:,long_index] += mmr_ti[:,long_index,i]
		t_svp_ti[:,long_index] += svp_ti[:,long_index,i]
		 
		t_mmr_fe[:,long_index] += mmr_fe[:,long_index,i]
		t_svp_fe[:,long_index] += svp_fe[:,long_index,i]
		 
		t_mmr_mg[:,long_index] += mmr_mg[:,long_index,i]
		t_svp_mg[:,long_index] += svp_mg[:,long_index,i]
		 
		t_mmr_cr[:,long_index] += mmr_cr[:,long_index,i]
		t_svp_cr[:,long_index] += svp_cr[:,long_index,i]
		 
		t_mmr_kc[:,long_index] += mmr_kc[:,long_index,i]
		t_svp_kc[:,long_index] += svp_kc[:,long_index,i]
	 
		t_mmr_al[:,long_index] += mmr_al[:,long_index,i]
		t_svp_al[:,long_index] += svp_al[:,long_index,i]


print(count[:], 'number per longitude', longitude)					

for il in range(ilong):
	t_tio2[:,:,il] /= count[il]
	t_mass_tio2[:,:,il] /= count[il] 
	t_het_fe[:,:,il] /= count[il]
	t_mass_het_fe[:,:,il] /= count[il]
	t_het_mg2sio4[:,:,il] /= count[il]
	t_mass_het_mg2sio4[:,:,il] /= count[il]
	t_pure_fe[:,:,il] /= count[il]
	t_mass_pure_fe[:,:,il] /= count[il]
	t_het_cr[:,:,il] /= count[il]
	t_mass_het_cr[:,:,il] /= count[il]
	t_cr[:,:,il] /= count[il]
	t_mass_cr[:,:,il] /= count[il]
	t_mns[:,:,il] /= count[il]
	t_mass_mns[:,:,il] /= count[il]
	t_na2s[:,:,il] /= count[il]
	t_mass_na2s[:,:,il] /= count[il]
	t_zns[:,:,il] /= count[il]
	t_mass_zns[:,:,il] /= count[il]
	t_kcl[:,:,il] /= count[il]
	t_mass_kcl[:,:,il] /= count[il]
	t_al2o3[:,:,il] /= count[il]
	t_mass_al2o3[:,:,il] /= count[il]
	
	t_mmr_ti[:,il] /= count[il]
	t_svp_ti[:,il] /= count[il]
	t_mmr_fe[:,il] /= count[il]
	t_svp_fe[:,il] /= count[il]
	t_mmr_mg[:,il] /= count[il]
	t_svp_mg[:,il] /= count[il]
	t_mmr_cr[:,il] /= count[il]
	t_svp_cr[:,il] /= count[il]
	t_mmr_kc[:,il] /= count[il]
	t_svp_kc[:,il] /= count[il]
	t_mmr_al[:,il] /= count[il]
	t_svp_al[:,il] /= count[il]

end_time = pytime()  # Record end time
elapsed_time = end_time - start_time 
print('Time taken: ',elapsed_time)

'''
num_density_ti = zeros((nz,ilong))
num_density_mg = zeros((nz,ilong))
num_density_fe = zeros((nz,ilong))
num_density_cr = zeros((nz,ilong))
num_density_mn = zeros((nz,ilong))
num_density_na = zeros((nz,ilong))
num_density_zn = zeros((nz,ilong))
num_density_kc = zeros((nz,ilong))
num_density_al = zeros((nz,ilong))

mass_density_ti = zeros((nz,ilong))
mass_density_mg = zeros((nz,ilong))
mass_density_fe = zeros((nz,ilong))
mass_density_cr = zeros((nz,ilong))
mass_density_mn = zeros((nz,ilong))
mass_density_na = zeros((nz,ilong))
mass_density_zn = zeros((nz,ilong))
mass_density_kc = zeros((nz,ilong))
mass_density_al = zeros((nz,ilong))

ps_ti_mass = zeros((nz,ilong))
ps_mg_mass = zeros((nz,ilong))
ps_fe_mass = zeros((nz,ilong))
ps_cr_mass = zeros((nz,ilong))
ps_mn_mass = zeros((nz,ilong))
ps_na_mass = zeros((nz,ilong))
ps_zn_mass = zeros((nz,ilong))
ps_kc_mass = zeros((nz,ilong))
ps_al_mass = zeros((nz,ilong))

for k in range(nz):
	for il in range(ilong):
		num_density_ti[k,il] = np.sum(t_tio2[k,:,il])
		num_density_mg[k,il] = np.sum(t_het_mg2sio4[k,:,il])
		num_density_fe[k,il] = np.sum(t_het_fe[k,:,il]+t_pure_fe[k,:,il])
		num_density_cr[k,il] = np.sum(t_cr[k,:,il]+t_het_cr[k,:,il])
		num_density_mn[k,il] = np.sum(t_mns[k,:,il])
		num_density_na[k,il] = np.sum(t_na2s[k,:,il])
		num_density_zn[k,il] = np.sum(t_zns[k,:,il])
		num_density_kc[k,il] = np.sum(t_kcl[k,:,il])
		num_density_al[k,il] = np.sum(t_al2o3[k,:,il])
		
		mass_density_ti[k,il] = np.sum(t_mass_tio2[k,:,il])
		mass_density_mg[k,il] = np.sum(t_mass_het_mg2sio4[k,:,il])
		mass_density_fe[k,il] = np.sum(t_mass_het_fe[k,:,il]+t_mass_pure_fe[k,:,il])
		mass_density_cr[k,il] = np.sum(t_mass_cr[k,:,il]+t_mass_het_cr[k,:,il])
		mass_density_mn[k,il] = np.sum(t_mass_mns[k,:,il])
		mass_density_na[k,il] = np.sum(t_mass_na2s[k,:,il])
		mass_density_zn[k,il] = np.sum(t_mass_zns[k,:,il])
		mass_density_kc[k,il] = np.sum(t_mass_kcl[k,:,il])
		mass_density_al[k,il] = np.sum(t_mass_al2o3[k,:,il])
			
		ps_ti_mass[k,il] = np.sum(t_mass_tio2[k,:,il])/(num_density_ti[k,il])
		ps_ti_mass[k,il] = (ps_ti_mass[k,il]/M_TIO2/(4./3.)/np.pi)**(1./3.)*1.e4 #microns	
		ps_mg_mass[k,il] = np.sum(t_mass_het_mg2sio4[k,:,il])/(num_density_mg[k,il])
		ps_mg_mass[k,il] = (ps_mg_mass[k,il]/M_MG2SIO4/(4./3.)/np.pi)**(1./3.)*1.e4
		ps_fe_mass[k,il] = np.sum(t_mass_het_fe[k,:,il]+t_mass_pure_fe[k,:,il])/(num_density_fe[k,il])
		ps_fe_mass[k,il] = (ps_fe_mass[k,il]/M_FE/(4./3.)/np.pi)**(1./3.)*1.e4
		ps_cr_mass[k,il] = np.sum(t_mass_cr[k,:,il]+t_mass_het_cr[k,:,il])/(num_density_cr[k,il])
		ps_cr_mass[k,il] = (ps_cr_mass[k,il]/M_CR/(4./3.)/np.pi)**(1./3.)*1.e4
		ps_mn_mass[k,il] = np.sum(t_mass_mns[k,:,il])/(num_density_mn[k,il])
		ps_mn_mass[k,il] = (ps_mn_mass[k,il]/M_MNS/(4./3.)/np.pi)**(1./3.)*1.e4
		ps_na_mass[k,il] = np.sum(t_mass_na2s[k,:,il])/(num_density_na[k,il])
		ps_na_mass[k,il] = (ps_na_mass[k,il]/M_NA2S/(4./3.)/np.pi)**(1./3.)*1.e4
		ps_zn_mass[k,il] = np.sum(t_mass_zns[k,:,il])/(num_density_zn[k,il])
		ps_zn_mass[k,il] = (ps_zn_mass[k,il]/M_ZNS/(4./3.)/np.pi)**(1./3.)*1.e4
		ps_kc_mass[k,il] = np.sum(t_mass_kcl[k,:,il])/(num_density_kc[k,il])
		ps_kc_mass[k,il] = (ps_kc_mass[k,il]/M_KCL/(4./3.)/np.pi)**(1./3.)*1.e4
		ps_al_mass[k,il] = np.sum(t_mass_al2o3[k,:,il])/(num_density_al[k,il])
		ps_al_mass[k,il] = (ps_al_mass[k,il]/M_AL2O3/(4./3.)/np.pi)**(1./3.)*1.e4


np.savetxt(path+'/num_density'+str(name)+'_ti.txt', num_density_ti.T)
np.savetxt(path+'/num_density'+str(name)+'_mg.txt', num_density_mg.T)
np.savetxt(path+'/num_density'+str(name)+'_fe.txt', num_density_fe.T)
np.savetxt(path+'/num_density'+str(name)+'_cr.txt', num_density_cr.T)
np.savetxt(path+'/num_density'+str(name)+'_mn.txt', num_density_mn.T)
np.savetxt(path+'/num_density'+str(name)+'_na.txt', num_density_na.T)
np.savetxt(path+'/num_density'+str(name)+'_zn.txt', num_density_zn.T)
np.savetxt(path+'/num_density'+str(name)+'_kc.txt', num_density_kc.T)
np.savetxt(path+'/num_density'+str(name)+'_al.txt', num_density_al.T)

np.savetxt(path+'/mass_density'+str(name)+'_ti.txt', mass_density_ti.T)
np.savetxt(path+'/mass_density'+str(name)+'_mg.txt', mass_density_mg.T)
np.savetxt(path+'/mass_density'+str(name)+'_fe.txt', mass_density_fe.T)
np.savetxt(path+'/mass_density'+str(name)+'_cr.txt', mass_density_cr.T)
np.savetxt(path+'/mass_density'+str(name)+'_mn.txt', mass_density_mn.T)
np.savetxt(path+'/mass_density'+str(name)+'_na.txt', mass_density_na.T)
np.savetxt(path+'/mass_density'+str(name)+'_zn.txt', mass_density_zn.T)
np.savetxt(path+'/mass_density'+str(name)+'_kc.txt', mass_density_kc.T)
np.savetxt(path+'/mass_density'+str(name)+'_al.txt', mass_density_al.T)


np.savetxt(path+'/mmr'+str(name)+'_ti.txt', t_mmr_ti.T)
np.savetxt(path+'/mmr'+str(name)+'_mg.txt', t_mmr_mg.T)
np.savetxt(path+'/mmr'+str(name)+'_fe.txt', t_mmr_fe.T)
np.savetxt(path+'/mmr'+str(name)+'_cr.txt', t_mmr_cr.T)
np.savetxt(path+'/mmr'+str(name)+'_kc.txt', t_mmr_kc.T)
np.savetxt(path+'/mmr'+str(name)+'_al.txt', t_mmr_al.T)

np.savetxt(path+'/svp'+str(name)+'_ti.txt', t_svp_ti.T)
np.savetxt(path+'/svp'+str(name)+'_mg.txt', t_svp_mg.T)
np.savetxt(path+'/svp'+str(name)+'_fe.txt', t_svp_fe.T)
np.savetxt(path+'/svp'+str(name)+'_cr.txt', t_svp_cr.T)
np.savetxt(path+'/svp'+str(name)+'_kc.txt', t_svp_kc.T)
np.savetxt(path+'/svp'+str(name)+'_al.txt', t_svp_al.T)



#t_mg2sio4 = t_het_mg2sio4
#t_fe = t_het_fe + t_pure_fe
#t_cr_out = t_het_cr + t_cr
#
#np.savetxt(path+'/num_t_'+str(name)+'_ti.txt', t_tio2.T)
#np.savetxt(path+'/num_t_'+str(name)+'_mg.txt', t_mg2sio4.T)
#np.savetxt(path+'/num_t_'+str(name)+'_fe.txt', t_fe.T)
#np.savetxt(path+'/num_t_'+str(name)+'_cr.txt', t_cr_out.T)
#np.savetxt(path+'/num_t_'+str(name)+'_mn.txt', t_mns.T)
#np.savetxt(path+'/num_t_'+str(name)+'_na.txt', t_na2s.T)
#np.savetxt(path+'/num_t_'+str(name)+'_zn.txt', t_zns.T)
#np.savetxt(path+'/num_t_'+str(name)+'_kc.txt', t_kcl.T)
#np.savetxt(path+'/num_t_'+str(name)+'_al.txt', t_al2o3.T)


print(count)

east_long = 90. #degrees
west_long = -90. #degrees

west_ti = zeros((nz,nbin))
west_mg = zeros((nz,nbin))
west_fe = zeros((nz,nbin))
west_cr = zeros((nz,nbin))
west_mn = zeros((nz,nbin))
west_na = zeros((nz,nbin))
west_zn = zeros((nz,nbin))
west_kc = zeros((nz,nbin))
west_al = zeros((nz,nbin))

east_ti = zeros((nz,nbin))
east_mg = zeros((nz,nbin))
east_fe = zeros((nz,nbin))
east_cr = zeros((nz,nbin))
east_mn = zeros((nz,nbin))
east_na = zeros((nz,nbin))
east_zn = zeros((nz,nbin))
east_kc = zeros((nz,nbin))
east_al = zeros((nz,nbin))

m_west_ti = zeros((nz,nbin))
m_west_mg = zeros((nz,nbin))
m_west_fe = zeros((nz,nbin))
m_west_cr = zeros((nz,nbin))
m_west_mn = zeros((nz,nbin))
m_west_na = zeros((nz,nbin))
m_west_zn = zeros((nz,nbin))
m_west_kc = zeros((nz,nbin))
m_west_al = zeros((nz,nbin))

m_east_ti = zeros((nz,nbin))
m_east_mg = zeros((nz,nbin))
m_east_fe = zeros((nz,nbin))
m_east_cr = zeros((nz,nbin))
m_east_mn = zeros((nz,nbin))
m_east_na = zeros((nz,nbin))
m_east_zn = zeros((nz,nbin))
m_east_kc = zeros((nz,nbin))
m_east_al = zeros((nz,nbin))

cunt_east = 0
cunt_west = 0

for il in range(ilong):
	if longitude[il] < west_long+20 and longitude[il] > west_long-20:
		cunt_west += 1
		
		west_ti += t_tio2[:,:,il]
		west_mg += t_het_mg2sio4[:,:,il]
		west_fe += t_het_fe[:,:,il] + t_pure_fe[:,:,il]
		west_cr += t_het_cr[:,:,il] + t_cr[:,:,il]
		west_mn += t_mns[:,:,il] 
		west_na += t_na2s[:,:,il]
		west_zn += t_zns[:,:,il]
		west_kc += t_kcl[:,:,il]
		west_al += t_al2o3[:,:,il]
		
		m_west_ti += t_mass_tio2[:,:,il]
		m_west_mg += t_mass_het_mg2sio4[:,:,il]
		m_west_fe += t_mass_het_fe[:,:,il] + t_mass_pure_fe[:,:,il]
		m_west_cr += t_mass_het_cr[:,:,il] + t_mass_cr[:,:,il]
		m_west_mn += t_mass_mns[:,:,il] 
		m_west_na += t_mass_na2s[:,:,il]
		m_west_zn += t_mass_zns[:,:,il]
		m_west_kc += t_mass_kcl[:,:,il]
		m_west_al += t_mass_al2o3[:,:,il]

	
	if longitude[il] < east_long+20 and longitude[il] > east_long-20:
		cunt_east += 1
		
		east_ti += t_tio2[:,:,il]
		east_mg += t_het_mg2sio4[:,:,il]
		east_fe += t_het_fe[:,:,il] + t_pure_fe[:,:,il]
		east_cr += t_het_cr[:,:,il] + t_cr[:,:,il]
		east_mn += t_mns[:,:,il] 
		east_na += t_na2s[:,:,il]
		east_zn += t_zns[:,:,il]
		east_kc += t_kcl[:,:,il]
		east_al += t_al2o3[:,:,il]
		
		m_east_ti += t_mass_tio2[:,:,il]
		m_east_mg += t_mass_het_mg2sio4[:,:,il]
		m_east_fe += t_mass_het_fe[:,:,il] + t_mass_pure_fe[:,:,il]
		m_east_cr += t_mass_het_cr[:,:,il] + t_mass_cr[:,:,il]
		m_east_mn += t_mass_mns[:,:,il] 
		m_east_na += t_mass_na2s[:,:,il]
		m_east_zn += t_mass_zns[:,:,il]
		m_east_kc += t_mass_kcl[:,:,il]
		m_east_al += t_mass_al2o3[:,:,il]
		
west_ti /= cunt_west
west_mg /= cunt_west
west_fe /= cunt_west
west_cr /= cunt_west
west_mn /= cunt_west
west_na /= cunt_west
west_zn /= cunt_west
west_kc /= cunt_west
west_al /= cunt_west

east_ti /= cunt_east
east_mg /= cunt_east
east_fe /= cunt_east
east_cr /= cunt_east
east_mn /= cunt_east
east_na /= cunt_east
east_zn /= cunt_east
east_kc /= cunt_east
east_al /= cunt_east

m_west_ti /= cunt_west
m_west_mg /= cunt_west
m_west_fe /= cunt_west
m_west_cr /= cunt_west
m_west_mn /= cunt_west
m_west_na /= cunt_west
m_west_zn /= cunt_west
m_west_kc /= cunt_west
m_west_al /= cunt_west

m_east_ti /= cunt_east
m_east_mg /= cunt_east
m_east_fe /= cunt_east
m_east_cr /= cunt_east
m_east_mn /= cunt_east
m_east_na /= cunt_east
m_east_zn /= cunt_east
m_east_kc /= cunt_east
m_east_al /= cunt_east

combined_east = np.vstack((east_ti,east_mg,east_fe,east_cr,east_mn,east_na,east_zn,east_kc,east_al))
combined_west = np.vstack((west_ti,west_mg,west_fe,west_cr,west_mn,west_na,west_zn,west_kc,west_al))

np.savetxt(path+'/east_cloud_properties'+str(name)+'_ti.txt', east_ti.T)
np.savetxt(path+'/east_cloud_properties'+str(name)+'_mg.txt', east_mg.T)
np.savetxt(path+'/east_cloud_properties'+str(name)+'_fe.txt', east_fe.T)
np.savetxt(path+'/east_cloud_properties'+str(name)+'_cr.txt', east_cr.T)
np.savetxt(path+'/east_cloud_properties'+str(name)+'_mn.txt', east_mn.T)
np.savetxt(path+'/east_cloud_properties'+str(name)+'_na.txt', east_na.T)
np.savetxt(path+'/east_cloud_properties'+str(name)+'_zn.txt', east_zn.T)
np.savetxt(path+'/east_cloud_properties'+str(name)+'_kc.txt', east_kc.T)
np.savetxt(path+'/east_cloud_properties'+str(name)+'_al.txt', east_al.T)

np.savetxt(path+'/west_cloud_properties'+str(name)+'_ti.txt', west_ti.T)
np.savetxt(path+'/west_cloud_properties'+str(name)+'_mg.txt', west_mg.T)
np.savetxt(path+'/west_cloud_properties'+str(name)+'_fe.txt', west_fe.T)
np.savetxt(path+'/west_cloud_properties'+str(name)+'_cr.txt', west_cr.T)
np.savetxt(path+'/west_cloud_properties'+str(name)+'_mn.txt', west_mn.T)
np.savetxt(path+'/west_cloud_properties'+str(name)+'_na.txt', west_na.T)
np.savetxt(path+'/west_cloud_properties'+str(name)+'_zn.txt', west_zn.T)
np.savetxt(path+'/west_cloud_properties'+str(name)+'_kc.txt', west_kc.T)
np.savetxt(path+'/west_cloud_properties'+str(name)+'_al.txt', west_al.T)

np.savetxt(path+'/mass_east_cloud_properties'+str(name)+'_ti.txt', m_east_ti.T)
np.savetxt(path+'/mass_east_cloud_properties'+str(name)+'_mg.txt', m_east_mg.T)
np.savetxt(path+'/mass_east_cloud_properties'+str(name)+'_fe.txt', m_east_fe.T)
np.savetxt(path+'/mass_east_cloud_properties'+str(name)+'_cr.txt', m_east_cr.T)
np.savetxt(path+'/mass_east_cloud_properties'+str(name)+'_mn.txt', m_east_mn.T)
np.savetxt(path+'/mass_east_cloud_properties'+str(name)+'_na.txt', m_east_na.T)
np.savetxt(path+'/mass_east_cloud_properties'+str(name)+'_zn.txt', m_east_zn.T)
np.savetxt(path+'/mass_east_cloud_properties'+str(name)+'_kc.txt', m_east_kc.T)
np.savetxt(path+'/mass_east_cloud_properties'+str(name)+'_al.txt', m_east_al.T)

np.savetxt(path+'/mass_west_cloud_properties'+str(name)+'_ti.txt', m_west_ti.T)
np.savetxt(path+'/mass_west_cloud_properties'+str(name)+'_mg.txt', m_west_mg.T)
np.savetxt(path+'/mass_west_cloud_properties'+str(name)+'_fe.txt', m_west_fe.T)
np.savetxt(path+'/mass_west_cloud_properties'+str(name)+'_cr.txt', m_west_cr.T)
np.savetxt(path+'/mass_west_cloud_properties'+str(name)+'_mn.txt', m_west_mn.T)
np.savetxt(path+'/mass_west_cloud_properties'+str(name)+'_na.txt', m_west_na.T)
np.savetxt(path+'/mass_west_cloud_properties'+str(name)+'_zn.txt', m_west_zn.T)
np.savetxt(path+'/mass_west_cloud_properties'+str(name)+'_kc.txt', m_west_kc.T)
np.savetxt(path+'/mass_west_cloud_properties'+str(name)+'_al.txt', m_west_al.T)

pressure_plot = 1e-1 #bar
press_index = (np.abs(p/1.e6 - pressure_plot)).argmin()

X1, Y1 = np.meshgrid(longitude,r[:,0]) 
X2, Y2 = np.meshgrid(longitude,r[:,1])
X3, Y3 = np.meshgrid(longitude,r[:,2])
X4, Y4 = np.meshgrid(longitude,r[:,3])
X5, Y5 = np.meshgrid(longitude,r[:,4])
X6, Y6 = np.meshgrid(longitude,r[:,5])
X7, Y7 = np.meshgrid(longitude,r[:,6])
X8, Y8 = np.meshgrid(longitude,r[:,7])
X9, Y9 = np.meshgrid(longitude,r[:,8])
X10, Y10 = np.meshgrid(longitude,r[:,9])
X11, Y11 = np.meshgrid(longitude,r[:,10])

print(shape(t_tio2[press_index,:,:] ))

do_sd_plots = 0

if do_sd_plots:

	for press_index,pp in enumerate(p):
		print(press_index)
		if pp/1.e6 < 1. and pp/1.e6 > 1e-3:
			
			fig, ax = subplots(nrows = 1, ncols = 1)
			im1 = ax.contourf(X1,Y1,t_tio2[press_index,:,:] , norm=LogNorm(),levels=levels_num, cmap='Blues', alpha=0.75)
			im2 = ax.contourf(X2,Y2,t_al2o3[press_index,:,:] , norm=LogNorm(),levels=levels_num, cmap='Greens', alpha=0.75)
			im4 = ax.contourf(X4,Y4,t_het_mg2sio4[press_index,:,:] , norm=LogNorm(),levels=levels_num, cmap='Purples', alpha=0.75)
			im3 = ax.contourf(X3,Y3,t_pure_fe[press_index,:,:] +t_het_fe[press_index,:,:] , norm=LogNorm(),levels=levels_num, cmap='Reds', alpha=0.75)
			im5 = ax.contourf(X5,Y5,t_cr[press_index,:,:] , norm=LogNorm(),levels=levels_num, cmap='Oranges', alpha=0.75)
			im6 = ax.contourf(X6,Y6,t_mns[press_index,:,:] , norm=LogNorm(),levels=levels_num, cmap='copper', alpha=0.75)
			im7 = ax.contourf(X7,Y7,t_na2s[press_index,:,:] , norm=LogNorm(),levels=levels_num, cmap='bone', alpha=0.75)
			im8 = ax.contourf(X10,Y10,t_kcl[press_index,:,:] , norm=LogNorm(),levels=levels_num, cmap='pink', alpha=0.75)
			im9 = ax.contourf(X11,Y11,t_zns[press_index,:,:] , norm=LogNorm(),levels=levels_num, cmap='Greys', alpha=0.75)
			cbar = plt.colorbar(im9)
			cbar.set_label(r'Number Density [cm$^{-3}$]', fontsize=16, rotation=-90, labelpad=30)
			ax.set_yscale('log')
			ax.set_xlim(min(longitude),max(longitude))
			ax.set_ylim(min(r[:,0]),max(r[:,0]))
			ax.set_xlabel(r'Longitude [degree]')
			ax.set_ylabel(r'Particle Radius [$\mu$m]')
			plt.title('Pressure = '+str(p[press_index]/1e6)+' bar')
			plt.savefig(path+'/num_sizedist_'+str(p[press_index]/1e6)+'_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
			#plt.savefig(path+'/num_sizedist_'+str(p[press_index]/1e6)+'_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
			plt.close()
			
			
			
			fig, ax = subplots(nrows = 1, ncols = 1)
			im1 = ax.contourf(X1,Y1,t_mass_tio2[press_index,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Blues', alpha=0.75)
			im2 = ax.contourf(X2,Y2,t_mass_al2o3[press_index,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Greens', alpha=0.75)
			im4 = ax.contourf(X4,Y4,t_mass_het_mg2sio4[press_index,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Purples', alpha=0.75)
			im3 = ax.contourf(X3,Y3,t_mass_pure_fe[press_index,:,:] +t_mass_het_fe[press_index,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Reds', alpha=0.75)
			im5 = ax.contourf(X5,Y5,t_mass_cr[press_index,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Oranges', alpha=0.75)
			im6 = ax.contourf(X6,Y6,t_mass_mns[press_index,:,:] , norm=LogNorm(),levels=levels_mass, cmap='copper', alpha=0.75)
			im7 = ax.contourf(X7,Y7,t_mass_na2s[press_index,:,:] , norm=LogNorm(),levels=levels_mass, cmap='bone', alpha=0.75)
			im8 = ax.contourf(X10,Y10,t_mass_kcl[press_index,:,:] , norm=LogNorm(),levels=levels_mass, cmap='pink', alpha=0.75)
			im9 = ax.contourf(X11,Y11,t_mass_zns[press_index,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Greys', alpha=0.75)
			cbar = plt.colorbar(im9)
			cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=16, rotation=-90, labelpad=30)
			ax.set_yscale('log')
			ax.set_xlim(min(longitude),max(longitude))
			ax.set_ylim(min(r[:,0]),max(r[:,0]))
			ax.set_xlabel(r'Longitude [degree]')
			ax.set_ylabel(r'Particle Radius [$\mu$m]')
			plt.title('Pressure = '+str(p[press_index]/1e6)+' bar')
			plt.savefig(path+'/mass_sizedist_'+str(p[press_index]/1e6)+'_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
			#plt.savefig(path+'/mass_sizedist_'+str(p[press_index]/1e6)+'_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
			plt.close()
			

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

print('Plotting!')


			
e_pressure_level = np.loadtxt(path+'/'+str(name)+'_opticaldepth_level_east.txt', skiprows = 1)
w_pressure_level = np.loadtxt(path+'/'+str(name)+'_opticaldepth_level_west.txt', skiprows = 1)

wave_plevel = e_pressure_level[0]

e_pressure_level_miri = e_pressure_level[1][7]
w_pressure_level_miri = w_pressure_level[1][7]

e_pressure_level_miri_plot = find_nearest(p/1.e6,e_pressure_level_miri)[1]
w_pressure_level_miri_plot = find_nearest(p/1.e6,w_pressure_level_miri)[1]


fig, ax = subplots(nrows = 1, ncols = 1)
im1 = ax.contourf(X1,Y1,t_mass_tio2[e_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Blues', alpha=0.75)
im2 = ax.contourf(X2,Y2,t_mass_al2o3[e_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Greens', alpha=0.75)
im4 = ax.contourf(X4,Y4,t_mass_het_mg2sio4[e_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Purples', alpha=0.75)
im3 = ax.contourf(X3,Y3,t_mass_pure_fe[e_pressure_level_miri_plot,:,:] +t_mass_het_fe[e_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Reds', alpha=0.75)
im5 = ax.contourf(X5,Y5,t_mass_cr[e_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Oranges', alpha=0.75)
im6 = ax.contourf(X6,Y6,t_mass_mns[e_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='copper', alpha=0.75)
im7 = ax.contourf(X7,Y7,t_mass_na2s[e_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='bone', alpha=0.75)
im8 = ax.contourf(X10,Y10,t_mass_kcl[e_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='pink', alpha=0.75)
im9 = ax.contourf(X11,Y11,t_mass_zns[e_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Greys', alpha=0.75)
cbar = plt.colorbar(im9)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=16, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xlim(min(longitude),max(longitude))
ax.set_ylim(min(r[:,0]),max(r[:,0]))
ax.set_xlabel(r'Longitude [degree]')
ax.set_ylabel(r'Particle Radius [$\mu$m]')
plt.title('Pressure = '+str(p[e_pressure_level_miri_plot]/1e6)+' bar')
plt.savefig(path+'/east_miri_tau_mass_sizedist_'+str(p[e_pressure_level_miri_plot]/1e6)+'_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
plt.close()

fig, ax = subplots(nrows = 1, ncols = 1)
im1 = ax.contourf(X1,Y1,t_mass_tio2[w_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Blues', alpha=0.75)
im2 = ax.contourf(X2,Y2,t_mass_al2o3[w_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Greens', alpha=0.75)
im4 = ax.contourf(X4,Y4,t_mass_het_mg2sio4[w_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Purples', alpha=0.75)
im3 = ax.contourf(X3,Y3,t_mass_pure_fe[w_pressure_level_miri_plot,:,:] +t_mass_het_fe[w_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Reds', alpha=0.75)
im5 = ax.contourf(X5,Y5,t_mass_cr[w_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Oranges', alpha=0.75)
im6 = ax.contourf(X6,Y6,t_mass_mns[w_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='copper', alpha=0.75)
im7 = ax.contourf(X7,Y7,t_mass_na2s[w_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='bone', alpha=0.75)
im8 = ax.contourf(X10,Y10,t_mass_kcl[w_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='pink', alpha=0.75)
im9 = ax.contourf(X11,Y11,t_mass_zns[w_pressure_level_miri_plot,:,:] , norm=LogNorm(),levels=levels_mass, cmap='Greys', alpha=0.75)
cbar = plt.colorbar(im9)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=16, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xlim(min(longitude),max(longitude))
ax.set_ylim(min(r[:,0]),max(r[:,0]))
ax.set_xlabel(r'Longitude [degree]')
ax.set_ylabel(r'Particle Radius [$\mu$m]')
plt.title('Pressure = '+str(p[e_pressure_level_miri_plot]/1e6)+' bar')
plt.savefig(path+'/west_miri_tau_mass_sizedist_'+str(p[e_pressure_level_miri_plot]/1e6)+'_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
plt.close()


	

X1, Y1 = np.meshgrid(longitude,p/1.e6)
X2, Y2 = np.meshgrid(longitude,p/1.e6)
X3, Y3 = np.meshgrid(longitude,p/1.e6)
X4, Y4 = np.meshgrid(longitude,p/1.e6)
X5, Y5 = np.meshgrid(longitude,p/1.e6)
X6, Y6 = np.meshgrid(longitude,p/1.e6)
X7, Y7 = np.meshgrid(longitude,p/1.e6)
X8, Y8 = np.meshgrid(longitude,p/1.e6)
X9, Y9 = np.meshgrid(longitude,p/1.e6)
X10, Y10 = np.meshgrid(longitude,p/1.e6)
X11, Y11 = np.meshgrid(longitude,p/1.e6)

fig, ax = subplots(nrows = 1, ncols = 1)
im1 = ax.contourf(X1,Y1,num_density_ti[:,:] , norm=LogNorm(),levels=levels_num*10, cmap='Blues', alpha=0.75)
im2 = ax.contourf(X2,Y2,num_density_al[:,:] , norm=LogNorm(),levels=levels_num*10, cmap='Greens', alpha=0.75)
im4 = ax.contourf(X4,Y4,num_density_mg[:,:] , norm=LogNorm(),levels=levels_num*10, cmap='Purples', alpha=0.75)
im3 = ax.contourf(X3,Y3,num_density_fe[:,:] , norm=LogNorm(),levels=levels_num*10, cmap='Reds', alpha=0.75)
im5 = ax.contourf(X5,Y5,num_density_cr[:,:] , norm=LogNorm(),levels=levels_num*10, cmap='Oranges', alpha=0.75)
im6 = ax.contourf(X6,Y6,num_density_mn[:,:] , norm=LogNorm(),levels=levels_num*10, cmap='copper', alpha=0.75)
im7 = ax.contourf(X7,Y7,num_density_na[:,:] , norm=LogNorm(),levels=levels_num*10, cmap='bone', alpha=0.75)
im8 = ax.contourf(X10,Y10,num_density_kc[:,:] , norm=LogNorm(),levels=levels_num*10, cmap='pink', alpha=0.75)
im9 = ax.contourf(X11,Y11,num_density_zn[:,:] , norm=LogNorm(),levels=levels_num*10, cmap='Greys', alpha=0.75)
cbar = plt.colorbar(im9)
cbar.set_label(r'Number Density [cm$^{-3}$]', fontsize=20, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xlim(min(longitude),max(longitude))
ax.set_ylim(1,min(p/1.e6))
ax.set_xlabel(r'Longitude [degree]')
ax.set_ylabel(r'Pressure [bar]')
plt.savefig(path+'/num_altitude_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
#plt.savefig(path+'/num_altitude_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
plt.close()

fig, ax = subplots(nrows = 1, ncols = 1)
im1 = ax.contourf(X1,Y1,mass_density_ti[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
im2 = ax.contourf(X2,Y2,mass_density_al[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
im4 = ax.contourf(X4,Y4,mass_density_mg[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
im3 = ax.contourf(X3,Y3,mass_density_fe[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
im5 = ax.contourf(X5,Y5,mass_density_cr[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
im6 = ax.contourf(X6,Y6,mass_density_mn[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='copper', alpha=0.75)
im7 = ax.contourf(X7,Y7,mass_density_na[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='bone', alpha=0.75)
im8 = ax.contourf(X10,Y10,mass_density_kc[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='pink', alpha=0.75)
im9 = ax.contourf(X11,Y11,mass_density_zn[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greys', alpha=0.75)
cbar = plt.colorbar(im9)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=20, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xlim(min(longitude),max(longitude))
ax.set_ylim(1,min(p/1.e6))
ax.set_xlabel(r'Longitude [degree]')
ax.set_ylabel(r'Pressure [bar]')
plt.savefig(path+'/mass_altitude_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
#plt.savefig(path+'/mass_altitude_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
plt.close()

fig, ax = subplots(nrows = 1, ncols = 1)
im1 = ax.contourf(X1,Y1,mass_density_ti[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Blues', alpha=0.75)
cbar = plt.colorbar(im1)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=20, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xlim(min(longitude),max(longitude))
ax.set_ylim(1,min(p/1.e6))
ax.set_xlabel(r'Longitude [degree]')
ax.set_ylabel(r'Pressure [bar]')
plt.savefig(path+'/mass_altitude_ti_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
#plt.savefig(path+'/mass_altitude_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
plt.close()

fig, ax = subplots(nrows = 1, ncols = 1)
im2 = ax.contourf(X2,Y2,mass_density_al[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Greens', alpha=0.75)
cbar = plt.colorbar(im2)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=20, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xlim(min(longitude),max(longitude))
ax.set_ylim(1,min(p/1.e6))
ax.set_xlabel(r'Longitude [degree]')
ax.set_ylabel(r'Pressure [bar]')
plt.savefig(path+'/mass_altitude_al_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
#plt.savefig(path+'/mass_altitude_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
plt.close()

fig, ax = subplots(nrows = 1, ncols = 1)
im4 = ax.contourf(X4,Y4,mass_density_mg[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Purples', alpha=0.75)
cbar = plt.colorbar(im4)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=20, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xlim(min(longitude),max(longitude))
ax.set_ylim(1,min(p/1.e6))
ax.set_xlabel(r'Longitude [degree]')
ax.set_ylabel(r'Pressure [bar]')
plt.savefig(path+'/mass_altitude_mg_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
#plt.savefig(path+'/mass_altitude_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
plt.close()

fig, ax = subplots(nrows = 1, ncols = 1)
im3 = ax.contourf(X3,Y3,mass_density_fe[:,:], norm=LogNorm(),levels=levels_mass*10, cmap='Reds', alpha=0.75)
cbar = plt.colorbar(im3)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=20, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xlim(min(longitude),max(longitude))
ax.set_ylim(1,min(p/1.e6))
ax.set_xlabel(r'Longitude [degree]')
ax.set_ylabel(r'Pressure [bar]')
plt.savefig(path+'/mass_altitude_fe_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
#plt.savefig(path+'/mass_altitude_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
plt.close()

fig, ax = subplots(nrows = 1, ncols = 1)
im5 = ax.contourf(X5,Y5,mass_density_cr[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Oranges', alpha=0.75)
cbar = plt.colorbar(im5)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=20, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xlim(min(longitude),max(longitude))
ax.set_ylim(1,min(p/1.e6))
ax.set_xlabel(r'Longitude [degree]')
ax.set_ylabel(r'Pressure [bar]')
plt.savefig(path+'/mass_altitude_cr_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
#plt.savefig(path+'/mass_altitude_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
plt.close()

fig, ax = subplots(nrows = 1, ncols = 1)
im5 = ax.contourf(X5,Y5,mass_density_mn[:,:],norm=LogNorm(),levels=levels_mass*10, cmap='Copper', alpha=0.75)
cbar = plt.colorbar(im5)
cbar.set_label(r'Mass Density [g cm$^{-3}$]', fontsize=20, rotation=-90, labelpad=30)
ax.set_yscale('log')
ax.set_xlim(min(longitude),max(longitude))
ax.set_ylim(1,min(p/1.e6))
ax.set_xlabel(r'Longitude [degree]')
ax.set_ylabel(r'Pressure [bar]')
plt.savefig(path+'/mass_altitude_mn_'+str(name)+'.jpg', format = 'jpg', bbox_inches='tight')
#plt.savefig(path+'/mass_altitude_'+str(name)+'.eps', format = 'eps', bbox_inches='tight')
plt.close()

'''