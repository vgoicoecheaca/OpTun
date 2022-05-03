'''
Used to calculate the transmisison and refleciton coefficients of 
the Acrylic, ITO system 
Credits to TMM software, github sbyrnes321/tmm

'''
import numpy as np
import matplotlib.pyplot as plt
import tmm
from scipy.interpolate import interp1d
import uproot as ur
import awkward as ak

plt.style.use('mystyle.mlstyle')

cali                 = False
other                = False
degree               = 180/np.pi
rad                  = np.pi/180
lAr_n,lAr_t          = 1.233, np.inf
gAr_n, gAr_t         = 1.008, np.inf
ito_t, ito_n         = 15,  2.08+1*-0.04j
acrylic_t, acrylic_n = 4.5*1000000,  1.4889      # 4.5mm
tpb_t,tpb_n          = 10, 1.5                   # Unkown? Taking that of g4ds
nsamples             = 200
ls                   = np.linspace(300,600,nsamples)
ks                   = (2*np.pi)/ls
c_i                  = ["i","c","i","c","i"]

#read the indices of refraction from ITO, interpolate a funciton for later use at different wavelengths
fito                     = np.loadtxt("ito.txt",skiprows=1)
ito_ns                   = np.zeros((len(fito),2)).astype(complex)
ito_ns[:,0],ito_ns[:,1]  = fito[:,0], fito[:,1] + fito[:,2]*1j
ito_nk_fn                = interp1d(ito_ns[:,0].real,ito_ns[:,1],kind='quadratic')

fsilicia                 = np.loadtxt("silicia.txt",skiprows=1,delimiter=",")
silicia_n_fn             = interp1d(fsilicia[:,0]*1000,fsilicia[:,1],kind='linear') #convert to nm

facrylic                 = np.loadtxt("acrylic.txt",skiprows=1)
acrylic_n_fn             = interp1d(facrylic[:,0],facrylic[:,1],kind='linear')


#grid_trans              = np.loadtxt("gridtrans.txt")
#plt.plot(grid_trans[:,0],grid_trans[:,1])
#plt.ylabel("T [\%]")
#plt.xlabel("Angle of Incidence [deg]")
#plt.grid(True)
#plt.show()
#import sys
#sys.exit()


#check the indices of refraciton interpolations
#plt.plot(ls,ito_nk_fn(ls).real,label = "ITO - refractive index",      linestyle = "solid",  color="blue")
#plt.plot(ls[ls>370],acrylic_n_fn(ls[ls>370]),label = "Acrylic - refractive index",  linestyle = "dashdot",color="blue")
#plt.plot(ls,silicia_n_fn(ls),label = "Silicia - refractive index",  linestyle = ":", color="blue")
#plt.plot(ls,ito_nk_fn(ls).imag,label = "ITO - exctinction coefficient",linestyle = "solid",  color="red")
#plt.xlabel("Wavelength [nm]")
#plt.legend(loc='best')
#plt.grid(True)
#plt.show()

################################################
#### Def funcitons for g4ds output data ########
################################################

#getting all information needed from direction (theta0) and position of deposit (for R,T)
fname           = "ITO_TRA.root" 
nevents         = 10000000
branches_names  = ["ph_z","pz","px"]
branches        = {} 

#for branch in branches_names:
#    branches[branch] = ur.open(fname)["dstree"][branch].array(library="ak")[:nevents]
#
##ph_pid_maxs = ak.max(branches["ph_pid"],axis=1)
#ph_z =  branches["ph_z"]
#
#def get_tra_mc(ph_z,px,pz):       
#    angles   = np.linspace(0,1.57,200)
#    thetas   = np.arctan(px/pz) 
#    bin_idxs = np.digitize(thetas,angles)   
#    rato     =  np.zeros((len(angles),4))
#    for i in range(len(angles)):
#        idxs = np.where(bin_idxs==i)[0]                              #indices of events that fall in this bin
#        if len(ph_z[idxs])!=0:
#            #zs = ph_z[idxs,-1]
#            zs = ak.flatten(ph_z[idxs])
#            rato[i,0] = len(zs[zs<3.505]) / len(zs)
#            rato[i,1] = len(zs[(zs>3.505) & (zs<3.605)]) / len(zs)
#            rato[i,2] = len(zs[zs>3.605]) / len(zs)
#    return rato, angles
#
#rato, angles_mc = get_tra_mc(ph_z,branches["px"],branches["pz"])
#r_mc,a_mc,t_mc,o_mc = rato[:,0] ,rato[:,1], rato[:,2] ,rato[:,3] 

#############################################
### Transmission, Reflection, Absorption ####
#############################################

if cali:
    print("Validation")
    ts = [np.inf,15,np.inf]
    ns = [1.47, ito_n,1.233]
    c_i = ["i","c","i"]
    i_ito = 1
else:
    ts = [gAr_t,ito_t,acrylic_t, ito_t, lAr_t]
    ns = [gAr_n,ito_n,acrylic_n, ito_n, lAr_n]
    c_i = ["i","c","i","c","i"]
    ts  = [gAr_t,tpb_n,ito_t,acrylic_t, lAr_t]
    ns  = [gAr_n,tpb_n,ito_n,acrylic_n, lAr_n]

def Lambda_normal():
    ls = np.linspace(300,550,nsamples) 
    ks = (2*np.pi)/ls
    #first at normal incidence
    R,T,A  = [0 for i in range(len(ls))],[0 for i in range(len(ls))],[0 for i in range(len(ls))]
    for i in range(len(ls)): 
        if cali:
            ns[i_ito] = ito_nk_fn(ls[i])                                                      #update the index of refraciton if ITO
        else:
            ns[1],ns[3] = ito_nk_fn(ls[i]) ,ito_nk_fn(ls[i])                               #update the index of refraciton if ITO
            ns[2] = acrylic_n if ls[i]<370 else acrylic_n_fn(ls[i])
        data = tmm.unpolarized_RT_inc(ns, ts, c_i, 0, ls[i])
        R[i] = data["R"]
        T[i] = data["T"]   
        A[i] = 1 - R[i] - T[i]

    plt.plot(ls, R, 'green',label="R") 
    plt.plot(ls, T, "blue",label="T")
    plt.plot(ls, A, "red", label="A")
    plt.xlabel(r'$\mathrm{\lambda} \mathrm{[nm]}$', loc='center')
    plt.ylabel("T,R,A [\%]",loc='center')
    plt.legend(loc='best')
    plt.grid(True)
    plt.ylim(0,1)
    plt.xlim(300,550)
    plt.show()

#Lambda_normal()

def ITOts():
    wl        = 420
    k         = 1/wl         
    R,T,A     = [0 for i in range(len(ls))],[0 for i in range(len(ls))],[0 for i in range(len(ls))]
    sample_ts = np.linspace(1,500,nsamples) 
    #ns[2]     = silicia_n_fn(wl) if cali else acrylic_n_fn(wl) 
    for i in range(len(sample_ts)):
        if cali:
            ns[i_ito] = ito_nk_fn(ls[i])
            ts[1] = sample_ts[i]
        else:
            ts[1],ts[3] = sample_ts[i] , sample_ts[i]
            ns[1],ns[3] = ito_nk_fn(ls[i]) ,ito_nk_fn(ls[i])
        data = tmm.unpolarized_RT_inc(ns, ts, c_i, 0, wl)
        R[i]  = data["R"]*100
        T[i]  = data["T"]*100
        A[i]  = (100 - R[i] - T[i])
 
    ts[1] = ito_t              # give back the ITO thicnkess value
    if not cali:
        ts[3] = ito_t
    plt.plot(sample_ts, R, 'green',label="R") 
    plt.plot(sample_ts, T, "blue",label="T")
    plt.plot(sample_ts, A, "red", label="A")
    plt.xlabel("ITO Thickness [nm]", loc='center')
    plt.ylabel("T,R,A [\%]",loc='center')
    plt.legend(loc='best')
    plt.grid(True)
    plt.ylim(0,100)
    plt.xlim(0,sample_ts[-1])
    plt.show()

#ITOts()

def angle():
    wl     = 420
    #ns[2]  = silicia_n_fn(wl) if cali else acrylic_n_fn(wl) 
    k = 1/wl       
    sample_angles = np.linspace(0,1.57,nsamples) 
    R,T,A  = [0 for i in range(len(ls))],[0 for i in range(len(ls))],[0 for i in range(len(ls))]
    for i in range(len(sample_angles)):
        if not cali:
            ns[1],ns[3] = ito_nk_fn(ls[i]) ,ito_nk_fn(ls[i])
            ns[2]  = acrylic_n_fn(wl) 
        else:
            ns[i_ito] = ito_nk_fn(ls[i])
        data = tmm.unpolarized_RT_inc(ns, ts, c_i, sample_angles[i], wl)
        R[i]  = data['R']*100
        T[i]  = data['T']*100
        A[i] = 100 - R[i] - T[i]

    #print(sample_angles[11],r_mc[11]) 
    plt.plot(sample_angles, R, 'green',label="R") 
    plt.plot(sample_angles, T, "blue",label="T")
    plt.plot(sample_angles, A, "red", label="A")
    #plt.plot(angles_mc[r_mc!=0],r_mc[r_mc!=0], "xg")
    #plt.plot(angles_mc[t_mc!=0],t_mc[t_mc!=0], "xb")
    #plt.plot(angles_mc[a_mc!=0],a_mc[a_mc!=0], "xr")
    plt.xlabel("Incident Angle [rad]", loc='center')
    plt.ylabel("T,R,A [\%]",loc='center')
    plt.xlim(0,sample_angles[-1])
    plt.ylim(0,100)
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()

angle()

def acrylic():
    for wl in [128,420]:
        #wl     = 128
        sample_angles = np.linspace(0,1.57,nsamples) 
        R,T,A  = [0 for i in range(len(ls))],[0 for i in range(len(ls))],[0 for i in range(len(ls))]
        for i in range(len(sample_angles)):
            data = tmm.unpolarized_RT_inc([1.4,1.5,1], [np.inf,1,np.inf], ["i","c","i"], sample_angles[i], wl)
            R[i]  = data['R']*100
            T[i]  = data['T']*100
            A[i] = 100 - R[i] - T[i]
#
        plt.plot(sample_angles, R, 'g' if wl ==420 else "xg",label="R - Vis" if wl ==420 else "R - UV") 
        plt.plot(sample_angles, T, "b" if wl ==420 else "xb",label="T - Vis" if wl ==420 else "T - UV")
        plt.plot(sample_angles, A, "r" if wl ==420 else "xr", label="A - Vis" if wl ==420 else "A - UV")

    #plt.plot(angles_mc[r_mc!=0],r_mc[r_mc!=0]*100, "xg")
    #plt.plot(angles_mc[t_mc!=0],t_mc[t_mc!=0]*100, "xb")
    #plt.plot(angles_mc[a_mc!=0],a_mc[a_mc!=0]*100, "xr")
    plt.xlabel("Incident Angle [rad]", loc='center')
    plt.ylabel("T,R,A [\%]",loc='center')
    plt.xlim(0,sample_angles[-1])
    plt.ylim(0,100)
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()

#acrylic()

# for TPB-Acrylic with LAr on both sides
ns = [1.23,1.22,1.48,1.23]
ts = [np.inf,10,4.5e6,np.inf]
cs = ["i","c","i","i"]

def custom_layers(ns,ts,cs,wl = 420):
    sample_angles = np.linspace(0,1.57,nsamples) 
    R,T,A  = [0 for i in range(len(ls))],[0 for i in range(len(ls))],[0 for i in range(len(ls))]
    for i in range(len(sample_angles)):
        data = tmm.unpolarized_RT_inc(ns, ts, cs,  sample_angles[i], wl)
        R[i]  = data['R']*100
        T[i]  = data['T']*100
        A[i] = 100 - R[i] - T[i]
    plt.plot(sample_angles, R, 'g' if wl ==420 else "xg",label="R - Vis" if wl ==420 else "R - UV") 
    plt.plot(sample_angles, T, "b" if wl ==420 else "xb",label="T - Vis" if wl ==420 else "T - UV")
    plt.plot(sample_angles, A, "r" if wl ==420 else "xr", label="A - Vis" if wl ==420 else "A - UV")
    plt.xlabel("Incident Angle [rad]", loc='center')
    plt.ylabel("T,R,A [\%]",loc='center')
    plt.xlim(0,sample_angles[-1])
    plt.ylim(0,100)
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()

custom_layers(ns,ts,cs)