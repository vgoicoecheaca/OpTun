from tkinter import W
import uproot as ur
import numpy as np
import ROOT as R
from array import array

class Helper():
    def __init__(self,manager):
        self.nevents = manager.config("run","nevents","int")
        self.channel_map          = manager.config('channels','channel_na','str').split(',')
        self.channel_config       = manager.config('channels',"channel_config","str").split(',')

    def read_branch(self,fname,branch):       
        return ur.open(fname)["dstree"][branch].array(library="np")[:self.nevents]
     
    def get_s1ns2(self,cl_q,ncl_s1,ncl_s2):
        print("Getting S1 & S2")
        cut_ncl = (ncl_s2 == 1) & (ncl_s1 == 1)

        return cl_q[cut_ncl][0], cl_q[cut_ncl][1]

    def get_tdrift(self,ncl,cl_startt,):
        print("Getting tdrift")

        return (cl_startt[ncl>1][1] - cl_startt[ncl>1][0])*2e-3

    def map_channel(self,event,chs):
          nx_sipms, ny_sipms = 4,6
          hits_top = np.zeros((nx_sipms*ny_sipms),dtype = np.float32)
          hits_bot = np.zeros((4),dtype = np.float32)
          for i in range(len(chs)):
            if chs[i][0] == 'F':
              hits_bot[self.channel_config.index(chs[i])] = event[i]
            #if chs[i] in self.channel_map and chs[i][0] != "F":
            else:
              hits_top[self.channel_config.index(chs[i])-4] = event[i] 
            #else:
              #hits_top[i] = 0  
          return hits_top, hits_bot

    def get_zlyield(self,pars,zs,s1):    
      ''' Only works for Kr83m peak now
      '''
      nbins, half_z_bin = pars["bins"],(pars["range"][1]-pars["range"][0])/(2*pars["bins"])
      hz = R.TH1F("hz","",pars["bins"],pars["range"][0],pars["range"][1])
      for i in range(len(zs)): hz.Fill(zs[i])
      fit_pars = np.zeros((nbins,4))        
      z_bins = np.arange(-2.5,2.5,5/nbins)        
      for i,z in enumerate(z_bins):        
          ul,bl = z + half_z_bin, z-half_z_bin        
          hs1 = R.TH1F("hs1","",pars["bins"],pars["fmin"],pars["fmax"])        
          s1_z = s1[(zs>bl) & (zs<ul)]        
          if len(s1_z) != 0:        
              for s in range(len(s1_z)): hs1.Fill(s1_z[s])        
              hs1.Fit("gaus","R","",pars["fmin"],pars["fmax"])         
              fg = hs1.GetFunction('gaus')        
              try:        
                  fit_pars[i,0],fit_pars[i,1],fit_pars[i,2]= fg.GetParameter(0), fg.GetParameter(1), fg.GetParameter(2)        
                  fit_pars[i,-1] = len(s1_z)
              except TypeError:        
                  pass        
              hs1.Reset("ICE")       
      s1lyields = fit_pars[:,1] / (43)
      s1lyields_err = fit_pars[:,2] / (43*fit_pars[:,-1])
      self.sigma_over_mu = np.abs(np.average(fit_pars[2:,2]/fit_pars[2:,1]))
      print("Sigma_over_mu", self.sigma_over_mu)
      zs = np.abs(np.linspace(pars["range"][0],pars["range"][1],nbins))/1.25
      zeros = [0 for i in range(nbins)]
      
      return s1lyields/s1lyields[12],zs,zeros,s1lyields_err/s1lyields[12]

    def detector_response(self,s1,pars):
      print("Hello there")
      mu_over_sigma = 0.064
      #mu_over_sigma = 950/300 # from S1 in plot
      sigma         = 9.8/mu_over_sigma
      sigma = 300
      s1_res   = s1.Clone("hr_res")
      s1_res.Reset()
      f_smear  = R.TF1('f1_smear',	'gaus')
      f_res    = R.TF1("f1_res","0.009 + 0.485/sqrt(x)", pars["range"][0], pars["range"][1])
      for i in range(pars["bins"]):
        bincenter  = s1.GetBinCenter(i)
        bincontent = s1.GetBinContent(i)
        binres     = f_res.Eval(bincenter)
        f_smear.SetParameters(1,bincenter,sigma)
        #f_smear.SetParameters(1,bincenter,bincenter*binres)
        s1_res.FillRandom("f1_smear", int(bincontent))

      return s1_res

    def normalize_hists(self,run,mc,bins,range):
      hr, br = np.histogram(run, bins, range)
      hm, _ = np.histogram(mc, bins, range)
      bc = br[:-1] + np.diff(br) / 2
      hr, hm = hr / np.sum(hr), hm / np.sum(hm)

      return hr, hm, bc 

    def descrete_mean_and_unc(self,edges,xdata,ydata):
        d = np.digitize(xdata, edges)
        means    = np.zeros(len(edges))
        stat_unc = np.zeros(len(edges))
        for i in range(len(edges)):
                means[i]    = np.average(ydata[d==i])
                stat_unc[i] = np.std(ydata[d==i])/np.sqrt(len(ydata[d==i]))  
                #stat_unc[i] = np.std(ydata[d==i]))
        return edges, means, stat_unc 
        #return edges[np.isfinite(means).astype(int)], means[np.isfinite(means).astype(int)], stat_unc[np.isfinite(means).astype(int)]
        #return edges[(np.isfinite(means)) & (np.isfinite(stat_unc))], means[(np.isfinite(means)) & (np.isfinite(stat_unc))], stat_unc[(np.isfinite(means)) & (np.isfinite(stat_unc))]
 
    def reorder(self,array,order):
      out = np.zeros(len(array))
      for i,idx in enumerate(order):
        out[i] = array[idx]     
      
      return out









