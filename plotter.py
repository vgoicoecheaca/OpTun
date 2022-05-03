import matplotlib.pyplot as plt
from branches import Branches
import numpy as np
from cmocean.cm import dense
import matplotlib as mpl
from scipy.interpolate import interp1d
import ROOT as R
from array import array

plt.style.use('mystyle.mlstyle')

class Plotter():
    def __init__(self,manager):
        self.m = manager

    def get_branches(self,branches): 
        self.branches = branches

    def check_s1s2(self,d,mode="s"):
        if "S1" not in self.m.config("run","mc","str")  and "S2" not in self.m.config("run","run","mc"):
            print("Change name of MC file, no S1 or S2")
            import sys
            sys.exit()
        val = False
        if mode != "s":
            for data in d:
               if "s2" in data and "S2" not in self.m.config("run","mc","str"):
                   val =  True
        else:
            if "s2" in d and "S2" not in self.m.config("run","mc","str"):
                val = True
        return val

    def oneDplotter(self,**pars): 
        d_s = pars["d"]
        if self.check_s1s2(d_s):
            return None
        d = self.branches[d_s] if "cl" not in d_s else self.branches[d_s][:,0]
        c = R.TCanvas()
        hr  = R.TH1F("hr","",pars["bins"],pars["range"][0],pars["range"][1])
        for i in range(len(d)): hr.Fill(d[i])
        if "fit" in pars.keys() and pars["fit"]=="gaus":
                binmax    = hr.GetMaximumBin()
                bincenter = hr.GetXaxis().GetBinCenter(binmax)
                fg        = R.TF1("fg","gaus",bincenter-3*np.sqrt(bincenter),bincenter+3*np.sqrt(bincenter))
                hr.Fit("fg","R")
                fg.Draw("same")
        if "d_mc" in pars.keys() and pars["d_mc"]:        
            d_s_mc = d_s.replace("run","mc")
            hmc = R.TH1F("hmc","",pars["bins"],pars["range"][0],pars["range"][1])
            for i in range(len(self.branches[d_s_mc])): hmc.Fill(self.branches[d_s_mc][i])
            #if "d_res" in pars.keys() and pars["d_res"]:
                #hmc = self.m.helper.detector_response(hmc,pars)
            hmc.SetLineColor(2)
            hmc.Scale(hr.Integral()/hmc.Integral())
            hmc.SetMarkerColor(2)
            hmc.Draw("hist")
            leg = R.TLegend(0.55,0.75,0.65,0.85)
            leg.AddEntry(hmc,"MC","l")
            leg.AddEntry(hr,"Run","lep")
            #leg.Draw()
            hmc.GetXaxis().SetTitle(pars["xlabel"])
            hmc.GetYaxis().SetTitle(pars["ylabel"])        
            #hr.Scale(1./hr.Integral(),"width")
            ymaxr = hr.GetMaximum() + hr.GetMaximum()*0.1
            hmc.SetAxisRange(0,ymaxr*1.2,"Y")
            chi2 = hr.Chi2Test(hmc,"UW CHI2")
            tlx,tly = pars["range"][1]*0.55, ymaxr*0.9
            tl = R.TLatex(tlx,tly,"\chi^{2} = "+str(round(chi2,2)))
            leg.AddEntry(tl,"","")
            leg.Draw()
        hr.Sumw2()
        hr.Draw("same")
        if "fmin" in pars.keys():
            hr.Fit("gaus","","",pars["fmin"],pars["fmax"])
            fg = hr.GetFunction('gaus')
            fg.SetLineColor(6)
            fg.SetLineWidth(1)
            pconst, pmu, psigma = fg.GetParameter(0), fg.GetParameter(1), fg.GetParameter(2)               
            fg.Draw("same")
        hr.GetXaxis().SetTitle(pars["xlabel"])
        hr.GetYaxis().SetTitle(pars["ylabel"])        
        if "grid" in pars.keys() and pars["grid"]:
            c.SetGrid()
        c.Update()
        c.SaveAs("plots/"+pars["title"]+".pdf") 
        c.Clear()
 
    def twoDplotter(self,**pars):
        d_s = pars["d"].split(",")
        if self.check_s1s2(d_s,mode="l"):
            return None
        d = [0,0] 
        for i,var in enumerate(d_s):                                    #only works for s1 this way, will fix later
            d[i] = self.branches[d_s[i]]
            if "cl" == var[:2]:
                d[i] = self.branches[d_s[i]][:,0]          
        c = R.TCanvas()
        hr = R.TH2F("hr","",pars["bins"],pars["range"][0][0],pars["range"][0][1],pars["bins"],pars["range"][1][0],pars["range"][1][1])
        for i in range(len(d[0])): hr.Fill(d[0][i],d[1][i] ) 
        hr.Draw("COLZ")  
        R.gStyle.SetPadRightMargin(0.16)
        hr.GetZaxis().SetTitle("Events")
        hr.GetXaxis().SetTitle(pars["xlabel"])
        hr.GetYaxis().SetTitle(pars["ylabel"])        
        c.SaveAs("plots/"+pars["title"]+".pdf") 
        c.Clear()

    def lyield_z(self,**pars):
        ''' Only for kr83m right now, if binning changes, change central bin for both mc and run
        '''
        d_s = pars["d"].split(",")
        s1 = np.array(self.branches[d_s[0]])
        zs  = np.array(self.branches[d_s[1]])
        s1lyields, z, zeros, s1lyields_err = self.m.helper.get_zlyield(pars,zs,s1)
        c = R.TCanvas("c","",700,900)         
        if "d_mc" in pars.keys() and pars["d_mc"]:        
            zs_mc , s1lyields_mc, s1lyields_err_mc = self.m.helper.descrete_mean_and_unc(z,self.branches["z_mc"],self.branches["s1npe_mc"]/self.branches["s1ene_mc"])   
            ratio = R.TGraph(pars["bins"])
            for i in range(pars["bins"]): ratio.SetPoint(i,z[i],s1lyields[i]/s1lyields_mc[i])
            z,s1lyields,  s1lyields_err           = array("f",z) ,   array("f",s1lyields), array("f",s1lyields_err)
            zeros, s1lyields_mc, s1lyields_err_mc = array("f",zeros),array("f",s1lyields_mc/s1lyields_mc[12]), array("f",s1lyields_err_mc/s1lyields_mc[12])
            p1,p2 = R.TPad("p1","p1",0.,0.3,1.,1.), R.TPad("p2","p2",0.,0.,1.,0.3)
            p2.Draw()
            p1.Draw()
            p2.SetTopMargin(0.001)
            p2.SetBottomMargin(0.3)
            p2.SetGrid()
            p1.SetBottomMargin(0.001)
            p1.cd()
            p1.SetGrid()
            gr_run = R.TGraphErrors(pars["bins"],z,s1lyields,zeros,s1lyields_err)
            gr_mc = R.TGraphErrors(pars["bins"],z,s1lyields_mc,zeros,s1lyields_err_mc)   
            gr_run.GetXaxis().SetTitle(pars["xlabel"])
            gr_run.GetYaxis().SetTitle(pars["ylabel"])           
            gr_run.GetYaxis().SetRangeUser(0.5,3)
            gr_run.Draw("AP")
            gr_mc.SetMarkerColor(2)
            gr_mc.SetLineColor(2)
            gr_mc.SetFillColor(2)
            gr_mc.Draw("P") 
            leg = R.TLegend(0.65,0.75,0.75,0.85)
            leg.AddEntry(gr_mc,"MC","lep")
            leg.AddEntry(gr_run,"Run","lep")
            leg.Draw()
            p2.cd()
            ratio.GetXaxis().SetTitle(pars["xlabel"])
            ratio.Draw("AP")
        else:
            gr_run = R.TGraphErrors(pars["bins"],z,s1lyields,zeros,s1lyields_err)
            gr_run.GetXaxis().SetTitle(pars["xlabel"])
            gr_run.GetYaxis().SetTitle(pars["ylabel"])        
            gr_run.Draw("ap")
            gr_run.GetXaxis().SetRangeUser(0,2.5)
        
        c.Update()
        c.SaveAs("plots/"+pars["title"]+".pdf") 
        c.Clear()


    def profiler(self,**pars):
        d_s = pars["d"].split(",")
        if self.check_s1s2(d_s,mode="l"):
            return None
        hr  = R.TProfile("hr" ,"",pars["bins"],pars["range"][0][0],pars["range"][0][1],pars["range"][1][0],pars["range"][1][1])
        hmc = R.TProfile("hmc","",pars["bins"],pars["range"][0][0],pars["range"][0][1],pars["range"][1][0],pars["range"][1][1])
        for i in range(len(self.branches[d_s[0]])): 
            hr.Fill(self.branches[d_s[0]][i],self.branches[d_s[1]][i]) 
        for i in range(len(self.branches[d_s[0].replace("run","mc")])): 
            hmc.Fill(self.branches[d_s[0].replace("run","mc")][i],self.branches[d_s[1].replace("run","mc")][i]) 

        c  = R.TCanvas()
        hr.SetMarkerStyle(7)
        hr.SetLineColor(1)
        hr.SetMarkerColor(1)
        hr.GetXaxis().SetTitle(pars["xlabel"])
        hr.GetYaxis().SetTitle(pars["ylabel"])
        hmc.SetLineColor(2)
        hmc.SetMarkerColor(2)
        hmc.SetMarkerStyle(7)
        if "grid" in pars.keys() and pars["grid"]:
            c.SetGrid()
        #ratio = R.TRatioPlot(hmc,hr)
        #ratio.Draw()
        #ratio.GetUpperPad().Clear()
        hr.Draw("PE")
        hmc.Draw("PE,same")
        hmc.GetYaxis().SetRangeUser(-0.3,0.3)
        hr.GetYaxis().SetRangeUser(-0.3,0.3)
        chi2 = hr.Chi2Test(hmc,"UW CHI2")
        leg = R.TLegend(0.55,0.75,0.65,0.85)
        leg.AddEntry(hmc,"MC","lep")
        leg.AddEntry(hr,"Run","lep")
        leg.SetTextSize(0.04)
        #t = R.TText(0.55,0.65,"Chi^2 = "+str(round(chi2,2)))
        tlx,tly = pars["range"][0][1]*0.55, pars["range"][1][1]*0.9
        tl = R.TLatex(tlx,tly,"\chi^{2} = "+str(round(chi2,2)))
        leg.AddEntry(tl,"","")
        leg.Draw() 
     
    #ratio.SetSeparationMargin(0.0)
    #ratio.GetUpperRefYaxis().SetTitle(pars["ylabel"])
    #ratio.GetLowerRefXaxis().SetTitle(pars["xlabel"])
    #ratio.GetLowerRefYaxis().SetTitle("MC / Data")
    #ratio.GetLowerRefYaxis().SetRangeUser(-1000,1000)
    #ratio.GetUpperRefYaxis().SetRangeUser(hr.GetMinimum()-0.1*hr.GetMinimum(),hr.GetMaximum()+0.1*hr.GetMaximum())
    #if "ylim" in pars.keys():
    #    ratio.GetUpperRefYaxis().SetRangeUser(pars["ylim"][0],pars["ylim"][1])
    #ratio.GetLowerRefXaxis().SetRangeUser(pars["range"][0][0],pars["range"][0][1])
        c.Update()
        c.SaveAs("plots/"+pars["title"]+".pdf") 
        c.Clear()

    def fs1_channel(self):
        #first calculate the s1 light frac per channel  (should it be max or just average over all events)
        #m.branches.add_branch("s1maxfrac_run",self.branches["s1npe_ch_top_run"].max(axis=1)/self.branches["s1npe_run"])
        '''
        In Maximos convention order should be 
            TL for G0 G1 G2 G4, TR for G0... ML, MR, BL, BR 
            converted to indices that is pixel_id
        Order of Run arrays: 
        Order of MC  arrays: loops from negative to positive x for every y  
        '''
        group_labels      = ["G0","G1","G2","G3"]
        pixel_labels      = ["TL","TR","ML","MR","BL","BR"]
        pixel_id          = [20,23,0,3,21,22,1,2,16,19,4,7,17,18,5,6,12,15,8,11,13,14,9,10]
        s1_ch_frac_maximo = [0.033,0.0332,0.033,0.036,0.0365,0.0395,0.0396,0.041,0.042,0.040,0.0415,0.0425,0.044,0.0445,0.045,0.046,0.041,0.044,0.0420,0.045,0.046,0.047,0.049,0.0495]
        first   = np.divide(self.branches["s1npe_ch_top_mc"],np.sum(self.branches["s1npe_ch_top_mc"],axis=1).reshape(len(self.branches["s1npe_ch_top_mc"]),1) )
        second = [np.average(first[:,i]) for i in range(24)]
        s1_ch_frac_mc  = self.m.helper.reorder(second,pixel_id)
        first   = np.divide(self.branches["s1npe_ch_top_run"],np.sum(self.branches["s1npe_ch_top_run"],axis=1).reshape(len(self.branches["s1npe_ch_top_run"]),1) )
        s1_ch_frac_data = [np.average(first[:,i]) for i in range(24)]
        f,ax = plt.subplots() 
        ax.plot([i for i in range(len(pixel_id))], s1_ch_frac_data, "--.b",label = "Run")
        #ax.plot([i for i in range(len(pixel_id))], s1_ch_frac_maximo, "--.k",label = "Maximo's Light Fraction")
        ax.plot([i for i in range(len(pixel_id))], s1_ch_frac_mc, "--.r",label = "MC")
        for i in range(1,6):
            ax.axvline(x=(i*4+0.5)-1,color = 'black',linestyle = '--')
        for i,label in enumerate(pixel_labels):
            plt.text(4*i+0.5,0.055,label)
        ax.set_xticks([i for i in range(len(pixel_id))])
        ax.set_xticklabels(group_labels*6)
        ax.set_ylabel(r"$f_{s1}$")
        ax.set_ylim(0.0245,0.06)
        plt.legend(loc="lower right")
        plt.savefig("plots/s1chanfrac.pdf",dpi=300)
        plt.clf()
        plt.close()
    
    def fs2_channel(self):
        #not sure about this order
        pixel_id = [i for i in range(24)]
        s2_ch_frac_maximo = [0.05652,0.0365,0.042,0.049,0.0465,0.0375,0.0389,0.046,0.044,0.031,0.0385,0.04,0.03750,0.0322,0.0295,0.0422,0.041,0.04,0.032,0.0442,0.0503,0.047,0.0445,0.0565]
        first   = np.divide(self.branches["s2npe_ch_top_run"],np.sum(self.branches["s2npe_ch_top_run"],axis=1).reshape(len(self.branches["s2npe_ch_top_run"]),1) )
        s2_ch_frac_data  = [np.average(first[:,i]) for i in range(24)]
        first   = np.divide(self.branches["s2npe_ch_top_mc"],np.sum(self.branches["s2npe_ch_top_mc"],axis=1).reshape(len(self.branches["s2npe_ch_top_mc"]),1) )
        second = [np.average(first[:,i]) for i in range(24)]
        s2_ch_frac_mc = self.m.helper.reorder(second,pixel_id)
        f,ax = plt.subplots() 
        ax.plot([i for i in range(24)], s2_ch_frac_maximo, "--.k",label = "Maximo's Light Fraction")
        ax.plot([i for i in range(24)], s2_ch_frac_mc, "--.r",label = "MC Light Fraction")
        ax.plot([i for i in range(24)], s2_ch_frac_data, "--.b",label = "Data Light Fraction")
        ax.set_ylabel(r"$f_{s2}$")
        ax.set_ylim(0.025,0.06)
        plt.legend(loc='best')
        plt.savefig("plots/s2chanfrac.png",dpi=300) 
        plt.clf()
        plt.close()

    def s1chanfrac(self,**pars):
        hr,hmc = R.TH1F("hr","",24,pars["range"][0],pars["range"][1]), R.TH1F("hmc","",24,pars["range"][0],pars["range"][1]) 
       # nevents,24
        s1chanfrac_run = np.zeros((len(self.branches["s1npe_ch_top_run"]),24))
        for i in range(len(self.branches["s1npe_ch_top_run"])):
            s1chanfrac_run[i] = self.branches["s1npe_ch_top_run"][i] / np.sum(self.branches["s1npe_ch_top_run"][i]) 
       #     hr.Fill(np.max(self.branches["s1npe_ch_top_run"][i])/np.sum(self.branches["s1npe_ch_top_run"][i]))
        s1chanfrac_mc = np.zeros((len(self.branches["s1npe_ch_top_mc"]),24))
        for i in range(len(self.branches["s1npe_ch_top_mc"])):
            s1chanfrac_mc[i] = self.branches["s1npe_ch_top_mc"][i] / np.sum(self.branches["s1npe_ch_top_mc"][i]) 
       #     hmc.Fill(np.max(self.branches["s1npe_ch_top_mc"][i])/np.sum(self.branches["s1npe_ch_top_mc"][i]))
        s1chanfrac_run,s1chanfrac_mc = np.average(s1chanfrac_run,axis=0),np.average(s1chanfrac_mc,axis=0)
        plt.plot(np.arange(24),s1chanfrac_run,'-xk',label="Run")
        plt.plot(np.arange(24),s1chanfrac_mc, '-xr',label = "MC")
        plt.xlabel("Channel ID")
        plt.ylabel("S1ChanFrac")
        plt.grid(True)
        plt.ylim(0.01,0.08)
        plt.legend(loc='best')
        plt.savefig("plots/"+pars["title"]+".pdf") 

    def r_legend(self,hmc,hr):
            legend = R.TLegend(0.7,0.6,0.85,0.75)
            legend.AddEntry(hmc,"MC")
            legend.AddEntry(hr,"Data")
            legend.SetLineWidth(0)
            return legend
    # overall max channel/s fraction

    # plotting the channel's max channel fraction

