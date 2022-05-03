import numpy as np
from manager import Manager
import warnings
warnings.filterwarnings("ignore")

m = Manager()

class FineTune():
    def __init__(self,fmc,fReD=None):  
        self.fmc      = fmc
        self.fred     = fReD
        #Read both the MC and the Run File 
        #s1min ,s1max  = 400,600
        s1min ,s1max = 340,460
        #Ws1min ,s1max = 2,10000
        self.branches = m.branches(self.fmc,"mc") 
        #m.branches.cut_branches((self.branches["s1pe_mc"]>0),data="mc")
        m.branches.cut_branches((self.branches["s1npe_mc"]>s1min) & (self.branches["s1npe_mc"]<s1max),data="mc")
        m.branches.cut_branches((self.branches["z_mc"]<2.5) & (self.branches["z_mc"]>-2.5),data="mc")
        
        m.plotter.get_branches(self.branches)

        #look at MC
        m.plotter.oneDplotter(d="s1ene_mc",               d_mc=False,bins=200,range=(0,100),                 xlabel="Ene [KeV]",  ylabel="Events",               title="s1ene_mc")
        m.plotter.oneDplotter(d="tpcene_mc",              d_mc=False,bins=200,range=(0,100),                 xlabel="Ene [KeV]",  ylabel="Events",               title="s1tpcene_mc")
        m.plotter.oneDplotter(d="s1npe_mc",               d_mc=False,bins=70,range=(300,500),               xlabel="S1 [PE]",    ylabel="Events",               title="s1_mc",fit='gaus')
        m.plotter.oneDplotter(d="s1lyield_mc",            d_mc=False,bins=100,range=(0,35),               xlabel="LYield [PE/keV]",    ylabel="Events",               title="s1lyield_mc")
        m.plotter.twoDplotter(d="s1ene_mc,s1npe_mc",      d_mc=False,bins=100,range=((0,300),(0,5000)),     xlabel="Energy [keV]",ylabel="S1 [NPE]",            title="enes1_mc")
        m.plotter.twoDplotter(d="s1npe_mc,s1maxfrac_mc",  d_mc=False,bins=100,range=((0,5000),(0.04,0.15)),     xlabel="S1 [NPE]",   ylabel="S1 Max Frac",          title="s1s1maxfrac_mc")
        m.plotter.twoDplotter(d="z_mc,s1maxfrac_mc",      d_mc=False,bins=100,range=((-2.5,2.5),(0.04,0.15)),   xlabel="Z [cm]",     ylabel="S1 Max Frac",          title="zs1maxfrac_mc")
        m.plotter.twoDplotter(d="z_mc,s1npe_mc",          d_mc=False,bins=100,range=((-2.5,2.5),(0,10000)), xlabel="Z [cm]",     ylabel="S1 [PE]",              title="zs1_mc")
        m.plotter.twoDplotter(d="z_mc,s1tba_mc",          d_mc=False,bins=100,range=((-2.5,2.5),(-0.4,0.4)),    xlabel="Z [cm]",     ylabel="S1 TBA",               title="zs1tba_mc")
        m.plotter.twoDplotter(d="s1tba_mc,s1npe_mc",      d_mc=False,bins=100,range=((-1,1),(0,10000)),     xlabel="S1 TBA",     ylabel="S1 [PE]",              title="s1tbas1_mc")
        m.plotter.twoDplotter(d="z_mc,s1lyield_mc",       d_mc=False,bins=100,range=((-2.5,2.5),(0,35)),    xlabel="Z [cm]",     ylabel="Light Yield [PE/keV]", title="zlyield_mc")

        if self.fred != None:
            self.branches     = m.branches(self.fred,"run")
            s1min ,s1max  = 450, 600
            #s1min ,s1max  = 1,10000
            #s1min, s1max = 0,500000
            m.branches.cut_branches(self.branches["s1maxfrac_run"]>0.001)
            m.branches.cut_branches((self.branches["z_run"]<2.5) & (self.branches["z_run"]>-2.5))
            m.branches.cut_branches((self.branches["s1npe_run"]>s1min) & (self.branches["s1npe_run"]<s1max))
            if "S2" in m.config("run","mc","str"):
                m.branches.cut_branches(self.branches["s2npe_run"]>2000)
                m.branches.cut_branches(self.branches["s2npe_mc"]>2000,data="mc")

            m.plotter.get_branches(self.branches)

            #look at data
            #m.plotter.lyield_z(d="s1npe_run,z_run" ,                d_mc =True,bins=50,range=(0,2),                   xlabel="Z / Z_{1/2}",     ylabel="Y / Y(z_{1/2})",       title="yield", fmin=420,fmax=570)       
            #m.plotter.lyield_z(d="s1npe_run,z_run" ,                d_mc =True,bins=50,range=(0,2),                   xlabel="Z / Z_{1/2}",     ylabel="Y / Y(z_{1/2})",       title="yield", fmin=420,fmax=570)       
            m.plotter.oneDplotter(d="s1npe_run",    d_res=True,     d_mc=True, bins=200,range=(0,4000),            xlabel="S1 [PE]",    ylabel="Events",title="s1")
            m.plotter.oneDplotter(d="s2npe_run",    d_res=True,     d_mc=False, bins=100,range=(0,20000),           xlabel="S2 [PE]",    ylabel="Events",title="s2")
            m.plotter.oneDplotter(d="tdrift_run",                   d_mc=False,bins=100,range=(0,100),             xlabel="tdrift[us]", ylabel="Events",       title="tdrift")
            m.plotter.twoDplotter(d="s1npe_run,s2npe_run",          d_mc=False,bins=100,range=((0,1000),(0,20000)),xlabel="S1 [PE]",    ylabel="S2 [Pe]",      title="s1s2") 
            m.plotter.twoDplotter(d="s1npe_run,s2npe_run",          d_mc=False,bins=100,range=((0,1000),(0,20000)),xlabel="S1 [PE]",    ylabel="S2 [Pe]",      title="s1s2")
            m.plotter.twoDplotter(d="tdrift_run,s2npe_run",         d_mc=False,bins=100,range=((0,60),  (0,20000)),xlabel="TDrift [us]",ylabel="S2 [Pe]",      title="tdrifts2")
            m.plotter.twoDplotter(d="s1npe_run,cl_qfp_run",         d_mc=False,bins=100,range=((0,5000), (0,1)),   xlabel="S1 [PE]",    ylabel="f400",         title="f400")
            m.plotter.oneDplotter(d="cl_npeaks_run",                d_mc=False,bins=100,range=(0,200),             xlabel="S1 NPeaks",  ylabel="Events",       title="s1npeaks")
            m.plotter.twoDplotter(d="s1maxfrac_run,cl_npeaks_run",  d_mc=False,bins=100,range=((0,0.05),(0,200)),  xlabel="S1 Max Frac",ylabel="S1 Peaks",     title="s1maxfracpeaks")
            m.plotter.twoDplotter(d="tdrift_run,s1npe_run",     d_run=False,bins=100,range=((0,60),(0,1000)),     xlabel="TDrift [us]",ylabel="S1 [NPE]",            title="tdrifts1_run")
            m.plotter.twoDplotter(d="tdrift_run,s1tba_run",         d_mc=False,bins=100,range=((0,60),(-0.4,0.4)), xlabel="TDrift [us]",ylabel="S1 TBA",       title="tdrifts1tba")
            m.plotter.twoDplotter(d="tdrift_run,s1maxfrac_run",     d_mc=False,bins=100,range=((0,60),(0,0.1)),    xlabel="TDrift [us]",ylabel="S1 Max Frac",  title="tdrifts1maxfrac")
            m.plotter.twoDplotter(d="z_run,s1maxfrac_run",          d_mc=False,bins=100,range=((-2.5,2.5),(0,0.1)),xlabel="Z [cm]",     ylabel="S1 Max Frac",  title="zs1maxfrac_run")
            m.plotter.twoDplotter(d="s1npe_run,s1maxfrac_run",      d_mc=False,bins=100,range=((0,5000),(0,0.1)),  xlabel="S1 [NPE]",   ylabel="S1 Max Frac",  title="s1s1maxfrac")
            m.plotter.twoDplotter(d="tdrift_run,cl_npeaks_run",     d_mc=False,bins=100,range=((0,60),(0,200)),    xlabel="TDrift [us]",ylabel="S1 NPeaks",    title="tdrifts1peaks")
            m.plotter.twoDplotter(d="z_run,s1tba_run",              d_mc=False,bins=100,range=((-3,3),(-0.5,0.5)), xlabel="Z [cm]",          ylabel="S1 TBA",               title="zs1tba_run")
            m.plotter.twoDplotter(d="z_run,s1npe_run",              d_mc=False,bins=100,range=((-2.5,2.5),(0,10000)),  xlabel="Z [cm]",          ylabel="S1 [PE]",              title="zs1_run")
            m.plotter.twoDplotter(d="s1tba_run,s1npe_run",          d_mc=False,bins=100,range=((-1,1),(0,10000)),      xlabel="S1 TBA",          ylabel="S1 [PE]",              title="s1tbas1_run")

       #    # compare to mc
            m.plotter.oneDplotter(d="s1topbot_run",                 d_mc=True,grid=True, bins=100,range=(0,2),                   xlabel="S1 Top / Bottom", ylabel="Events",        title="s1topbot")
            m.plotter.oneDplotter(d="s2topbot_run",                 d_mc=True,grid=True, bins=50, range=(0.5,5),                 xlabel="S2 Top / Bottom", ylabel="Events",        title="s2topbot")
            m.plotter.oneDplotter(d="s1tba_run",                    d_mc=True,grid=True, bins=100,range=(-1,1),                  xlabel="S1 TBA",          ylabel="Events",        title="s1tba")
            m.plotter.oneDplotter(d="s2tba_run",                    d_mc=True,grid=True, bins=50, range=(-0.3,0.7),              xlabel="S2 TBA",          ylabel="Events",        title="s2tba")
            m.plotter.oneDplotter(d="s1maxfrac_run",                d_mc=True,grid=True, bins=50, range=(0.04,0.15),                 xlabel="S1 Max Fraction", ylabel="Events",        title="s1maxfrac")
            m.plotter.oneDplotter(d="s2maxfrac_run",                d_mc=True,grid=True, bins=50, range=(0,0.5),                 xlabel="S2 Max Fraction", ylabel="Events",        title="s2maxfrac")        
            m.plotter.profiler(d="z_run,s1maxfrac_run" ,            d_mc=True,grid=True, bins=50, range=((-2.5,2.5),(0.04,0.15)),       xlabel="Z [cm]",     ylabel="S1 Max Frac",    title="zs1maxfrac",ylim = (0.04,0.015))
            m.plotter.profiler(d="z_run,s1tba_run",                 d_mc=True,grid=True, bins=50,range=((-2.5,2.5),(-0.3,0.3)),      xlabel="Z [cm]",     ylabel="S1 TBA",               title="zs1tba",ylim=(-0.3,0.3))
            m.plotter.profiler(d="z_run,s2tba_run",                 d_mc=True,grid=True, bins=50,range=((-2.5,2.5),(-0.1,1)),        xlabel="Z [cm]",     ylabel="S2 TBA",               title="zs2tba",ylim=(-0.1,1))

            #m.plotter.s1chanfrac(range = (0.04,0.15),xlabel="S1ChanFrac",ylabel="S1ChanFrac",title = "s1chanfrac")
            #m.plotter.fs1_channel()
            #m.plotter.fs2_channel()

##### Run ####
fmc  = m.config('run','mc','str')
fred = m.config('run','path','str')+m.config('run','run','str')
tune = FineTune(fmc,fred)