import sys
import uproot as ur
import numpy as np

class Branches():
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)

    def __init__(self,manager):
        self.m              = manager
        self.branches_run   = manager.config("branches","run","str").split(",")
        self.branches_mc    = manager.config("branches","mc","str").split(",")
        self.branches       = {}
        self.nevents        = manager.config("run","nevents","int")
        self.channel_map    = manager.config("channels","channel_na","str").split(",")

    def __repr__(self) -> str:
        return self.branches.keys()

    def __call__(self,file,mode):    
        self.mode = mode 
        branches_names = self.branches_mc if mode == "mc" else self.branches_run
        for branch in branches_names:
            print("Reading file, branch",file,branch)
            self.branches[branch] = self.read_branch(file,branch) 

        if self.mode == "mc":
            self.change_branch(self.branches_mc,[branch+"_mc" for branch in self.branches_mc])
            # add up the bottom channels to match the 4 channel configration in ReD
            self.branches["s1npe_ch_bot_mc"]  = self.sum_npes_bot(self.branches["s1npe_ch_bot_mc"])
            self.branches["s2npe_ch_bot_mc"]  = self.sum_npes_bot(self.branches["s2npe_ch_bot_mc"])
            self.add_branch("s1lyield_"+self.mode,self.branches["s1npe_mc"]/self.branches["tpcene_mc"])

        if self.mode == "run":
            self.change_branch(self.branches_run,[branch+"_run" for branch in self.branches_run])
            s1,s2, tdrift, cut_ncl  = self.get_s1ns2(self.branches["cl_q_run"],self.branches["ncl_s1_run"],self.branches["ncl_s2_run"],self.branches["cl_startt_run"],self.branches["cl_endt_run"])
            self.cut_branches(cut_ncl)
            self.add_branch("s1npe_run", s1)
            self.add_branch("s2npe_run", s2)
            self.add_branch("tdrift_run",tdrift)
            self.branches["cl_startt_run"] = self.turn_array(self.branches["cl_startt_run"])
            self.branches["cl_qfp_run"]    = self.turn_array(self.branches["cl_qfp_run"])
            self.branches["cl_npeaks_run"] = self.turn_array(self.branches["cl_npeaks_run"])

            pe_ch_top, pe_ch_bot = self.get_pe_ch(self.branches["pk_t_run"],self.branches["pk_npe_run"],self.branches["pk_ch_run"],self.branches["cl_startt_run"],self.branches["cl_endt_run"])
            self.add_branch("s1npe_ch_top_run",pe_ch_top[:,:,0])
            self.add_branch("s1npe_ch_bot_run",pe_ch_bot[:,:,0])
            self.add_branch("s2npe_ch_top_run",pe_ch_top[:,:,1])
            self.add_branch("s2npe_ch_bot_run",pe_ch_bot[:,:,1])
            z = self.get_zpos(self.branches["tdrift_run"])
            self.add_branch("z_run",z)

        #find the tba branches
        for s in ["s1","s2"]:
            self.branches[s+"npe_ch_top_"+self.mode][self.branches[s+"npe_ch_top_"+self.mode].sum(axis=1)!=0]
            self.branches[s+"npe_ch_bot_"+self.mode][self.branches[s+"npe_ch_bot_"+self.mode].sum(axis=1)!=0]

        stba = self.get_stba(self.mode)
        self.add_branch("s1tba_"+self.mode,stba[:,0])
        self.add_branch("s2tba_"+self.mode,stba[:,1])       
        self.cut_branches(self.branches["s1npe_ch_bot_"+self.mode].sum(axis=1)>0, data=self.mode)
        self.cut_branches(self.branches["s2npe_ch_bot_"+self.mode].sum(axis=1)>0, data=self.mode)
        self.cut_branches(self.branches["s1npe_ch_top_"+self.mode].sum(axis=1)>0, data=self.mode)
        self.cut_branches(self.branches["s2npe_ch_top_"+self.mode].sum(axis=1)>0, data=self.mode)        
        self.add_branch("s1topbot_"+self.mode,self.branches["s1npe_ch_top_"+self.mode].sum(axis=1)/self.branches["s1npe_ch_bot_"+self.mode].sum(axis=1))
        self.add_branch("s2topbot_"+self.mode,self.branches["s2npe_ch_top_"+self.mode].sum(axis=1)/self.branches["s2npe_ch_bot_"+self.mode].sum(axis=1))
        self.add_branch("s1maxfrac_"+self.mode,self.branches["s1npe_ch_top_"+self.mode].max(axis=1)/np.sum((self.branches["s1npe_ch_top_"+self.mode]),axis=1))
        self.add_branch("s2maxfrac_"+self.mode,self.branches["s2npe_ch_top_"+self.mode].max(axis=1)/np.sum((self.branches["s2npe_ch_top_"+self.mode]),axis=1))
        
        return self.branches

    def add_branch(self,branch,data):
        self.branches[branch] = data

    def read_branch(self,fname,branch):       
        return ur.open(fname)["dstree"][branch].array(library="np")[:self.nevents]

    def cut_branches(self,cut,data="run"):
        for key in self.branches.keys():
            if data in key:
                self.branches[key] = self.branches[key][cut]

    def change_branch(self,olds, news):
        for i in range(len(olds)):
            self.branches[news[i]] = self.branches[olds[i]]
            del self.branches[olds[i]]

    def get_s1ns2(self,cl_q,ncl_s1,ncl_s2,cl_startt,cl_endt):
        print("Getting S1 & S2 & tdrift")
        cut_ncl = (ncl_s2 == 1) & (ncl_s1 == 1)
        cl_q = np.stack(cl_q[cut_ncl],axis=1).T
        cl_startt = np.stack(cl_startt[cut_ncl],axis=1).T
        cl_endt   = np.stack(cl_endt[cut_ncl],axis=1).T
        tdrift = (cl_startt[:,1]-cl_startt[:,0]) *2e-3

        return cl_q[:,0], cl_q[:,1], tdrift, cut_ncl

    def turn_array(self,branch): return np.stack(branch,axis=1).T

    def get_zpos(self,tdrift):    
        field    = 200 
        charge   = 1.6*10**(-19)
        m        = 9.31*10**(-31)
        speed_e  = ((2*charge*field)/m)**0.5
        return -((speed_e * tdrift * 1e-9 * 10) - 2.5)                # convert to s, cm

    def get_stba(self, mode):    
        if mode == "run":
            keys = ["s1npe_ch_top_run","s1npe_ch_bot_run","s2npe_ch_top_run","s2npe_ch_bot_run"]
        else:
            keys = ["s1npe_ch_top_mc","s1npe_ch_bot_mc","s2npe_ch_top_mc","s2npe_ch_bot_mc"] 
        stba = np.zeros((len(self.branches[keys[0]]),2))
        stba[:,0] = (np.sum(self.branches[keys[0]],axis=1) - np.sum(self.branches[keys[1]],axis=1) ) / (np.sum(self.branches[keys[0]],axis=1) + np.sum(self.branches[keys[1]],axis=1) ) 
        stba[:,1] = (np.sum(self.branches[keys[2]],axis=1) - np.sum(self.branches[keys[3]],axis=1) ) / (np.sum(self.branches[keys[2]],axis=1) + np.sum(self.branches[keys[3]],axis=1) )

        return stba

    def get_pe_ch(self,pk_t,pk_npe,pk_ch,cl_startt,cl_endt):
        print("Getting PEs per channel, for top bottom, S1 & S2")
        pe_ch_top = np.zeros((len(pk_t),24,2),dtype = np.float32)
        pe_ch_bot = np.zeros((len(pk_t),4,2),dtype = np.float32)
        for event in range(len(pk_t)):
            if event % 10000==0:
                print(event,"/",len(pk_t))
            for s in range(2):                      # do s1 and s2
                try:
                    t0 = cl_startt[event][s]        # this is 0 for s1, 1 for s2
                    t1 = cl_endt[event][s]
                    mask = (pk_t[event] > t0) & (pk_t[event] < t1)   
                except:
                    ValueError("Value too large")
                    continue 
                pnpe = np.array(pk_npe[event][mask])  # this is the same as cl_p for s2
                ch = np.array(np.unique(pk_ch[event][mask]))
                if len(ch) == 0:
                    continue
                ch_w = np.zeros(ch.shape,dtype=np.float)
                for i,c in np.ndenumerate(ch):
                    ch_w[i] = sum(pnpe[pk_ch[event][mask]==c]) 
                if ch_w.shape[0] == 28:
                    pe_ch_top[event,:,s], pe_ch_bot[event,:,s] = self.m.helper.map_channel(ch_w,self.channel_map)
        
        return pe_ch_top,pe_ch_bot
    
    def sum_npes_bot(self,bot):
        new_bot = np.zeros((len(bot),4))        #number of events, 4 channles
        for i in range(len(bot)):
            new_bot[i,0] = np.sum(bot[i,list(map(int,self.m.config("channels","bottom_q1","str").split(",")))])
            new_bot[i,1] = np.sum(bot[i,list(map(int,self.m.config("channels","bottom_q2","str").split(",")))]) 
            new_bot[i,2] = np.sum(bot[i,list(map(int,self.m.config("channels","bottom_q3","str").split(",")))])
            new_bot[i,3] = np.sum(bot[i,list(map(int,self.m.config("channels","bottom_q4","str").split(",")))])
        
        return new_bot
    
