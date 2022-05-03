import numpy as np
import ROOT as R 
import uproot as ur

#files = ["Kr83m_ReD_S1_ESRAcrylic_p8.root","Kr83m_ReD_S1_ESRAcrylic_p9.root"]
#labels = ["0.8","0.9"]

#files = ["Kr83m_ReD_S1_UVAbs_p5.root","Kr83m_ReD_S1_UVAbs_1.root","Kr83m_ReD_S1_UVAbs_1p5.root","Kr83m_ReD_S1_UVAbs_2.root"]
#labels = ["0.5m","1m","1.5m","2m"]

#files = ["Kr83m_ReD_S1_ESRAlpha_p05.root","Kr83m_ReD_S1_ESRAlpha_p15.root"]
#labels = ["0.05 rad","0.15 rad"]

#files = ["Kr83m_ReD_S1_TPBAlpha_p05.root","Kr83m_ReD_S1_TPBAlpha_p15.root"]
#labels = ["0.05 rad","0.15 rad"]

# This one no impact
#files = ["Kr83m_ReD_S1_GridUV_p3.root","Kr83m_ReD_S1_GridUV_p5.root","Kr83m_ReD_S1_GridVis_p5.root"]
#labels = ["0.3","0.5","0.5 Vis"]

#files = ["Kr83m_ReD_S1_ArTPBVisTran_0p4.root","Kr83m_ReD_S1_ArTPBVisTran_0p65.root","Kr83m_ReD_S1_ArTPBVisTran_0p7.root","Kr83m_ReD_S1_ArTPBVisTran_0p9.root","Kr83m_ReD_S1_lessbot.root"]
#labels = ["0.4","0.65","0.7","0.9","less"]

files  = ["Kr83m_ReD_S1_best.root","Kr83m_ReD_S1.root"]
labels = ["Old","New"]

#files = ["Kr83m_ReD_S1.root","Kr83m_ReD_S1_old.root"]
#labels = ["After Fix","Before Fix"]

def read_file(file,branch):
    return ur.open(file)["dstree"][branch].array(library="np")

branches = ["s1npe_ch_top","s1npe","s1npe_ch_bot"]

file_branches = {}
for n,file in enumerate(files):
    for branch in branches:
        file_branches[branch+"_"+str(n)] = read_file(file,branch)

c = R.TCanvas()


h = {}
leg = R.TLegend(0.65,0.75,0.75,0.85)
for n,file in enumerate(files):
    h[str(n)] = R.TH1F("h","",100,0.005,0.1)
    for i in range(len(file_branches["s1npe"+"_"+str(n)])):
        if file_branches["s1npe"+"_"+str(n)][i] >0:
            s1maxfrac = np.max(file_branches["s1npe_ch_top"+"_"+str(n)][i])/file_branches["s1npe"+"_"+str(n)][i]
            h[str(n)].Fill(s1maxfrac)
    h[str(n)].SetLineColor(1+n)
    h[str(n)].SetLineWidth(2)
    h[str(n)].Scale(1/h[str(n)].Integral())
    h[str(n)].GetXaxis().SetTitle("S1 Max Frac")
    h[str(n)].GetXaxis().SetTitle("")
    h[str(n)].GetYaxis().SetRangeUser(0,0.07)
    h[str(n)].Draw("hist" if n==0 else "hist same")
    leg.AddEntry(h[str(n)],labels[n],"l")

leg.Draw()
c.Update()
c.SaveAs("plots/s1maxfrac_mc_stack.pdf")

h = {}
leg = R.TLegend(0.65,0.75,0.75,0.85)
for n,file in enumerate(files):
    h[str(n)] = R.TH1F("h","",100,-1,1)
    for i in range(len(file_branches["s1npe"+"_"+str(n)])):
        if file_branches["s1npe"+"_"+str(n)][i] >0:
            s1tba = (np.sum(file_branches["s1npe_ch_top"+"_"+str(n)][i]) - np.sum(file_branches["s1npe_ch_bot"+"_"+str(n)][i])) / (np.sum(file_branches["s1npe_ch_bot"+"_"+str(n)][i]) + np.sum(file_branches["s1npe_ch_top"+"_"+str(n)][i]) )
            h[str(n)].Fill(s1tba)
    h[str(n)].SetLineColor(1+n)
    h[str(n)].SetLineWidth(2)
    h[str(n)].Scale(1/h[str(n)].Integral())
    h[str(n)].GetXaxis().SetTitle("S1TBA")
    h[str(n)].GetXaxis().SetTitle("")
    h[str(n)].GetYaxis().SetRangeUser(0,0.07)
    h[str(n)].Draw("hist" if n==0 else "hist same")
    leg.AddEntry(h[str(n)],labels[n],"l")

leg.Draw()
c.Update()
c.SaveAs("plots/s1tba_mc_stack.pdf")
