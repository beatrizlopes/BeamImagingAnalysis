
import ROOT
from lib.plot.plot import ColorBase
ColorBase.kBird()

Fill=4266
BCID=51
Beam=1
Direction="Y"

tfile = ROOT.TFile.Open('/afs/cern.ch/user/b/bribeiro/CMSSW_9_4_11/src/BeamImagingAnalysis/beamimaging/Fill{}_rereco_many_medium.root'.format(Fill))
ttree = tfile.Get('Beam{0}Move{1}_bunch{2}Add'.format(Beam,Direction,BCID))
text = ROOT.TLatex()
text.SetNDC()
text.SetTextFont(62)
text.SetTextSize(0.0375)
text.SetTextAlign(13)
for i in range(19):
    hist = ROOT.TH2F('hist{}'.format(i), '', 50, 0.0, 0.12, 50, 0.04, 0.16)
    hist.GetZaxis().SetRangeUser(1.0, 400.0)
    ttree.Draw('vtx_y:vtx_x>>hist{}'.format(i), 'vtx_nTrk>14 && scanstep=={}'.format(i+1), 'goff')
    plot = ColorBase(hist)
    plot._xtitle = 'x [cm]'
    plot._ytitle = 'y [cm]'
    plot.draw()
    plot.SetLogz()
    text.DrawLatex(0.15, 0.88, 'Scan {0}{1}, step {2}'.format(Beam,Direction,i))
    plot.Update()
    plot.SaveAs('Fill{0}_rereco_{1}{2}_bcid{3}_step{4}.png'.format(Fill,Beam,Direction,BCID,i))
    plot.Close()
print 'convert -delay 20 -loop 0 Fill{0}_rereco_{1}{2}_bcid{3}_step{0..18}.png Fill{0}_rereco_{1}{2}_bcid{3}_allsteps.gif'.format(Fill,Beam,Direction,BCID)
