"""Provides classes for plotting residuals.

ResidualPlot: Plot two-dimensional residuals.
RadialResidualPlot: Plot residuals projected to radial axis.
AngularResidualPlot: Plot residuals projected to angular axis.
"""

from math import atan, pi

from ROOT import TF1, TLatex, TPad, TProfile
from ROOT.TH1 import kPoisson

from lib.plot.plot import ColorBase, PlotBase, SingleHistBase

class ResidualPlot(ColorBase):
    """Plot two-dimensional residuals.

    __init__: Initialize.
    """

    def __init__(self, hist, fill=4954, workinprogress=False):
        """Initialize a two-dimensional residual plot."""
        ColorBase.__init__(self, hist, 'res', fill, workinprogress)
        self._xtitle = 'x [cm]'
        self._ytitle = 'y [cm]'
        self._ztitle = 'Pulls'
        self.zrange(-5.0, 5.0)
        for xbin in range(self.xaxis().GetNbins()+2):
            for ybin in range(self.yaxis().GetNbins()+2):
                if self.get_bin(xbin, ybin) < -4.9999:
                    self.set_bin(-4.9999, xbin, ybin)
                if self.get_bin(xbin, ybin) == 0.0:
                    self.set_bin(-10.0, xbin, ybin)

class RadialResidualPlot(SingleHistBase):
    """Plot residuals projected to one radial dimension.

    __init__: Initialize.
    chisq: Return the chi-squared summed over radial bins.
    """

    def __init__(
        self, data, model, nbins=None, maxr=None, calcr=None, skipzero=False,
        fill=4954, workinprogress=False
    ):
        """Project residual histogram to one-dimensional plot."""
        res = self.create_residual(
            data, model, nbins=nbins, maxr=maxr, calcr=calcr, skipzero=skipzero
        )
        SingleHistBase.__init__(self, res, 'radialRes', fill, workinprogress)
        self._xtitle = 'r [cm]'
        self._ytitle = 'Pull'

    @staticmethod
    def create_residual(
        data, model, nbins=None, maxr=None, calcr=None, skipzero=False
    ):
        """Create radial residual histogram."""
        if nbins is None:
            nbins = data.GetXaxis().GetNbins()
        if maxr is None:
            maxr = data.GetXaxis().GetXmax()
        if calcr is None:
            calcr = lambda x, y: (x**2+y**2)**0.5
        data.SetBinErrorOption(kPoisson)
        res = TProfile('radialRes', '', nbins, 0.0, maxr)
        for xbin in range(data.GetXaxis().GetNbins()):
            for ybin in range(data.GetYaxis().GetNbins()):
                d = data.GetBinContent(xbin+1, ybin+1)
                if skipzero and d == 0.0:
                    continue
                m = model.GetBinContent(xbin+1, ybin+1)
                if m < d:
                    e = data.GetBinErrorLow(xbin+1, ybin+1)
                else:
                    e = data.GetBinErrorUp(xbin+1, ybin+1)

                x = data.GetXaxis().GetBinCenter(xbin+1)
                y = data.GetYaxis().GetBinCenter(ybin+1)

                if e!=0.:
                    res.Fill(calcr(x, y), (d-m)/e)
        return res

    def chisq(self):
        """Return the chi-squared summed over the radial bins."""
        result = 0.0
        for i in range(self.xaxis().GetNbins()):
            result += self.get_bin(i+1) ** 2
        return result

class AngularResidualPlot(SingleHistBase):
    """Plot residuals projected to one angular dimension.

    __init__: Initialize.
    chisq: Return the chi-squared summed over radial bins.
    """

    def __init__(
        self, data, model, nbins=None, calcphi=None, skipzero=False, fill=4954,
        workinprogress=False
    ):
        """Project residual histogram to one-dimensional plot."""
        res = self.create_residual(
            data, model, nbins=nbins, calcphi=calcphi, skipzero=skipzero
        )
        SingleHistBase.__init__(self, res, 'angularRes', fill, workinprogress)
        self._xtitle = '#phi [rad]'
        self._ytitle = 'Pull'

    @staticmethod
    def create_residual(data, model, nbins=None, calcphi=None, skipzero=False):
        """Create angular residual histogram."""
        if nbins is None:
            nbins = data.GetXaxis().GetNbins()
        if calcphi is None:
            calcphi = lambda x, y: atan(x/y)
        data.SetBinErrorOption(kPoisson)
        res = TProfile('angularRes', '', nbins, -0.5*pi, 0.5*pi)
        for xbin in range(data.GetXaxis().GetNbins()):
            for ybin in range(data.GetYaxis().GetNbins()):
                d = data.GetBinContent(xbin+1, ybin+1)
                if skipzero and d == 0.0:
                    continue
                m = model.GetBinContent(xbin+1, ybin+1)
                if m < d:
                    e = data.GetBinErrorLow(xbin+1, ybin+1)
                else:
                    e = data.GetBinErrorUp(xbin+1, ybin+1)
                x = data.GetXaxis().GetBinCenter(xbin+1)
                y = data.GetYaxis().GetBinCenter(ybin+1)
                if e!=0.:
                    res.Fill(calcphi(x, y), (d-m)/e)

        return res

    def chisq(self):
        """Return the chi-squared summed over the radial bins."""
        result = 0.0
        for i in range(self.xaxis().GetNbins()):
            result += self.get_bin(i+1) ** 2
        return result

class CombinedResidualPlot(PlotBase):
    """"Plot residuals projected to radial and angular dimension.

    __init__: Initialize.
    draw: Draw the plots.
    """

    def __init__(
        self, data, model, nbins=None, maxr=None, calcr=None, calcphi=None,
        skipzero=False, fill=4954, workinprogress=False
    ):
        """Project residual histogram to one-dimensional plots."""
        PlotBase.__init__(
            self, name='combinedRes', fill=fill, workinprogress=workinprogress
        )
        self._rad = RadialResidualPlot.create_residual(
            data, model, nbins=nbins, maxr=maxr, calcr=calcr, skipzero=skipzero
        )
        self._ang = AngularResidualPlot.create_residual(
            data, model, nbins=nbins, calcphi=calcphi, skipzero=skipzero
        )
        self._rtitle = 'r [cm]'
        self._y1title = 'Pull'
        self._ptitle = '#phi [rad]'
        self._y2title = 'Pull'
        if maxr is None:
            self._rrange = None
        else:
            self._rrange = (0, maxr)
        self._prange = (-pi/2, pi/2)
        self._y1range = None
        self._y2range = None

    def SetName(self, name):
        PlotBase.SetName(self, name)
        self._rad.SetName('{0}_rad'.format(name))
        self._ang.SetName('{0}_ang'.format(name))

    def draw(self, canvas=None):
        if canvas is None:
            canvas = self
        canvas.Divide(1, 2)
        canvas.cd(1).SetMargin(0.1, 0.1, 0.05, 0.2)
        self._rad.Draw()
        canvas.cd(2).SetMargin(0.1, 0.1, 0.2, 0.05)
        self._ang.Draw()
        canvas.cd()
        canvas.Update()
        for c, axis in [
            ('r', self._rad.GetXaxis()), ('y1', self._rad.GetYaxis()),
            ('p', self._ang.GetXaxis()), ('y2', self._ang.GetYaxis())
        ]:
            axis.SetLabelSize(0.05)
            axis.SetTitleSize(0.07)
            axis.SetLabelOffset(0.01)
            rng = getattr(self, '_{0}range'.format(c))
            if rng is not None:
                axis.SetRangeUser(rng[0], rng[1])
            ttl = getattr(self, '_{0}title'.format(c))
            if c.startswith('y') and ttl is not None:
                axis.SetTitle(ttl)
            if c.startswith('y'):
                axis.SetTitleOffset(0.6)
                axis.SetNdivisions(505)
        canvas.cd()
        self._pad = TPad('newpad', '', 0, 0, 1, 1)
        self._pad.SetFillStyle(4000)
        self._pad.Draw()
        self._pad.cd()
        if self._drawCMS:
            self.draw_cms()
        text = TLatex()
        text.SetNDC()
        text.SetTextFont(42)
        text.SetTextSize(0.035)
        text.SetTextAlign(13)
        text.SetTextAngle(90)
        if self._rtitle is not None:
            text.DrawLatexNDC(0.905, 0.54, self._rtitle)
        if self._ptitle is not None:
            text.DrawLatexNDC(0.905, 0.115, self._ptitle)
        for obj in self._container_draw:
            obj.Draw('SAME')
        for c in (canvas.cd(1), canvas.cd(2), canvas):
            c.Modified()
            c.Update()
        canvas.cd()
