"""Provides function for preparation of Beam Imaging data.

make_filelist: Scan raw files for needed data.
make_trees: Create Beam Imaging trees from raw data.
make_histograms_forbcid: Create Beam Imaging histograms from created trees for each bcid.
make_histograms: Create Beam Imaging histograms from created trees summing over bcid.
"""
import re
import pprint
import json
from os import listdir, stat
from os.path import exists

from ROOT import gDirectory, TFile, TH1F, TH2F, TLatex, gROOT, kTRUE
from lib.plot.plot import ColorBase
ColorBase.kBird()

gROOT.SetBatch(kTRUE)

from lib.io import RootChain, RootTree

def make_filelist(
    directories, times, treename='lumi/tree', timestamp='timeStamp_begin',
    verbose=True
):
    """Scan raw files for data belonging to Beam Imaging scans.

    directories: List of directory paths (string) to be searched.
    times: Dictionary with scan names and 2-tuples (begin & end timestamp).
    treename: Name of the tree that contains the data.
    timestamp: Name of the tree's branch that contains the time.
    returns dictionary with scan names and filelists (list of strings).
    """
    filelist = {scan: [] for scan in times}
    for directory in directories:
        if verbose:
            print '<<< Enter directory: {0}'.format(directory)
        files = [
            f for f in listdir(directory) if f.endswith('.root')
            if stat('{0}/{1}'.format(directory, f)).st_size > 0.0
        ]
        for filename in files:
            tfile = TFile('{0}/{1}'.format(directory, filename))
            if not tfile:
                continue
            tree = tfile.Get(treename)
            if not tree:
                continue
            for scan, (begin, end) in times.iteritems():
                n = tree.GetEntries(
                    '{0}>={1} && {0}<={2}'.format(timestamp, begin, end)
                )
                if n > 0:
                    filelist[scan].append(
                        '{0}/{1}'.format(directory, filename)
                    )
    return filelist

def make_trees(
    filelist, times, bcids, mintrk=0, verbose=True, noerror=False,
    qualities=None
):
    """Run over raw data and collect Beam Imaging data in trees.

    filelist: List of files (string) to be searched.
    times: List of 2-tuples (begin & end timestamp of scan steps).
    bcids: List of bunch crossings (int).
    mintrk: Minimal number of tracks (int) for the event selection (default: 0).
    verbose: Set true to print progress to stdout.
    noerror: Set true to not compute resolutions.
    qualities: tuple of boolean fields that are required.
    returns list of ROOT trees.
    """
    chain = RootChain('lumi/tree')
    chain.add_files(filelist)
    chain.add_fields([
        ('nVtx', 'i', 1),
        ('vtx_nTrk', 'i', 200),
        ('vtx_x', 'f', 200),
        ('vtx_y', 'f', 200),
        ('timeStamp_begin', 'I', 1)
    ])
    if 0 not in bcids:
        chain.add_fields([('bunchCrossing', 'i', 1),])
    if not noerror:
        chain.add_fields([
            ('vtx_xError', 'f', 200),
            ('vtx_yError', 'f', 200),
        ])
    if qualities is None:
        qualities = ('vtx_isGood', '!vtx_isFake')
    for quality in qualities:
        chain.add_fields([(quality if quality[0]!='!' else quality[1:], 'b', 200),])
    trees = {b: RootTree('bunch{0}Add'.format(b)) for b in bcids}
    for tree in trees.itervalues():
        tree.branch_f('vtx_x')
        tree.branch_f('vtx_y')
        tree.branch_f('vtx_xError')
        tree.branch_f('vtx_yError')
        tree.branch_f('vtx_nTrk')
        tree.branch_i('scanstep')
        tree.branch_i('timestamp')
    for event in chain.events(verbose):
        if event['nVtx'] <= 0:
            continue
        if -1 not in bcids:
            bcid = event['bunchCrossing']
            if bcid not in bcids:
                continue
        else:
            bcid = -1
        for scanstep, (begin, end) in enumerate(times):
            if event['timeStamp_begin'] <= begin:
                continue
            if event['timeStamp_begin'] >= end:
                continue
            break
        else:
            continue
        trees[bcid].set('scanstep', scanstep)
        trees[bcid].set('timestamp', event['timeStamp_begin'])
        for vtx in range(event['nVtx']):
            for quality in qualities:
                if quality[0]!='!':
                    if not event[quality][vtx]:
                        continue
                else:
                    if event[quality[1:]][vtx]:
                        continue
            if event['vtx_nTrk'][vtx] < mintrk:
                continue
            trees[bcid].set('vtx_x', event['vtx_x'][vtx])
            trees[bcid].set('vtx_y', event['vtx_y'][vtx])
            trees[bcid].set('vtx_xError', 0.0 if noerror else event['vtx_xError'][vtx])
            trees[bcid].set('vtx_yError', 0.0 if noerror else event['vtx_yError'][vtx])
            trees[bcid].set('vtx_nTrk', event['vtx_nTrk'][vtx])
            trees[bcid].Fill()
    return trees

def get_beamcorrection(
    name, correctionsfile
):
    """Run over created corrections file and get beam-beam corrections.

    trees: Dictionary of trees to be used for histogram creation (key unused).
    correlationfile: Json file with the beam-beam corrections.

    returns dictionary with beam corrections for a specific tree and the value of bcid.
    """

    num_beam, bunch = re.findall('\d+', name)
    if "X" in name:
        axis_beam = num_beam + 'X'
    elif "Y" in name:
        axis_beam = num_beam + 'Y'

    corr = {}
    with open(correctionsfile) as json_data:
        data = json.load(json_data)
        for Scan in data:
            for dataScan in data[Scan]:
                if dataScan['ScanName'] == axis_beam and dataScan['ScanPointNumber'] == 1:
                    corr = data[Scan]

    return corr, bunch


def get_betacorrection(
    name, betacorrectionsfile
):
    """Run over created corrections file and get dynamic beta corrections.

    trees: Dictionary of trees to be used for histogram creation (key unused).
    correlationfile: Json file with the beam-beam corrections.

    returns dictionary with dynamic beta corrections for a specific tree and the value of bcid.
    """

    num_beam, bunch = re.findall('\d+', name)
    if "X" in name:
        axis_beam = num_beam + 'X'
    elif "Y" in name:
        axis_beam = num_beam + 'Y'

    with open(betacorrectionsfile) as json_data:
        data = json.load(json_data)

    return data[axis_beam][bunch]


def make_histograms_forbcid(
    trees, nbins, mintrk, scaling=1.0, verbose=False, extracond=None, beamcorrectionsfile=None, betacorrectionsfile=None
):
    """Run over created trees and select Beam Imaging data for histograms for each single bcid.

    trees: Dictionary of trees to be used for histogram creation (key unused).
    nbins: Number of bins in each dimension.
    mintrk: Minimal number of tracks (int) for the event selection.
    scaling: Float by which the histograms are rescaled (x=xraw/scaling).
    extracond: Additional condition to be applied to selected data.
    returns dictionary with histograms.
    """
    if extracond is None:
        extracond = ''
    else:
        extracond = ' && ({0})'.format(extracond)
    hists = {}

    for i, tree in enumerate(trees.itervalues()):
        print(tree)
        name = 'hist_{0}'.format(tree.GetName())
        print(name)
        condition = 'vtx_nTrk>={0}{1}'.format(mintrk, extracond)
        print(condition)
        hist2 = TH2F(name, name, nbins, -10.0, 10.0, nbins, -10.0, 10.0)
        n1,n2 = 0,0
        draw1 = '(vtx_y)/{0}:(vtx_x)/{0}>>hnew{1}'.format(scaling, i)
        n = tree.Draw(draw1, condition, 'goff')
        print(n)
        hist1 = gDirectory.Get('hnew{0}'.format(i))
        offx = round(hist1.GetMean(1), 2)
        offy = round(hist1.GetMean(2), 2)
        if verbose:
            print '<<< {0} entries with offset {1}, {2}'.format(n, offx, offy)

        text = TLatex()
        text.SetNDC()
        text.SetTextFont(62)
        text.SetTextSize(0.0375)
        text.SetTextAlign(13)
        #Apply corrections for Beam-Beam effects
        if beamcorrectionsfile and betacorrectionsfile:
            treename = tree.GetName()
            beamcorrections, bunch = get_beamcorrection(treename, beamcorrectionsfile)
            betacorrections = get_betacorrection(treename, betacorrectionsfile)

            for beamdata in beamcorrections:

                for scan in range(19):
                    if beamdata['ScanPointNumber']-1 != scan:
                        continue
                    for betadata in betacorrections:
                        if betadata[0] == beamdata['ScanPointNumber']-1:
                            betacorr = betadata[1]

                    #Get the corrections according to the value of scanstep and bcid
                    x_corr, y_corr = 0. , 0.
                    for bcid in beamdata['corr_Xcoord']:
                        if bunch == bcid:
                            x_corr = beamdata['corr_Xcoord'][bcid]
                            y_corr = beamdata['corr_Ycoord'][bcid]

                    #Apply extraconditions for the scanstep value
                    ScanPointcondition = '({0} && scanstep == {1})*{2}'.format(condition,scan,betacorr )
                    #ScanPointcondition = '({0} && scanstep == {1})*{2}'.format(condition,scan,1 )
                    
                    #Draw corrected vtx_y:vtx_x shifted by the mean values
                    #draw2 = '(vtx_y+({4}/2))/{0}-{1}:(vtx_x+({5}/2))/{0}-{2}>>+{3}' \
                    #        .format(scaling, offy, offx, name, y_corr*10, x_corr*10)
                    draw2 = '(vtx_y+({4}/2))/{0}-{1}:(vtx_x+({5}/2))/{0}-{2}>>+{3}' \
                            .format(scaling, offy, offx, name, 0, 0)

                    tree.Draw(draw2, ScanPointcondition, 'goff')

                    #draw step by step plots
                    hist = TH2F('hist{}'.format(scan), '', 50, 0.0, 0.12, 50, 0.04, 0.16)
                    hist.GetZaxis().SetRangeUser(1.0, 400.0)
                    tree.Draw('vtx_y:vtx_x>>hist{}'.format(scan), 'vtx_nTrk>14 && scanstep=={}'.format(scan+1), 'goff')
                    plot = ColorBase(hist)
                    plot._xtitle = 'x [cm]'
                    plot._ytitle = 'y [cm]'
                    plot.draw()
                    plot.SetLogz()
                    text.DrawLatex(0.15, 0.88, '{0}, step {1}'.format(tree.GetName(),scan))
                    plot.Update()
                    plot.SaveAs('Fill4266_rereco_{0}_step{1}.png'.format(tree.GetName(),scan))
                    plot.Close()
        else:

            draw2 = 'vtx_y/{0}-{1}:vtx_x/{0}-{2}>>{3}' \
                .format(scaling, offy, offx, name)
            tree.Draw(draw2, condition, 'goff')

        hists[name] = hist2

    return hists

def make_histograms(
    hists, mergehistname, nbins
):
    """Run over created trees and select Beam Imaging data for histograms merging all bcids together

    hists: Dictionary of histograms to be used for histogram creation (key unused).
    mergehistname: Name of the final histogram where all the histograms for each bcid are merged together
    nbins: Number of bins in each dimension.
    returns dictionary with histograms.
    """

    mergedhist = TH2F(mergehistname, mergehistname, nbins, -10.0, 10.0, nbins, -10.0, 10.0)
    for i, hist in enumerate(hists.itervalues()):
        name = hist.GetName()
        if mergehistname in name:
            mergedhist.Add(hist)
    return mergedhist

def make_vdmhistos(
    trees, nbins, mintrk, stepsize,
    scaling=1.0, crange=(-10.0, 10.0), verbose=False
):
    """Run over created trees, select VdM data and use it for BI histograms.

    trees: Dictionary of trees to be used for histogram creation (key unused).
    nbins: Number of bins in each dimension.
    mintrk: Minimal number of tracks (int) for the event selection.
    stepsize: 2-tuple of step sizes of beam that is supposed to rest.
    scaling: Float by which the histograms are rescaled (x=xraw/scaling).
    crange: 2-tuple of limits of the coordinates.
    returns dictionary with histograms.
    """
    hists = {}
    for i, tree in enumerate(trees.itervalues()):
        name = 'hist_{0}'.format(tree.GetName())
        condition = 'vtx_nTrk>={0}'.format(mintrk)#' && scanstep>=2'

        draw1 = '(vtx_y{2:+.8f}*scanstep)/{0}):(vtx_x{1:+.8f}*scanstep)/{0}>>hnew{3}' \
                .format(scaling, -1.0*stepsize[0], -1.0*stepsize[1], i)
        n = tree.Draw(draw1, condition, 'goff')
        hist1 = gDirectory.Get('hnew{0}'.format(i))

        offx = round(hist1.GetMean(1), 2)
        offy = round(hist1.GetMean(2), 2)
        if verbose:
            print '<<< {0} entries with offset {1}, {2}'.format(n, offx, offy)

        draw2 = '(vtx_y{2:+.8f}*scanstep)/{0}{4:+f}:(vtx_x{1:+.8f}*scanstep)/{0}{3:+f}>>{5}' \
                .format(scaling, -1.0*stepsize[0], -1.0*stepsize[1], -1.0*offx, -1.0*offy, name)
        hist2 = TH2F(
            name, name, nbins, crange[0], crange[1], nbins, crange[0], crange[1]
        )
        tree.Draw(draw2, condition, 'goff')

        hists[name] = hist2
    return hists
