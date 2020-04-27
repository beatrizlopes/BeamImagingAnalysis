from json import dump, load
from os import listdir, mkdir
from os.path import exists
from string import replace
from sys import argv

from ROOT import TChain, TH1F

from lib.io import BareRootFile, NamedFloat, NamedString, Timestamp
from lib.prepare import make_histograms, make_histograms_forbcid, make_vdmhistos

def prepare_histograms(
    configfile, outputpath, suffix, nbins, mintrk, scaling=1.0, verbose=False,
    stepsize=None, stepsize1Y=None, stepsize2X=None, stepsize2Y=None,
    extracond=None, singlepair=False, beamcorrectionsfile=None, betacorrectionsfile=None
):
    scans = ['1X', '1Y', '2X', '2Y']
    if singlepair:
        sourcescan = {'1X':'1X', '1Y':'1Y', '2X':'1X', '2Y':'1Y'}
    else:
        sourcescan = {'1X':'1X', '1Y':'1Y', '2X':'2X', '2Y':'2Y'}
    with open(configfile) as f:
        config = load(f)
    fill = config['fill']
    version = config['version']
    bcids = config['bcids']
    heavyion = bool('heavyion' in config and config['heavyion'])
    name = 'Fill{0}_{1}'.format(fill, version)
    xerror = TH1F('xerror', '', 100, 0.0, 0.03)
    yerror = TH1F('yerror', '', 100, 0.0, 0.03)
    histograms = []
    if stepsize is None:
        steps = False
    elif None in (stepsize1Y, stepsize2X, stepsize2Y):
        steps = {
            '1X': (stepsize[0], 0.0), '1Y': (0.0, stepsize[1]),
            '2X': (-stepsize[0], 0.0), '2Y': (0.0, -stepsize[1])
        }
    else:
        steps = {
            '1X': stepsize, '1Y': stepsize1Y, '2X': stepsize2X, '2Y': stepsize2Y
        }
    if heavyion:
        crange = (-20.0, 20.0)
    else:
        crange = (-10.0, 10.0)
    for scan in scans:
        filename = '{0}/{1}_{2}.root'.format(outputpath, name, sourcescan[scan])
        if not exists(filename):
            filename = '{0}/{1}_{2}_file*.root'.format(outputpath, name, sourcescan[scan])
        trees = {}
        for bcid in bcids:
            treename = 'Beam{0}Move{1}_bunch{2}Add' \
                       .format(sourcescan[scan][0], sourcescan[scan][1], bcid)
            trees[bcid] = TChain(treename)
            trees[bcid].Add(filename)
            print trees[bcid].GetEntries()
            if singlepair and scan != sourcescan[scan]:
                trees[bcid].SetName(replace(trees[bcid].GetName(), 'Beam1', 'Beam2'))
        if stepsize:
            hists = make_vdmhistos(
                trees, nbins, mintrk, steps[scan],
                scaling=scaling, crange=crange, verbose=verbose
            )
        else:
            hists = make_histograms_forbcid(
                trees, nbins, mintrk, scaling=scaling, verbose=verbose,
                extracond=extracond, beamcorrectionsfile=beamcorrectionsfile, betacorrectionsfile=betacorrectionsfile
            )
            mergehistname = 'Beam{0}Move{1}' \
                       .format(sourcescan[scan][0], sourcescan[scan][1])
            mergedhists = make_histograms(
                hists, mergehistname, nbins
            )
            mergedhists.SetDirectory(0)
            histograms.append(mergedhists)
        for hist in hists.itervalues():
            hist.SetDirectory(0)
            histograms.append(hist)

        if not singlepair or not scan.startswith('2'):
            condition = 'vtx_nTrk>={0}'.format(mintrk)
            errx = TH1F('errx', '', 100, 0.0, 0.03)
            erry = TH1F('erry', '', 100, 0.0, 0.03)
            for tree in trees.itervalues():
                tree.Draw('vtx_xError>>+errx', condition, 'goff')
                tree.Draw('vtx_yError>>+erry', condition, 'goff')
            xerror.Add(errx)
            yerror.Add(erry)
            errx.Delete()
            erry.Delete()
    output = {
        'fill': fill,
        'name': '{0}_{1}'.format(name, suffix),
        'bcids': bcids,
        'datapath': outputpath,
        'dataname': '{0}_{1}'.format(name, suffix),
        'vtxresx': round(xerror.GetMean(), 6),
        'vtxresy': round(yerror.GetMean(), 6),
        'scaling': scaling,
        'mintrk': mintrk,
        'nbins': nbins,
        'heavyion': bool('heavyion' in config and config['heavyion'])
    }
    if extracond is not None:
        output['extracond'] = extracond
    filename1 = 'res/hist/{0}.json'.format(output['name'])
    if not exists('res'):
        mkdir('res')
    if not exists('res/hist'):
        mkdir('res/hist')
    with open(filename1, 'w') as f:
        dump(output, f, indent=4, separators=(',',': '))
    filename2 = '{0}/{1}.root'.format(outputpath, output['name'])
    with BareRootFile(filename2, 'RECREATE') as f:
        for hist in histograms:
            hist.Write()
        xerror.Write()
        yerror.Write()
        Timestamp().Write()
        for n, v in (
            ('stepsize1X', stepsize), ('stepsize1Y', stepsize1Y),
            ('stepsize2X', stepsize2X), ('stepsize2Y', stepsize2Y)
        ):
            if v is not None:
                NamedFloat(n+'_x', v[0]).Write()
                NamedFloat(n+'_y', v[1]).Write()
        NamedString('name', output['name']).Write()
    return filename2

def main():
    if len(argv) < 2 or not argv[1] or not exists(argv[1]):
        raise RuntimeError('Specify 1st argument: JSON config file.')
    configfile = argv[1]
    if len(argv) < 3 or not argv[2] or not exists(argv[2]):
        raise RuntimeError('Specify 2nd argument: Output directory.')
    outputpath = argv[2]
    if outputpath.endswith('/'):
        outputpath = outputpath[:-1]
    if len(argv) < 4 or not argv[3] or argv[3] not in ['many', 'some', 'few', 'fewer']:
        raise RuntimeError('Specify 3rd argument: Bins (many, some, or few).')
    binning = argv[3]
    nbins = {'many': 760, 'some': 190, 'few': 95, 'fewer': 50}[binning]
    if len(argv) < 5 or not argv[4] or argv[4] not in ['l', 'm', 't', 'vt', 'et']:
        raise RuntimeError(
            'Specify 4th argument: Track selection (l, m, t, vt, et).'
        )
    selection = {
        'l': 'loose', 'm': 'medium', 't': 'tight', 'vt': 'verytight',
        'et': 'extratight'
    }[argv[4]]
    mintrk = {
        'l': 10, 'm': 14, 't': 18, 'vt': 24, 'et': 50
    }[argv[4]]
    if len(argv) < 6 or not argv[5]:
        raise RuntimeError('Specify 5th argument: scaling (float).')
    try:
        scaling = float(argv[5])
    except ValueError:
        raise RuntimeError('Optional 5th argument: scaling (float).')
    suffix = '{0}_{1}'.format(binning, selection)
    if len(argv) < 7 or not argv[6] or not argv[6] in (
        'vdm', 'drift', 'extra', 'corr', 'singlevdm', 'mixed',
    ):
        prepare_histograms(
            configfile, outputpath, suffix, nbins, mintrk,
            scaling=scaling, verbose=True
        )
    elif argv[6] == 'extra':
        extracond = 'scanstep>2 && scanstep<16'
        suffix = '{0}_extra'.format(suffix)
        prepare_histograms(
            configfile, outputpath, suffix, nbins, mintrk,
            scaling=scaling, verbose=True, extracond=extracond
        )
    elif argv[6] == 'corr':
        if not argv[7] :
            raise RuntimeError('Specify 7th argument: Json file with beam beam correlations .')
        if not argv[8] :
            raise RuntimeError('Specify 8th argument: Json file with dynamic beta correlations .')
        prepare_histograms(
            configfile, outputpath, suffix, nbins, mintrk,
            scaling=scaling, verbose=True, beamcorrectionsfile = argv[7], betacorrectionsfile = argv[8]
        )
    elif argv[6] == 'vdm':
        if len(argv) < 9 or not argv[7] or not argv[8]:
            raise RuntimeError('Specify 7th, 8th argument: VdM scan step size.')
        try:
            stepsize = float(argv[7]), float(argv[8])
        except ValueError:
            raise RuntimeError('Specify 7th, 8th argument: VdM scan step size.')
        prepare_histograms(
            configfile, outputpath, suffix, nbins, mintrk,
            scaling=scaling, verbose=True, stepsize=stepsize
        )
    elif argv[6] == 'singlevdm':
        if len(argv) < 9 or not argv[7] or not argv[8]:
            raise RuntimeError('Specify 7th, 8th argument: VdM scan step size.')
        try:
            stepsize = float(argv[7]), float(argv[8])
        except ValueError:
            raise RuntimeError('Specify 7th, 8th argument: VdM scan step size.')
        prepare_histograms(
            configfile, outputpath, suffix, nbins, mintrk,
            scaling=scaling, verbose=True, stepsize=stepsize, singlepair=True
        )
    else: # 'drift', 'mixed'
        if (
            len(argv) < 15 or not argv[7] or not argv[8] or not argv[9]
            or not argv[10] or not argv[11] or not argv[12] or not argv[13]
            or not argv[14]
        ):
            raise RuntimeError('Specify 7th to 14th argument: Scan step sizes.')
        try:
            stepsize = float(argv[7]), float(argv[8])
            stepsize1Y = float(argv[9]), float(argv[10])
            stepsize2X = float(argv[11]), float(argv[12])
            stepsize2Y = float(argv[13]), float(argv[14])
        except:
            raise RuntimeError('Specify 7th to 14th argument: Scan step sizes.')
        suffix = '_'.join((suffix, argv[6]))
        prepare_histograms(
            configfile, outputpath, suffix, nbins, mintrk, scaling=scaling,
            verbose=True, stepsize=stepsize, stepsize1Y=stepsize1Y,
            stepsize2X=stepsize2X, stepsize2Y=stepsize2Y
        )

if __name__ == '__main__':
    main()
