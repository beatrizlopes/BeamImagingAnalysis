from itertools import product
from json import load
from sys import argv

from lib.submit import QSubmitter as submit

ma = 'beatriz.lopes@desy.de'
mo = False
number = 20 #100
test = False

# for shapeFitter, computeCorr, integrateResiduals
switchConfignummodel = False
configs = [
    '6868_rereco_second',
    '6868_rereco_second_verytight',
    '6868_rereco_second_drift1',
    '6868_rereco_second_drift2',
]
modelversion = [
    ('SupDG', 'v4'),
    ('SupDG', 'v3'),
    ('SupG', 'v2'),
    ('SG', 'v1'), ('noCorr', 'v1'), ('DG', 'v1'), ('TG', 'v1'),
]
confignummodel = [
    ('6868_second', 2, 'TG', 'v1'),
    #('6868_second', 4, 'TG', 'v1'),
]

# for closureTest
n = [10, 7, 5, 3]
toymodels = ['SupDG']
fitmodels = ['SupDG']
toyfitmodels = product(toymodels, fitmodels)

def shapeFitter(config, i, model, bcid, json, version):
    args = [json, str(bcid), model]
    names = ['sF'] + config.split('_') + [str(i), model, version]
    submit(
        'res/jobs/shapeFitter_{0}.sh'.format(version), args=args, names=names,
        ma=ma, mo=mo, rt=5, test=test
    )

def computeCorr(config, i, model, bcid, json, version):
    args = [json, str(bcid), model, version]
    names = ['cC'] + config.split('_') + [str(i), model, version]
    submit(
        'res/jobs/computeCorr.sh', args=args, names=names, ma=ma, mo=mo, rt=3
    )

def integrateResiduals(config, i, model, bcid, json, version):
    args = [json, str(bcid), model, version]
    names = ['iR'] + config.split('_') + [str(i), model, version]
    submit(
        'res/jobs/integrateResiduals.sh', args=args, names=names, ma=ma, mo=mo,
        vmem='5G', rt=12
    )

for arg in argv[1:]:
    if arg in ('shapeFitter', 'computeCorr', 'integrateResiduals'):
        thejob = locals()[arg]
        if not switchConfignummodel:
            for config in configs:
                json = '{0}/res/hist/Fill{1}.json'.format(submit.cwd(), config)
                with open(json) as f:
                    bcids = load(f)['bcids']
                for i, bcid in enumerate(bcids):
                    for model, version in modelversion:
                        thejob(config, i, model, bcid, json, version)
        else:
            for config, i, model, version in confignummodel:
                json = '{0}/res/hist/Fill{1}.json'.format(submit.cwd(), config)
                with open(json) as f:
                    bcid = load(f)['bcids'][i]
                thejob(config, i, model, bcid, json, version)

    elif arg.startswith('makeHtml'):
        fill = arg[8:]
        names = ['html', fill]
        submit(
            'res/jobs/makeHtml_{0}.sh'.format(fill), names=names, ma=ma, mo=mo,
            rt=1
        )

    elif arg in ('shapeFitter_v3', 'shapeFitter_v4'):
        for config in configs:
            json = '{0}/res/hist/Fill{1}.json'.format(submit.cwd(), config)
            with open(json) as f:
                bcids = load(f)['bcids']
            for i, bcid in enumerate(bcids):
                mymo = mo
                for model, version in modelversion:
                    if version != arg[-2:]:
                        continue
                    for j in range(10, 10+number):
                        if j == 5:
                            mymo = False
                        args = [json, str(bcid), model, str(j)]
                        names = ['sF', str(j)] + config.split('_') + [
                            str(i), model
                        ]
                        submit(
                            'res/jobs/{0}.sh'.format(arg), args=args,
                            names=names, ma=ma, mo=mymo, rt=5
                        )

    elif arg == 'closureTest':
        for i, then in enumerate(n):
            for toy, fit in toyfitmodels:
                mymo = mo
                for j in range(number):
                    if j == 5:
                        mymo = False
                    numb = 100+i*number+j
                    args = [
                        toy, fit, 'basic_{0}'.format(numb), str(0.7), str(n[i])
                    ]
                    names = ['cT', toy, fit, str(numb), version]
                    submit(
                        'res/jobs/closureTest_{0}.sh'.format(version),
                        args=args, names=names, ma=ma, mo=mymo, rt=then*8,
                        vmem='5G'
                    )
