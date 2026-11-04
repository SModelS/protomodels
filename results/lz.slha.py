smodelsOutput = {
'OutputStatus' : {
    'sigmacut' : 0.005,
    'minmassgap' : 5.0,
    'minmassgapisr' : 1.0,
    'maxcond' : 0.2,
    'ncpus' : 1,
    'model' : 'share.models.mssm',
    'promptwidth' : 1e-11,
    'stablewidth' : 1e-25,
    'checkinput' : True,
    'doinvisible' : True,
    'docompress' : True,
    'testcoverage' : True,
    'computestatistics' : True,
    'combinesrs' : True,
    'pyhfbackend' : 'pytorch',
    'reportallsrs' : False,
    'file status' : 1,
    'decomposition status' : 1,
    'warnings' : 'Input file ok',
    'input file' : './lz.slha',
    'database version' : '3.2.1',
    'smodels version' : '3.2.1'
},
'ExptRes' : [
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.009763414,
        'upper limit (fb)' : 0.037,
        'expected upper limit (fb)' : 0.038,
        'TxNames' : ['TDTM2F'],
        'FinalStates' : ['PV > (MET,*anySM),(MET)'],
        'Mass (GeV)' : [('C1+', 120.3), ('N1/N2', 120.0), ('N1', 120.0)],
        'AnalysisID' : 'ATLAS-SUSY-2018-19',
        'DataSetID' : 'SR_EW',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 136.0,
        'dataType' : 'efficiencyMap',
        'r' : 0.2638761,
        'r_expected' : 0.256932,
        'Width (GeV)' : [
            ('C1+', 1e-15),
            ('N1/N2', 9.085e-27),
            ('N1', 'stable')
        ],
        'nll' : 2.265302,
        'nll_min' : 2.058186,
        'nll_SM' : 2.058186
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.01441532,
        'upper limit (fb)' : 0.0602,
        'expected upper limit (fb)' : 0.0983,
        'TxNames' : ['TDTM2F'],
        'FinalStates' : ['PV > (MET,*anySM),(MET)'],
        'Mass (GeV)' : [('C1+', 120.3), ('N1/N2', 120.0), ('N1', 120.0)],
        'AnalysisID' : 'CMS-EXO-19-010',
        'DataSetID' : 'SR_nlay6p',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 101.0,
        'dataType' : 'efficiencyMap',
        'r' : 0.2394572,
        'r_expected' : 0.1466462,
        'Width (GeV)' : [
            ('C1+', 1e-15),
            ('N1/N2', 9.085e-27),
            ('N1', 'stable')
        ],
        'nll' : 4.993759,
        'nll_min' : 4.501325,
        'nll_SM' : 4.501325
    },
    {
        'maxcond' : 0.0,
        'theory prediction (fb)' : 0.01623656,
        'upper limit (fb)' : 0.22,
        'expected upper limit (fb)' : 0.28,
        'TxNames' : ['TDTM2F'],
        'FinalStates' : ['PV > (MET,*anySM),(MET)'],
        'Mass (GeV)' : [('C1+', 120.3), ('N1/N2', 120.0), ('N1', 120.0)],
        'AnalysisID' : 'ATLAS-SUSY-2016-06',
        'DataSetID' : 'SR_EW',
        'AnalysisSqrts (TeV)' : 13.0,
        'lumi (fb-1)' : 36.1,
        'dataType' : 'efficiencyMap',
        'r' : 0.07380255,
        'r_expected' : 0.05798772,
        'Width (GeV)' : [
            ('C1+', 1e-15),
            ('N1/N2', 9.085e-27),
            ('N1', 'stable')
        ],
        'nll' : 4.367143,
        'nll_min' : 4.277565,
        'nll_SM' : 4.277565
    }
],
'Total xsec for missing topologies (fb)' : 0.09289028,
'missing topologies' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.09289028,
        'SMS' : 'PV > (jet,jet,MET),(ta,ta,MET)'
    }
],
'Total xsec for missing topologies with displaced decays (fb)' : 0.277587,
'missing topologies with displaced decays' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.185058,
        'SMS' : 'PV > (jet,jet,MET),(l,l,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.09252899,
        'SMS' : 'PV > (jet,jet,MET),(ta,ta,MET)'
    }
],
'Total xsec for missing topologies with prompt decays (fb)' : 0.0,
'missing topologies with prompt decays' : [],
'Total xsec for topologies outside the grid (fb)' : 0.4644458,
'topologies outside the grid' : [
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.2786653,
        'SMS' : 'PV > (MET),(jet,jet,MET)'
    },
    {
        'sqrts (TeV)' : 13.0,
        'weight (fb)' : 0.1857806,
        'SMS' : 'PV > (jet,jet,MET),(l,l,MET)'
    }
]
}