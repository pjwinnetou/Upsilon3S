#!/bin/bash
#DrawYieldHist(const double ptMin = 0., const double ptMax = 50., const double ptCut = 3.5, int cBinLow = 0, int cBinHigh = 180, in    t Trig = "S1", bool trigSwitch = true)

#0-50GeV
#Single mu pt 3.5, S1, S2
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 0, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 87, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 0, 86, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 0, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 107, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 0, 106, "S2", true)'

#root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 87, 180, "S1", false)'
#root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 107, 180, "S2", false)'
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 0, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 0, 106, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 0, 86, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 87, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 50., 3.5, 107, 180, "S13", false)'

root -l -b -q 'DrawYieldHist.C(0., 50., 4, 0, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 50., 4, 87, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 50., 4, 0, 86, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 50., 4, 0, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(0., 50., 4, 107, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(0., 50., 4, 0, 106, "S2", true)'

#root -l -b -q 'DrawYieldHist.C(0., 50., 4, 87, 180, "S1", false)'
#root -l -b -q 'DrawYieldHist.C(0., 50., 4, 107, 180, "S2", false)'
root -l -b -q 'DrawYieldHist.C(0., 50., 4, 0, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 50., 4, 0, 106, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 50., 4, 0, 86, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 50., 4, 87, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 50., 4, 107, 180, "S13", false)'


#0-10GeV
#Single mu pt 3.5, S1, S2
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 0, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 87, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 0, 86, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 0, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 107, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 0, 106, "S2", true)'

#root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 87, 180, "S1", false)'
#root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 107, 180, "S2", false)'
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 0, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 0, 106, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 0, 86, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 87, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 10.,3.5, 107, 180, "S13", false)'

root -l -b -q 'DrawYieldHist.C(0., 10.,4, 0, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 10.,4, 87, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 10.,4, 0, 86, "S1", true)'
root -l -b -q 'DrawYieldHist.C(0., 10.,4, 0, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(0., 10.,4, 107, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(0., 10.,4, 0, 106, "S2", true)'

#root -l -b -q 'DrawYieldHist.C(0., 10.,4, 87, 180, "S1", false)'
#root -l -b -q 'DrawYieldHist.C(0., 10.,4, 107, 180, "S2", false)'
root -l -b -q 'DrawYieldHist.C(0., 10.,4, 0, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 10.,4, 0, 106, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 10.,4, 0, 86, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 10.,4, 87, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(0., 10.,4, 107, 180, "S13", false)'

#10-50GeV
#Single mu pt 3.5, S1, S2
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 0, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 87, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 0, 86, "S1", true)'
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 0, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 107, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 0, 106, "S2", true)'

#root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 87, 180, "S1", false)'
#root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 107, 180, "S2", false)'
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 0, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 0, 106, "S13", false)'
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 0, 86, "S13", false)'
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 87, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(10., 50., 3.5, 107, 180, "S13", false)'

root -l -b -q 'DrawYieldHist.C(10., 50., 4, 0, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(10., 50., 4, 87, 180, "S1", true)'
root -l -b -q 'DrawYieldHist.C(10., 50., 4, 0, 86, "S1", true)'
root -l -b -q 'DrawYieldHist.C(10., 50., 4, 0, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(10., 50., 4, 107, 180, "S2", true)'
root -l -b -q 'DrawYieldHist.C(10., 50., 4, 0, 106, "S2", true)'

#root -l -b -q 'DrawYieldHist.C(10., 50., 4, 87, 180, "S1", false)'
#root -l -b -q 'DrawYieldHist.C(10., 50., 4, 107, 180, "S2", false)'
root -l -b -q 'DrawYieldHist.C(10., 50., 4, 0, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(10., 50., 4, 0, 106, "S13", false)'
root -l -b -q 'DrawYieldHist.C(10., 50., 4, 0, 86, "S13", false)'
root -l -b -q 'DrawYieldHist.C(10., 50., 4, 87, 180, "S13", false)'
root -l -b -q 'DrawYieldHist.C(10., 50., 4, 107, 180, "S13", false)'

