from __future__ import print_function
# investigate shower development based on RecHits and SimClusters
import ROOT
import os
import numpy as np


from HGCalImagingAlgo import *
from TestBeamNtupleDataFormat import HGCalNtuple

### Basic setup for testing (sensor dependence, treshold w.r.t noise, cartesian coordinates for multi-clustering
dependSensor = True
# 2D clustering parameters
deltac = [2., 2., 5.] # in cartesian coordiantes in cm, per detector
ecut = 3 # relative to the noise
# multi-clustering parameters
multiclusterRadii = [2., 5., 5.] # in cartesian coordiantes in cm, per detector
minClusters = 3 # request at least minClusters+1 2D clusters
# verbosity, allowed events/layers for testing/histograming, etc.
allowedRangeLayers = [] # layer considered for histograming e.g. [10, 20], empty for none
allowedRangeEvents = list(range(0,3)) # event numbers considered for histograming, e.g. [0,1,2], empty for none
verbosityLevel = 2 # 0 - only basic info (default); 1 - additional info; 2 - detailed info printed, histograms produced

def getRecHitDetIds(rechits):
    recHitsList = []
    for rHit in rechits:
        recHitsList.append(rHit.detid())
    # print("RecHits -"*10)
    # print(recHitsList)
    recHits = np.array(recHitsList)
    return recHits

def getHitList(simClus, recHitDetIds):
    sClusHitsList = []
    for DetId in simClus.hits():
        sClusHitsList.append(DetId)
    sClusHits = np.array(sClusHitsList)
    # thanks to http://stackoverflow.com/questions/11483863/python-intersection-indices-numpy-array
    recHitIndices = np.nonzero(np.in1d(recHitDetIds, sClusHits))
    return recHitIndices

# get list of rechist associated to sim-cluster hits
def getRecHitsSimAssoc(rechits_raw, simcluster):
    # get sim-cluster associations
    nSimClus = 0
    simClusHitAssoc = []
    recHitDetIds = getRecHitDetIds(rechits_raw)
    for simClusIndex, simClus in enumerate(simcluster):
        simClusHitAssoc.append(getHitList(simClus, recHitDetIds))
        nSimClus += 1
    # get list of rechist associated to simhits
    rHitsSimAssoc = [[] for k in range(0,nSimClus)]
    for simClusIndex, simCl in enumerate(simcluster):
        if (verbosityLevel>=1): print( "Sim-cluster index: ",simClusIndex, ", pT: ",simCl.pt(), ", E: ",simCl.energy(), ", phi: ",simCl.phi(), ", eta: ",simCl.eta())
        # loop over sim clusters and then rechits
        rHitsSimAssocTemp = []
        for hitIndexArray in simClusHitAssoc[simClusIndex]:
            for hitIndex in hitIndexArray:
                thisHit = rechits_raw[hitIndex]
                if(not recHitAboveThreshold(thisHit, ecut, dependSensor)[1]): continue
                # independent of sim cluster, after cleaning
                rHitsSimAssocTemp.append(thisHit)
        rHitsSimAssoc[simClusIndex]= rHitsSimAssocTemp
    return rHitsSimAssoc

# 3D histograming of rechist associated to sim-cluster hits ("event displays")
def histRecHitsSimAssoc(rHitsSimAssoc, currentEvent, histDict, tag = "rHitsAssoc_", zoomed = False):
    # sanity check
    if (histDict == None): return

    # define event-level hists
    if (zoomed): # zoomed for testing/convenience around the eta/phi of most energetic hit
        rs = sorted(range(len(rHitsSimAssoc[0])), key=lambda k: rHitsSimAssoc[0][k].energy(), reverse=True) # indices sorted by decreasing rho
        c_phi = rHitsSimAssoc[0][rs[0]].phi()
        c_eta = rHitsSimAssoc[0][rs[0]].eta()
        histDict[tag+"map_lay_x_y_evt{}".format(currentEvent)]  = ROOT.TH3F(tag+"map_lay_x_y_evt{}".format(currentEvent), tag+"map_lay_x_y_evt{};z (cm);x (cm);y (cm)".format(currentEvent), 55, 0, 55, 50, c_phi-0.5, c_phi+0.5, 50, c_eta-0.5, c_eta+0.5) # 3D rechists associated to sim-cluster (with ecut cleaning)
    else:
        histDict[tag+"map_lay_x_y_evt{}".format(currentEvent)]  = ROOT.TH3F(tag+"map_lay_x_y_evt{}".format(currentEvent), tag+"map_lay_x_y_evt{};z (cm);x (cm);y (cm)".format(currentEvent), 55, 0, 55, 32, -8, 8, 32, -8, 8) # 3D rechists associated to sim-cluster (with ecut cleaning)

    for simClusIndex in range(0,len(rHitsSimAssoc)):
        # define sim-cluster-level hists
        histDict[tag+"map_lay_x_y_evt{}_sim{}".format(currentEvent, simClusIndex)]  = ROOT.TH3F(tag+"map_lay_x_y_evt{}_sim{}".format(currentEvent, simClusIndex), tag+"map_lay_x_y_evt{}_sim{};z (cm);x (cm);y (cm)".format(currentEvent, simClusIndex), 55, 0, 55, 32, -8, 8, 32, -8, 8)
        # loop over assoc. rec hits
        for thisHit in rHitsSimAssoc[simClusIndex]:
            histDict[tag+"map_lay_x_y_evt{}".format(currentEvent)].Fill(abs(thisHit.z()), thisHit.phi(), thisHit.eta()) # for each sim cluster, after cleaning
            if(thisHit.energy < ecut): continue
            histDict[tag+"map_lay_x_y_evt{}_sim{}".format(currentEvent, simClusIndex)].Fill(abs(thisHit.z()), thisHit.phi(), thisHit.eta()) # independent of sim cluster, before cleaning

    return histDict

# 2D histograming of rechists in the chosen layerts, given by allowedRangeLayers
def histRecHits(rHits, currentEvent, histDict, tag = "rHits_", zoomed = False):
    # sanity check
    if (histDict == None): return

    # define hists per layer
    for layer in range(1, 41):
        if (layer in allowedRangeLayers): # testing limitation
            if (zoomed): # zoomed for testing/convenience around the eta/phi of most energetic hit
                rs = sorted(range(len(rHits)), key=lambda k: rHits[k].energy(), reverse=True) # indices sorted by decreasing rho
                c_phi = rHits[rs[0]].phi()
                c_eta = rHits[rs[0]].eta()
                histDict[tag+"eng_eta_phi_evt{}_lay{}".format(currentEvent, layer)]  = ROOT.TH2F(tag+"eng_eta_phi_evt{}_lay{}".format(currentEvent, layer), tag+"eng_eta_phi_evt{}_lay{};x (cm);y (cm)".format(currentEvent, layer), 40, c_eta-0.1, c_eta+0.1, 40, c_phi-0.1, c_phi+0.1) # 2D energy-weighted-map of raw rechists (with ecut cleaning)
            else:
                histDict[tag+"eng_eta_phi_evt{}_lay{}".format(currentEvent, layer)]  = ROOT.TH2F(tag+"eng_eta_phi_evt{}_lay{}".format(currentEvent, layer), tag+"eng_eta_phi_evt{}_lay{};x (cm);y (cm)".format(currentEvent, layer), 32, -8, 8, 32, -8, 8) # 2D energy-weighted-map of raw rechists (with ecut cleaning)

    # loop over all raw rechits and fill per layer
    for rHit in rHits:
        if (rHit.layer() in allowedRangeLayers): # testing limitation
            if(rHit.energy() < ecut): continue
            histDict[tag+"eng_eta_phi_evt{}_lay{}".format(currentEvent, rHit.layer())].Fill(rHit.eta(), rHit.phi(), rHit.energy())

    return histDict

# 2D histograming of the clustered rechist with stand-alone algo, weighted by energy
def histHexelsClustered(hexelsClustered, currentEvent, histDict, tag = "clustHex_", zoomed = False):
    # sanity check
    if (histDict == None): return

    # define event-level hists
    if (zoomed): # zoomed for testing/convenience around the eta/phi of most energetic hit
        rs = sorted(range(len(hexelsClustered)), key=lambda k: hexelsClustered[k].weight, reverse=True) # indices sorted by decreasing rho
        c_phi = hexelsClustered[rs[0]].phi
        c_eta = hexelsClustered[rs[0]].eta
        histDict[tag+"eng_phi_eta_evt{}".format(currentEvent)]  = ROOT.TH3F(tag+"eng_phi_eta_evt{}".format(currentEvent), tag+"eng_phi_eta_evt{};z (cm);x (cm);y (cm)".format(currentEvent), 55, 0, 55, 80, c_phi-0.8, c_phi-0.8, 80, c_eta-0.8, c_eta-0.8) # 3D rechists clustered with algo (with ecut cleaning)
    else:
        histDict[tag+"eng_phi_eta_evt{}".format(currentEvent)]  = ROOT.TH3F(tag+"eng_phi_eta_evt{}".format(currentEvent), tag+"eng_phi_eta_evt{};z (cm);x (cm);y (cm)".format(currentEvent), 55, 0, 55, 32, -8, 8, 32, -8, 8) # 3D rechists clustered with algo (with ecut cleaning)

    # loop over all clustered rechist
    for iNode in hexelsClustered:
        histDict[tag+"eng_phi_eta_evt{}".format(currentEvent)].Fill(abs(iNode.z), iNode.x, iNode.y, iNode.weight)

    return histDict

# 1D histograming of given list of values
def histValue1D(fValues, histDict, tag = "hist1D_", title = "hist 1D", axunit = "a.u.", binsRangeList = [10, -1, 1], ayunit = "a.u."):
    # sanity check
    if (histDict == None): return

    # define event-level hists
    histDict[tag]  = ROOT.TH1F(tag, title+";"+axunit+";"+ayunit, binsRangeList[0], binsRangeList[1], binsRangeList[2])
    histDict[tag].GetYaxis().SetTitleOffset(histDict[tag].GetYaxis().GetTitleOffset()*1.5)
    # loop over all values
    if (verbosityLevel>=3): print( "tag: ", tag, ", fValues: ", fValues)
    for value in fValues:
        histDict[tag].Fill(value)
    return histDict

# print/save all histograms
def histPrintSaveAll(histDict, outDir):
    imgType = "pdf"
    canvas = ROOT.TCanvas(outDir, outDir, 500, 500)
    if (verbosityLevel>=3): print( "histDict.items(): ", histDict.items())
    for key, item in histDict.items():
        # do not save empty histograms
        if (type(item) == ROOT.TH1F) or (type(item) == ROOT.TH2F) or (type(item) == ROOT.TH3F):
            if item.GetEntries() == 0:
                continue
        ROOT.gStyle.SetPalette(ROOT.kBird)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPadTopMargin(0.05)
        ROOT.gStyle.SetPadBottomMargin(0.12)
        ROOT.gStyle.SetPadLeftMargin(0.15)
        ROOT.gStyle.SetPadRightMargin(0.02)
        if type(item) == ROOT.TH1F:
            item.Draw("hist0")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        if type(item) == ROOT.TH2F:
            item.Draw("colz")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        elif type(item) == ROOT.TH3F:
            item.Draw("box")
            canvas.SaveAs("{}/{}.{}".format(outDir, key, imgType))
        else:
            continue
    return

def main():
    # init output stuff
    outDir = "testClustering"
    if not os.path.exists(outDir): os.makedirs(outDir)
    histDict = {}

    # get sample/tree
    # please give an CMSSW930 NTUP root file.
    #########################################
    # ntuple = HGCalNtuple("root://eoscms.cern.ch//eos/cms/store/cmst3/group/hgcal/CMG_studies/Production/FlatRandomEGunProducer_pdgid211_E20_cmssw93X_withPRs_20170817/NTUP/partGun_PDGid211_x100_E20.0To20.0_NTUP_1.root") # CMSSW_9_3_0_pre3 with some pre4 PRs on top
    ntuple = HGCalNtuple("../ntuple_100.root") # cmssw921 with all recent fixes as of June 12

    # prepare some lists for comparions
    tot_nClust2D_rerun = []
    clusters2D_eng_rerun = []
    clusters2DMultiSelected_eng_rerun = []

    # start event loop
    for event in ntuple:
        if (not event.entry() in allowedRangeEvents): continue # checking external condition
        if (verbosityLevel>=1): print( "\nCurrent event: ", event.entry())

        # get collections of raw rechits, sim clusters, 2D clusters, multi clusters, etc.
        recHitsRaw = event.recHits()

        rHitsCleaned = [rechit for rechit in recHitsRaw if recHitAboveThreshold(rechit, ecut, dependSensor)[1]]

        ### Imaging algo run as stand-alone (python)
        # instantiate the stand-alone clustering implemented in HGCalImagingAlgo
        HGCalAlgo = HGCalImagingAlgo(ecut = ecut, deltac = deltac, multiclusterRadii = multiclusterRadii, minClusters = minClusters, dependSensor = dependSensor, verbosityLevel = 0)
        # produce 2D clusters with stand-alone algo, out of all raw rechits
        clusters2D_rerun = HGCalAlgo.makeClusters(recHitsRaw) # nested list of "hexels", per layer, per 2D cluster
        # produce multi-clusters with stand-alone algo, out of all 2D clusters
        multiClustersList_rerun = HGCalAlgo.make3DClusters(clusters2D_rerun) # flat list of multi-clusters (as basic clusters)

        # get for testing: flat list of 2D clustered, and flat list of clustered non-halo "hexeles" (from stand-alone algo)
        clusters2DList_rerun = HGCalAlgo.getClusters(clusters2D_rerun, verbosityLevel = 0) # flat list of 2D clusters (as basic clusters)
        hexelsClustered_rerun = [iNode for bClust in clusters2DList_rerun for iNode in bClust.thisCluster if not iNode.isHalo]  # flat list of clustered "hexeles", without the "halo" hexels

        ### Produce some basic histograms for each event (2D/3D view of associated sim-clusters, selected rec-hits, etc.)
        if (verbosityLevel>=2):
            # histograming of raw rechist (with ecut cleaning)
            histDict = histRecHits(rHitsCleaned, event.entry(), histDict, tag = "rHitsCleaned_", zoomed = True)
            # histograming of clustered hexels
            histDict = histHexelsClustered(hexelsClustered_rerun, event.entry(), histDict, tag = "clustHex_", zoomed = False)

        rHitsClustdDID = [iNode.detid for iNode in hexelsClustered_rerun] # list of detids for clustered hexels
        # print some info if requested
        if (verbosityLevel>=1):
            print( "num of rechits clustered with imaging algo. : ", len (rHitsClustdDID))

        ### Compare stand-alone and reco-level clustering
        clusters2DListMultiSelected_rerun = [cls for multiCluster in multiClustersList_rerun for cls in multiCluster.thisCluster]
        # print more details if requested
        if (verbosityLevel>=1):
            ls = sorted(range(len(clusters2DListMultiSelected_rerun)), key=lambda k: clusters2DListMultiSelected_rerun[k].thisCluster[0].layer, reverse=False) # indices sorted by increasing layer number
            for index in range(len(multiClustersList_rerun)): print( "Multi-cluster (RE-RUN) index: ", index, ", No. of 2D-clusters = ", len(multiClustersList_rerun[index].thisCluster), ", Energy  = ", multiClustersList_rerun[index].energy, ", x = ", multiClustersList_rerun[index].x, ", y = ", multiClustersList_rerun[index].y, ", z = ", multiClustersList_rerun[index].z )
            print( "num of clusters2D re-run: ", len(clusters2DListMultiSelected_rerun))
            print( "num of multi-cluster re-run: ", len(multiClustersList_rerun))

        ### Produce some basic histograms with general info (one per sample)
        if (verbosityLevel>=2):
            # number of 2D clusters from algo at re-run step
            tot_nClust2D_rerun.append(len(clusters2DListMultiSelected_rerun))
            clusters2D_eng_rerun.extend([clusters2DList_rerun[k].energy for k in range(0,len(clusters2DList_rerun))]) # eng re-run

    # histograms - 2D clusters counting
    histDict = histValue1D(tot_nClust2D_rerun, histDict, tag = "tot_nClust2D_rerun", title = "Rerun: total Num(2D clusters)",  axunit = "N_{cl.2D}",    binsRangeList = [100, 0, 100], ayunit = "total occurences")
    # histograms - 2D clusters energy spectra
    histDict = histValue1D(clusters2D_eng_rerun,   histDict, tag = "Clust2D_Eng_Rerun",   title = "Rerun E(all 2D clusters)",    axunit = "E_{cl.2D}[GeV]",  binsRangeList = [100, 0, 1000], ayunit = "N(2D clusters)")
    # print/save histograms
    histPrintSaveAll(histDict, outDir)

if __name__ == '__main__':
    main()
