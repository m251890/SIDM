"""Define all available histograms

All hists are defined as Histogram objects whose axes are given as a list of Axis objects, which
bundle a hist.axis with a function that defines how the axis will be filled. The underlying
hist.Hists storage is weight unless otherwise specified.
"""

# python
import math
import importlib
# columnar analysis
import hist
import awkward as ak
# local
from sidm.tools import histogram as h
from sidm.tools.utilities import dR, lxy, matched
from sidm.definitions.objects import derived_objs
# always reload local modules to pick up changes during development
importlib.reload(h)


counter_defs = {
    "Total LJs": lambda objs: ak.count(objs["ljs"].pt),
    "Gen As to muons": lambda objs: ak.count(objs["genAs_toMu"].pt),
    "Gen As to electrons": lambda objs: ak.count(objs["genAs_toE"].pt),
    "Matched gen As to muons": lambda objs: ak.count(derived_objs["genAs_toMu_matched_lj"](objs, 0.4).pt),
    "Matched gen As to electrons": lambda objs: ak.count(derived_objs["genAs_toE_matched_lj"](objs, 0.4).pt),
}

hist_defs = {
    # pv
    "pv_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 100, name="pv_n"),
                   lambda objs, mask: ak.num(objs["pvs"])),
        ],
    ),
    "pv_ndof": h.Histogram(
        [
            h.Axis(hist.axis.Regular(25, 0, 100, name="pv_ndof"),
                   lambda objs, mask: objs["pvs"].ndof),
        ],
    ),
    "pv_z": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, -50, 50, name="pv_z"),
                   lambda objs, mask: objs["pvs"].z),
        ],
    ),
    "pv_rho": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, -0.5, 0.5, name="pv_rho"),
                   lambda objs, mask: objs["pvs"].rho),
        ]
    ),
    
    #GSFelectron: Plottting electron ID varaiables, specifically for electrons within .5 of a dark photon. Also plotting 2D hists of the leading electron ID variables in barrel within Delta R < .5 of a dark photon vs the lxy of the dark photon 
    
     "electron_GsfEleDEtaInSeedCut": h.Histogram(
        [ 
            h.Axis(hist.axis.Regular(60, 0, .0065, name="electron_GsfEleDEtaInSeedCut"),
                   lambda objs, mask: objs["electrons"][(dR(objs["electrons"], objs["genAs_toE"]) < 0.5)].GsfEleDEtaInSeedCut_0),
        ],
    ),
     "electron_GsfEleDEtaInSeedCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(500, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]), 
            #
            h.Axis(hist.axis.Regular(60, 0, .0065, name="electron_GsfEleDEtaInSeedCut"),
                   lambda objs, mask: objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5][mask, 0:1].GsfEleDEtaInSeedCut_0)  
             
        ],
         evt_mask = lambda objs: ak.num(objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5])>0

     ),
     "electron_GsfEleDPhiInCut": h.Histogram(
        [ 
            h.Axis(hist.axis.Regular(60, 0, .0450, name="electron_GsfEleDPhiInCut"),
                   lambda objs, mask: objs["electrons"][(dR(objs["electrons"], objs["genAs_toE"]) < 0.5)].GsfEleDPhiInCut_0),
        ],
    ),
    "electron_GsfEleDPhiInCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(500, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]), 
            
            h.Axis(hist.axis.Regular(80, 0, .09, name="electron_GsfEleDPhiInCut"),
                   lambda objs, mask: objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5][mask, 0:1].GsfEleDPhiInCut_0)  
             
        ],
         evt_mask = lambda objs: ak.num(objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5])>0

     ),
    "electron_GsfEleEInverseMinusPInverseCut": h.Histogram(
        [
            h.Axis(hist.axis.Regular(45, 0, .3, name="electron_GsfEleEInverseMinusPInverseCut"),
                   lambda objs, mask: objs["electrons"][(dR(objs["electrons"], objs["genAs_toE"]) < 0.5)].GsfEleEInverseMinusPInverseCut_0),
        ],
    ), 
    "electron_GsfEleEInverseMinusPInverseCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(500, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]), 
            
            h.Axis(hist.axis.Regular(80, 0, .3, name="electron_GsfEleEInverseMinusPInverseCut"),
                   lambda objs, mask: objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5][mask, 0:1].GsfEleEInverseMinusPInverseCut_0)  
             
        ],
         evt_mask = lambda objs: ak.num(objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5])>0
     ),
    "electron_GsfEleRelPFIsoScaledCut": h.Histogram(
        [
            h.Axis(hist.axis.Regular(40, 0, .2, name="electron_GsfEleRelPFIsoScaledCut"),
                   lambda objs, mask: objs["electrons"][(dR(objs["electrons"], objs["genAs_toE"]) < 0.5)].GsfEleRelPFIsoScaledCut_0),
        ],
    ),
     "electron_GsfEleRelPFIsoScaledCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(500, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]), 
            
            h.Axis(hist.axis.Regular(40, 0, .2, name="electron_GsfEleRelPFIsoScaledCut"),
                   lambda objs, mask: objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5][mask, 0:1].GsfEleRelPFIsoScaledCut_0)  
             
        ],
         evt_mask = lambda objs: ak.num(objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5])>0
     ),
    "electron_GsfEleFull5x5SigmaIEtaIEtaCut": h.Histogram(
        [ 
            h.Axis(hist.axis.Regular(35, 0, .0450, name="electron_GsfEleFull5x5SigmaIEtaIEtaCut"),
                   lambda objs, mask: objs["electrons"][(dR(objs["electrons"], objs["genAs_toE"]) < 0.5)].GsfEleFull5x5SigmaIEtaIEtaCut_0),
        ],
    ),
    "electron_GsfEleFull5x5SigmaIEtaIEtaCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(500, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]), 
            
            h.Axis(hist.axis.Regular(60, 0, .02, name="electron_GsfEleFull5x5SigmaIEtaIEtaCut"),
                   lambda objs, mask: objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5][mask, 0:1].GsfEleFull5x5SigmaIEtaIEtaCut_0)  
             
        ],
         evt_mask = lambda objs: ak.num(objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5])>0
     ),
    "electron_GsfEleConversionVetoCut": h.Histogram(
        [ 
            h.Axis(hist.axis.Regular(2, 0, 2, name="electron_GsfEleConversionVetoCut"),
                   lambda objs, mask: objs["electrons"][(dR(objs["electrons"], objs["genAs_toE"]) < 0.5)].GsfEleConversionVetoCut_0),
        ],
    ),
    "electron_GsfEleConversionVetoCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(500, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]), 
            
            h.Axis(hist.axis.Regular(2, 0, 2, name="electron_GsfEleConversionVetoCut"),
                   lambda objs, mask: objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5][mask, 0:1].GsfEleConversionVetoCut_0)  
             
        ],
         evt_mask = lambda objs: ak.num(objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5])>0
     ),
    "electron_GsfEleHadronicOverEMEnergyScaledCut": h.Histogram(
         [
             h.Axis(hist.axis.Regular(30, 0, .15, name="electron_GsfEleHadronicOverEMEnergyScaledCut"), 
                    lambda objs, mask: objs["electrons"][(dR(objs["electrons"], objs["genAs_toE"]) < 0.5)].GsfEleHadronicOverEMEnergyScaledCut_0),
         ],
     ),
    "electron_GsfEleHadronicOverEMEnergyScaledCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(500, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]), 
            
            h.Axis(hist.axis.Regular(30, 0, .15, name="electron_GsfEleHadronicOverEMEnergyScaledCut"),
                   lambda objs, mask: objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5][mask, 0:1].GsfEleHadronicOverEMEnergyScaledCut_0)  
             
        ],
         evt_mask = lambda objs: ak.num(objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5])>0
     ),
    "electron_GsfEleMissingHitsCut": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="electron_GsfEleMissingHitsCut"),
                   lambda objs, mask: objs["electrons"][(dR(objs["electrons"], objs["genAs_toE"]) < 0.5)].GsfEleMissingHitsCut_0),
        ],
    ),
       "electron_GsfEleMissingHitsCut2d": h.Histogram(
        [  # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(500, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])[mask]), 
            
            h.Axis(hist.axis.Regular(10, 0, 10, name="electron_GsfEleMissingHitsCut"),
                   lambda objs, mask: objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5][mask, 0:1].GsfEleMissingHitsCut_0)  
             
        ],
         evt_mask = lambda objs: ak.num(objs["electrons"][dR(objs["electrons"], objs["genAs_toE"]) < 0.5])>0
     ),
  

    # pfelectron
    "electron_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="electron_n"),
                   lambda objs, mask: ak.num(objs["electrons"])),
        ],
    ),
    "electron_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="electron_pt"),
                   lambda objs, mask: objs["electrons"].pt),
        ],
    ),
    "electron_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="electron_eta"),
                   lambda objs, mask: objs["electrons"].eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="electron_phi"),
                   lambda objs, mask: objs["electrons"].phi),
        ],
    ),
    "electron_nearGenA_n": h.Histogram(
        [
            # number of electrons within dR=0.5 of a genA that decays to electrons
            h.Axis(hist.axis.Integer(0, 10, name="electron_nearGenA_n"),
                   lambda objs, mask: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5))),
        ],
    ),
  
     "electron_nearGenA_n_genA_lxy": h.Histogram(
        [
            # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])),
            # number of electrons within dR=0.5 of a genA that decays to electrons
            h.Axis(hist.axis.Integer(0, 4, name="electron_nearGenA_n", label="$N_{e}$"),
                   lambda objs, mask: ak.num(objs["electrons"][dR(objs["electrons"],
                                                              objs["genAs_toE"]) < 0.5])),
        ],
    ),
    
    # pfelectron-genA
    "electron_nearGenA_n_genA_lxy": h.Histogram(
        [
            # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])),
            # number of electrons within dR=0.5 of a genA that decays to electrons
            h.Axis(hist.axis.Integer(0, 4, name="electron_nearGenA_n", label="$N_{e}$"),
                   lambda objs, mask: ak.num(matched(objs["electrons"], objs["genAs_toE"], 0.5))),
        ],
    ),
    # pfelectron-genElectron
    "electron_genE_dR": h.Histogram(
        [
            # dR(e, nearest gen e)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="electron_genE_dR"),
                   lambda objs, mask: dR(objs["electrons"], objs["genEs"]))
        ],
    ),
    # pfphoton
    "photon_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="photon_n"),
                   lambda objs, mask: ak.num(objs["photons"])),
        ],
    ),
    "photon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="photon_pt"),
                   lambda objs, mask: objs["photons"].pt),
        ],
    ),
    "photon_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="photon_eta"),
                   lambda objs, mask: objs["photons"].eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="photon_phi"),
                   lambda objs, mask: objs["photons"].phi),
        ],
    ),
    "photon_nearGenA_n": h.Histogram(
        [
            # number of photons within dR=0.5 of a genA that decays to electrons
            h.Axis(hist.axis.Integer(0, 10, name="photon_nearGenA_n"),
                   lambda objs, mask: ak.num(matched(objs["photons"], objs["genAs_toE"], 0.5))),
        ],
    ),
    # pfphoton-genA
    "photon_nearGenA_n_genA_lxy": h.Histogram(
        [
            # lxy of dark photon that decays to electrons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toE"])),
            # number of photons within dR=0.5 of a genA that decays to electrons
            h.Axis(hist.axis.Integer(0, 4, name="photon_nearGenA_n", label="$N_{\gamma}$"),
                   lambda objs, mask: ak.num(matched(objs["photons"], objs["genAs_toE"], 0.5))),
        ],
    ),
    # pfphoton-genElectron
    "photon_genE_dR": h.Histogram(
        [
            # dR(photon, nearest gen e)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="photon_genE_dR"),
                   lambda objs, mask: dR(objs["photons"], objs["genEs"]))
        ],
    ),
    # pfmuon
    "muon_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="muon_n"),
                   lambda objs, mask: ak.num(objs["muons"])),
        ],
    ),
    "muon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="muon_pt"),
                   lambda objs, mask: objs["muons"].pt),
        ],
    ),
    "muon_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="muon_eta"),
                   lambda objs, mask: objs["muons"].eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="muon_phi"),
                   lambda objs, mask: objs["muons"].phi),
        ],
    ),
    "muon_absD0": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="muon_absD0", label=r"Muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["muons"].d0)),
        ],
    ),
    "muon_absD0_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 10, name="muon_absD0_lowRange",
                                     label=r"Muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["muons"].d0)),
        ],
    ),
    "muon_nearGenA_n": h.Histogram(
        [
            # number of muons within dR=0.5 of a genA that decays to muons
            h.Axis(hist.axis.Integer(0, 10, name="muon_nearGenA_n"),
                   lambda objs, mask: ak.num(matched(objs["muons"], objs["genAs_toMu"], 0.5))),
        ],
    ),
    # pfmuon-genA
    "muon_nearGenA_n_genA_lxy": h.Histogram(
        [
            # lxy of dark photon that decays to muons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toMu"])),
            # number of muons within dR=0.5 of a genA that decays to muons
            h.Axis(hist.axis.Integer(0, 4, name="muon_nearGenA_n", label="$N_{\mu^{PF}}$"),
                   lambda objs, mask: ak.num(matched(objs["muons"], objs["genAs_toMu"], 0.5))),
        ],
    ),
    # pfmuon-genMuon
    "muon_genMu_dR": h.Histogram(
        [
            # dR(mu, nearest gen mu)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="muon_genMu_dR"),
                   lambda objs, mask: dR(objs["muons"], objs["genMus"]))
        ],
    ),
    # dsamuon
    "dsaMuon_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="dsaMuon_n"),
                   lambda objs, mask: ak.num(objs["dsaMuons"])),
        ],
    ),
    "dsaMuon_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="dsaMuon_pt"),
                   lambda objs, mask: objs["dsaMuons"].pt),
        ],
    ),
    "dsaMuon_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="dsaMuon_eta"),
                   lambda objs, mask: objs["dsaMuons"].eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="dsaMuon_phi"),
                   lambda objs, mask: objs["dsaMuons"].phi),
        ],
    ),
    "dsaMuon_absD0": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="dsaMuon_absD0",
                                     label=r"DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"].d0)),
        ],
    ),
    "dsaMuon_absD0_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 10, name="dsaMuon_absD0_lowRange",
                                     label=r"DSA muon $|d_0|$ [cm]"),
                   lambda objs, mask: abs(objs["dsaMuons"].d0)),
        ],
    ),
    "dsaMuon_nearGenA_n": h.Histogram(
        [
            # number of muons within dR=0.5 of a genA that decays to muons
            h.Axis(hist.axis.Integer(0, 10, name="dsaMuon_nearGenA_n"),
                   lambda objs, mask: ak.num(matched(objs["dsaMuons"], objs["genAs_toMu"], 0.5))),
        ],
    ),
    # dsamuon-genA
    "dsaMuon_nearGenA_n_genA_lxy": h.Histogram(
        [
            # lxy of dark photon that decays to dsaMuons
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs_toMu"])),
            # number of dsaMuons within dR=0.5 of a genA that decays to muons
            h.Axis(hist.axis.Integer(0, 4, name="dsaMuon_nearGenA_n", label="$N_{\mu^{DSA}}$"),
                   lambda objs, mask: ak.num(matched(objs["dsaMuons"], objs["genAs_toMu"], 0.5))),
        ],
    ),
    # dsaMuon-genMuon
    "dsaMuon_genMu_dR": h.Histogram(
        [
            # dR(dsa mu, nearest gen mu)
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="dsaMuon_genMu_dR"),
                   lambda objs, mask: dR(objs["dsaMuons"], objs["genMus"]))
        ],
    ),
    # lj
    "lj_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_n"),
                   lambda objs, mask: ak.num(objs["ljs"])),
        ],
    ),
    "lj_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 100, name="lj_pt", label="Lepton jet pT [GeV]"),
                   lambda objs, mask: objs["ljs"].pt),
        ],
    ),
    "lj_pfIsolation05": h.Histogram(
        [
            h.Axis(hist.axis.Regular(80, 0, 0.8, name="lj_pfIsolation05",
                                     label="Lepton jet isolation"),
                   lambda objs, mask: objs["ljs"].pfIsolation05),
        ],
    ),
    "lj0_pfIsolation05": h.Histogram(
        [
            h.Axis(hist.axis.Regular(80, 0, 0.8, name="lj_pfIsolation05",
                                     label="Leading lepton jet isolation"),
                   lambda objs, mask: objs["ljs"][mask, 0].pfIsolation05),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 0,
    ),
    "lj1_pfIsolation05": h.Histogram(
        [
            h.Axis(hist.axis.Regular(80, 0, 0.8, name="lj_pfIsolation05",
                                     label="Subleading lepton jet isolation"),
                   lambda objs, mask: objs["ljs"][mask, 1].pfIsolation05),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_pfIsolationPtNoPU05": h.Histogram( # not in v2 ntuples
        [
            h.Axis(hist.axis.Regular(80, 0, 0.8, name="lj_pfIsolationPtNoPU05",
                                     label="Lepton jet isolation"),
                   lambda objs, mask: objs["ljs"].pfIsolationPtNoPU05),
        ],
    ),
    "lj_pfIsolationPt05": h.Histogram( # not in v2 ntuples
        [
            h.Axis(hist.axis.Regular(80, 0, 0.8, name="lj_pfIsolationPt05",
                                     label="Lepton jet isolation"),
                   lambda objs, mask: objs["ljs"].pfIsolationPt05),
        ],
    ),
    "lj_pfIsolation07": h.Histogram(
        [
            h.Axis(hist.axis.Regular(80, 0, 0.8, name="lj_pfIsolation07",
                                     label="Lepton jet isolation"),
                   lambda objs, mask: objs["ljs"].pfIsolation07),
        ],
    ),
    "lj_pfIsolationPtNoPU07": h.Histogram( # not in v2 ntuples
        [
            h.Axis(hist.axis.Regular(80, 0, 0.8, name="lj_pfIsolationPtNoPU07",
                                     label="Lepton jet isolation"),
                   lambda objs, mask: objs["ljs"].pfIsolationPtNoPU07),
        ],
    ),
    "lj_pfIsolationPt07": h.Histogram( # not in v2 ntuples
        [
            h.Axis(hist.axis.Regular(80, 0, 0.8, name="lj_pfIsolationPt07",
                                     label="Lepton jet isolation"),
                   lambda objs, mask: objs["ljs"].pfIsolationPt07),
        ],
    ),
    "lj_pfiso": h.Histogram(
        [
            h.Axis(hist.axis.Regular(80, 0, 0.8, name="lj_pfiso",
                                     label="Lepton jet isolation"),
                   lambda objs, mask: objs["ljs"].pfiso),
        ],
    ),
    "lj0_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 100, name="lj0_pt",
                                     label="Leading lepton jet pT [GeV]"),
                   lambda objs, mask: objs["ljs"][mask, 0].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 0,
    ),
    "lj1_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 100, name="lj1_pt",
                                     label="Subleading lepton jet pT [GeV]"),
                   lambda objs, mask: objs["ljs"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj0_e": h.Histogram(
        [
            h.Axis(hist.axis.Regular(350, 0, 700, name="lj_e",
                                     label="Leading lepton jet E [GeV]"),
                   lambda objs, mask: objs["ljs"][mask, 0].energy),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 0,
    ),
    "lj1_e": h.Histogram(
        [
            h.Axis(hist.axis.Regular(350, 0, 700, name="lj_e",
                                     label="Subleading lepton jet E [GeV]"),
                   lambda objs, mask: objs["ljs"][mask, 1].energy),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj0_dRSpread": h.Histogram(
        [
            h.Axis(hist.axis.Regular(250, 0, 1.0, name="lj0_dRSpread",
                                     label="Leading lepton jet dRSpread"),
                   lambda objs, mask: objs["ljs"][mask, 0].dRSpread),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 0,
    ),
    "lj1_dRSpread": h.Histogram(
        [
            h.Axis(hist.axis.Regular(250, 0, 1.0, name="lj1_dRSpread",
                                     label="Subleading lepton jet dRSpread"),
                   lambda objs, mask: objs["ljs"][mask, 1].dRSpread),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="lj_eta"),
                   lambda objs, mask: objs["ljs"].eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="lj_phi"),
                   lambda objs, mask: objs["ljs"].phi),
        ],
    ),
    "egm_lj_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 100, name="egm_lj_pt",
                                     label="EGM-type lepton jet pT [GeV]"),
                   lambda objs, mask: derived_objs["egm_ljs"](objs).pt),
        ],
    ),
    "mu_lj_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 100, name="mu_lj_pt",
                                     label="Mu-type lepton jet pT [GeV]"),
                   lambda objs, mask: derived_objs["mu_ljs"](objs).pt),
        ],
    ),
    "lj_electronN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_electronN"),
                   lambda objs, mask: objs["ljs"].electron_n),
        ],
    ),
    "lj_photonN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_photonN"),
                   lambda objs, mask: objs["ljs"].photon_n),
        ],
    ),
    "lj_electronPhotonN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_electronPhotonN"),
                   lambda objs, mask: objs["ljs"].electron_n + objs["ljs"].photon_n),
        ],
    ),
    "lj_muonN": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="lj_muonN"),
                   lambda objs, mask: objs["ljs"].muon_n),
        ],
    ),
    # ljsource
    "ljsource_n": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 10, name="ljsource_n"),
                   lambda objs, mask: ak.num(objs["ljsources"])),
        ],
    ),
    "ljsource_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 100, name="ljsource_pt",
                                     label="Lepton jet source pT [GeV]"),
                   lambda objs, mask: objs["ljsources"].pt),
        ],
    ),
    "ljsource_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="ljsource_eta"),
                   lambda objs, mask: objs["ljsources"].eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="ljsource_phi"),
                   lambda objs, mask: objs["ljsources"].phi),
        ],
    ),
    "ljsource_charge": h.Histogram(
        [
            h.Axis(hist.axis.Integer(-1, 1, name="ljsource_charge"),
                   lambda objs, mask: objs["ljsources"].charge),
        ],
    ),
    "ljsource_type": h.Histogram(
        [
            h.Axis(hist.axis.IntCategory([2, 3, 4, 8], name="lj_type"),
                   lambda objs, mask: objs["ljsources"]["type"]), # avoid ak.Array.type
        ],
    ),
    # pfelectron-lj
    "electron_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="electron_lj_dR"),
                   lambda objs, mask: dR(objs["electrons"], objs["ljs"]))
        ],
    ),
    "electron_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 1.0, name="electron_lj_dR_lowRange"),
                   lambda objs, mask: dR(objs["electrons"], objs["ljs"]))
        ],
    ),
    # pfphoton-lj
    "photon_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="photon_lj_dR"),
                   lambda objs, mask: dR(objs["photons"], objs["ljs"]))
        ],
    ),
    "photon_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 1.0, name="photon_lj_dR_lowRange"),
                   lambda objs, mask: dR(objs["photons"], objs["ljs"]))
        ],
    ),
    "photon_lj_dR_reallyLowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 0.1, name="photon_lj_dR_reallyLowRange"),
                   lambda objs, mask: dR(objs["photons"], objs["ljs"]))
        ],
    ),
    # pfmuon-lj
    "muon_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="muon_lj_dR"),
                   lambda objs, mask: dR(objs["muons"], objs["ljs"]))
        ],
    ),
    "muon_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 1.0, name="muon_lj_dR_lowRange"),
                   lambda objs, mask: dR(objs["muons"], objs["ljs"]))
        ],
    ),
    # dsamuon-lj
    "dsaMuon_lj_dR": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(100, 0, 2*math.pi, name="dsaMuon_lj_dR"),
                   lambda objs, mask: dR(objs["dsaMuons"], objs["ljs"]))
        ],
    ),
    "dsaMuon_lj_dR_lowRange": h.Histogram(
        [
            # dR(e, nearest LJ)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="dsaMuon_lj_dR_lowRange"),
                   lambda objs, mask: dR(objs["dsaMuons"], objs["ljs"]))
        ],
    ),
    # lj-lj
    "lj_lj_absdphi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, 0, 2*math.pi, name="ljlj_absdphi"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].phi
                                          - objs["ljs"][mask, 0].phi)),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_invmass": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2000, name="ljlj_mass", label=r"InvMass($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, :2].sum().mass),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    "lj_lj_invmass_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="ljlj_mass", label=r"InvMass($LJ_{0}$, $LJ_{1}$)"),
                   lambda objs, mask: objs["ljs"][mask, :2].sum().mass),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    # ABCD plane
    "abcd_lj_lj_dphi_vs_lj0_pfIsolationPt05": h.Histogram( # not in v2 ntuples
        [
            h.Axis(hist.axis.Regular(200, 0, 2*math.pi, name="ljlj_absdphi",
                                     label=r"Lepton jet pair |$\Delta\phi$|"),
                   lambda objs, mask: abs(objs["ljs"][mask, 1].phi
                                          - objs["ljs"][mask, 0].phi)),
            h.Axis(hist.axis.Regular(80, 0, 0.8, name="lj_pfIsolationPt05",
                                     label="Leading lepton jet isolation"),
                   lambda objs, mask: objs["ljs"][mask, 0].pfIsolationPt05),
        ],
        evt_mask=lambda objs: ak.num(objs["ljs"]) > 1,
    ),
    # gen
    "gen_abspid": h.Histogram(
        [
            h.Axis(hist.axis.Integer(0, 40, name="gen_abspid"),
                   lambda objs, mask: abs(objs["gens"].pid)),
        ],
    ),
    # genelectron
    "genE_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="genE_n"),
                   lambda objs, mask: ak.num(objs["genEs"])),
        ],
    ),
    "genE_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE_pt",
                                     label=r"Gen-level electron $p_{T}$ [GeV]"),
                   lambda objs, mask: abs(objs["genEs"].pt)),
        ],
    ),
    "genE_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(70, 0, 700, name="genE_pt",
                                     label=r"Gen-level electron $p_{T}$ [GeV]"),
                   lambda objs, mask: abs(objs["genEs"].pt)),
        ],
    ),
    "genE0_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE0_pt",
                                     label=r"Leading gen-level electron $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genEs"][mask, 0].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 0,
    ),
    "genE0_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(70, 0, 700, name="genE_pt",
                                     label=r"Leading gen-level electron $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genEs"][mask, 0].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 0,
    ),
    "genE1_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE1_pt",
                                     label=r"Sub-leading gen-level electron $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genEs"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE1_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(70, 0, 700, name="genE_pt",
                                     label=r"Sub-leading gen-level electron $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genEs"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="genE_eta", label=r"Gen-level electron $\eta$"),
                   lambda objs, mask: objs["genEs"].eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="genE_phi",
                                     label=r"Gen-level electron \phi"),
                   lambda objs, mask: objs["genEs"].phi),
        ],
    ),
    # genelectron-genelectron
    "genE_genE_dR": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(200, 0, 1.0, name="genE_genE_dR",
                                     label=r"$\Delta R$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE_genE_dR_lowRange": h.Histogram(
        [
            # dR(subleading gen E, leading gen E)
            h.Axis(hist.axis.Regular(100, 0, 0.5, name="genE_genE_dR_lowRange",
                                     label=r"$\Delta R$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: objs["genEs"][mask, 1].delta_r(objs["genEs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    "genE_genE_dEta": h.Histogram(
        [
            # abs(dEta(subleading gen E, leading gen E))
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genE_genE_dEta",
                                     label=r"$\Delta\, \eta$($e_0^{gen}$, $e_1^{gen}$)"),
                   lambda objs, mask: abs(objs["genEs"][mask, 1].eta
                                          - objs["genEs"][mask, 0].eta)),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genE_genE_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genE_genE_pt"),
                   lambda objs, mask: objs["genEs"][mask, :2].sum().pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genEs"]) > 1,
    ),
    # genmuon
    "genMu_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="genMu_n"),
                   lambda objs, mask: ak.num(objs["genMus"])),
        ],
    ),
    "genMu_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 200, name="genMu_pt",
                                     label=r"Gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: abs(objs["genMus"].pt)),
        ],
    ),
    "genMu_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 700, name="genMu_pt",
                                     label=r"Gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: abs(objs["genMus"].pt)),
        ],
    ),
    "genMu0_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 200, name="genMu0_pt",
                                     label=r"Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 0].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 0
    ),
    "genMu0_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(70, 0, 700, name="genMu0_pt",
                                     label=r"Leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 0].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 0,
    ),
    "genMu1_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genMu1_pt",
                                     label=r"Sub-leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu1_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(70, 0, 700, name="genMu1_pt",
                                     label=r"Sub-leading gen-level muon $p_{T}$ [GeV]"),
                   lambda objs, mask: objs["genMus"][mask, 1].pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="genMu_eta", label=r"Gen-level muon $\eta$"),
                   lambda objs, mask: objs["genMus"].eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="genMu_phi",
                                     label=r"Gen-level muon \phi"),
                   lambda objs, mask: objs["genMus"].phi),
        ],
    ),
    # genmuon-genmuon
    "genMu_genMu_dR": h.Histogram(
        [
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genMu_genMu_dR",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(
                       objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_genMu_dR_lowRange": h.Histogram(
        [
            # dR(subleading gen Mu, leading gen Mu)
            h.Axis(hist.axis.Regular(100, 0, 0.5, name="genMu_genMu_dR_lowRange",
                                     label=r"$\Delta R$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: objs["genMus"][mask, 1].delta_r(
                       objs["genMus"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_genMu_dEta": h.Histogram(
        [
            # abs(dEta(subleading gen Mu, leading gen Mu))
            h.Axis(hist.axis.Regular(50, 0, 1.0, name="genMu_genMu_dEta",
                                     label=r"$\Delta\, \eta$($\mu_0^{gen}$, $\mu_1^{gen}$)"),
                   lambda objs, mask: abs(objs["genMus"][mask, 1].eta
                                          - objs["genMus"][mask, 0].eta)),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    "genMu_genMu_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genMu_genMu_pt"),
                   lambda objs, mask: objs["genMus"][mask, :2].sum().pt),
        ],
        evt_mask=lambda objs: ak.num(objs["genMus"]) > 1,
    ),
    # dsamuon-genmuon, dR 0.4 window
    "dsaMuon_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="dsaMuon_genMu_ptRatio"),
                   lambda objs, mask: objs["dsaMuons"].pt
                       / objs["dsaMuons"].nearest(objs["genMus"], threshold=0.4).pt),
        ],
    ),
    "dsaMuon0_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="dsaMuon0_genMu_ptRatio"),
                   lambda objs, mask: (objs["dsaMuons"][mask, 0:1].pt
                       / objs["dsaMuons"][mask, 0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ak.num(matched(objs["dsaMuons"] ,objs["genMus"], 0.4)) > 0,
    ),
    "dsaMuon1_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="dsaMuon1_genMu_ptRatio"),
                   lambda objs, mask: (objs["dsaMuons"][mask, 1:2].pt
                       / objs["dsaMuons"][mask, 1:2].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ak.num(matched(objs["dsaMuons"], objs["genMus"], 0.4)) > 1,
    ),
    # pfmuon-genmuon, dR 0.4 window
    "pfMuon_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="pfMuon_genMu_ptRatio"),
                   lambda objs, mask: objs["muons"].pt
                       / objs["muons"].nearest(objs["genMus"], threshold=0.4).pt),
        ],
    ),
    "pfMuon0_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="pfMuon0_genMu_ptRatio"),
                   lambda objs, mask: (objs["muons"][mask,0:1].pt
                       / objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ak.num(matched(objs["muons"], objs["genMus"], 0.4)) > 0,
    ),
    "pfMuon1_genMu_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="pfMuon1_genMu_ptRatio"),
                   lambda objs, mask: (objs["muons"][mask,1:2].pt
                       / objs["muons"][mask,1:2].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ak.num(matched(objs["muons"], objs["genMus"], 0.4)) > 1,
    ),
    # gen dark photons (A)
    "genA_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="genA_n"),
                   lambda objs, mask: ak.num(objs["genAs"])),
        ],
    ),
    "genA_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genA_pt",
                                     label=r"Dark photon $p_{T}$ [GeV]"),
                   lambda objs, mask: abs(objs["genAs"].pt)),
        ],
    ),
    "genA_pt_highRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(140, 0, 700, name="genA_pt",
                                     label=r"Dark photon $p_{T}$ [GeV]"),
                   lambda objs, mask: abs(objs["genAs"].pt)),
        ],
    ),
    "genA_eta_phi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(50, -3, 3, name="genA_eta", label=r"Dark photon $\eta$"),
                   lambda objs, mask: objs["genAs"].eta),
            h.Axis(hist.axis.Regular(50, -1*math.pi, math.pi, name="genA_phi",
                                     label=r"Dark photon \phi"),
                   lambda objs, mask: objs["genAs"].phi),
        ],
    ),
    "genA_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs"]) ),
        ],
    ),
    "genA_lxy_lowRange": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 10, name="genA_lxy",
                                     label=r"Dark photon $L_{xy}$ [cm]"),
                   lambda objs, mask: lxy(objs["genAs"]) ),
        ],
    ),
    "genAs_toMu_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="genAs_toMu_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toMu"]) ),
        ],
    ),
    "genAs_toMu_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 700, name="genAs_toMu_pt",
                                    label=r"Dark photon (matched to Muon) $p_{T}$ [GeV]"),
                   lambda objs, mask: abs(objs["genAs_toMu"].pt) ),
        ],
    ),
    "genAs_toE_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="genAs_toE_lxy"),
                   lambda objs, mask: lxy(objs["genAs_toE"]) ),
        ],
    ),
    "genAs_toE_pt": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 700, name="genAs_toE_pt",
                                    label=r"Dark photon (matched to EGM) $p_{T}$ [GeV]"),
                   lambda objs, mask: abs(objs["genAs_toE"].pt) ),
        ],
    ),
    "genAs_toMu_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="genAs_toMu_n"),
                   lambda objs, mask: ak.num(objs["genAs_toMu"]) ),
        ],
    ),
    "genA_pt_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 200, name="genA_pt"),
                   lambda objs, mask: abs(objs["genAs"].pt)),
            h.Axis(hist.axis.Regular(250, 0, 500, name="genA_lxy"),
                   lambda objs, mask: lxy(objs["genAs"])),
        ],
    ),
    # genA-genA
    "genA_genA_dphi": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, math.pi, name="genA_genA_dphi",
                                     label=r"$\Delta\phi$ between dark photons"),
                   lambda objs, mask: objs["genAs"][mask, 1].delta_phi(objs["genAs"][mask, 0])),
        ],
        evt_mask=lambda objs: ak.num(objs["genAs"]) > 1,
    ),
    # genA-LJ
    "genA_lj_dR": h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(200, 0, 2*math.pi, name="genA_lj_dR"),
                   lambda objs, mask: dR(objs["genAs"], objs["ljs"]))
        ],
    ),
    "genA_lj_dR_lowRange": h.Histogram(
        [
            # dR(A, nearest LJ)
            h.Axis(hist.axis.Regular(200, 0, 1.0, name="genA_lj_dR_lowRange"),
                   lambda objs, mask: dR(objs["genAs"], objs["ljs"]))
        ],
    ),
    "genA_matched_lj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="matched_genA_lxy"),
                   lambda objs, mask: lxy(derived_objs["genAs_matched_lj"](objs, 0.4)) ),
        ],
    ),
    "genA_toMu_matched_lj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_toMu_matched_lj_lxy"),
                   lambda objs, mask: lxy(derived_objs["genAs_toMu_matched_lj"](objs, 0.4)) ),
        ],
    ),
    "genA_toE_matched_lj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_toE_matched_lj_lxy"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_lj"](objs, 0.4)) ),
        ],
    ),
    "genA_matched_muLj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_matched_muLj_lxy"),
                   lambda objs, mask: lxy(derived_objs["genAs_matched_muLj"](objs, 0.4)) ),
        ],
    ),
    "genA_toMu_matched_muLj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_toMu_matched_muLj_lxy"),
                   lambda objs, mask: lxy(derived_objs["genAs_toMu_matched_muLj"](objs, 0.4)) ),
        ],
    ),
    "genA_matched_egmLj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_matched_egmLj_lxy"),
                   lambda objs, mask: lxy(derived_objs["genAs_matched_egmLj"](objs, 0.4)) ),
        ],
    ),
    "genA_toE_matched_egmLj_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 500, name="genA_toE_matched_egmLj_lxy"),
                   lambda objs, mask: lxy(derived_objs["genAs_toE_matched_egmLj"](objs, 0.4)) ),
        ],
    ),
    "genA_matched_lj_n": h.Histogram(
        [
            h.Axis(hist.axis.Regular(10, 0, 10, name="genA_matched_n"),
                   lambda objs, mask: ak.num(derived_objs["genAs_matched_lj"](objs, 0.4)) ),
        ],
    ),
    # genA - LJ 0.4 matching radius, pT Ratios
    "genA_lj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_lj_ptRatio",
                   label="Lepton Jet pT / (closest) Dark Photon pT"),
                   lambda objs, mask: objs["ljs"].pt
                       / objs["ljs"].nearest(objs["genAs"], threshold=0.4).pt),
        ],
    ),
    "genA_egmLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_egmLj_ptRatio",
                   label="EGM Lepton Jet pT / (closest) Dark Photon pT"),
                   lambda objs, mask: derived_objs["egm_ljs"](objs).pt
                       / derived_objs["egm_ljs"](objs).nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_oneElectronLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_oneElectronLj_ptRatio",
                   label="(1) Electron Lepton Jet / (closest) Dark Photon pT"),
                   lambda objs, mask: derived_objs["electron_ljs"](objs,1).pt
                       / derived_objs["electron_ljs"](objs,1).nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_twoElectronLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_twoElectronLj_ptRatio",
                   label="(2) Electron Lepton Jet / (closest) Dark Photon pT"),
                   lambda objs, mask: derived_objs["electron_ljs"](objs,2).pt
                       / derived_objs["electron_ljs"](objs,2).nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_onePhotonLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_onePhotonLj_ptRatio",
                   label="(1) Photon Lepton Jet / (closest) Dark Photon pT"),
                   lambda objs, mask: derived_objs["photon_ljs"](objs,1).pt
                       / derived_objs["photon_ljs"](objs,1).nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_twoPhotonLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_twoPhotonLj_ptRatio",
                   label="(2) Photon Lepton Jet / (closest) Dark Photon pT"),
                   lambda objs, mask: derived_objs["photon_ljs"](objs,2).pt
                       / derived_objs["photon_ljs"](objs,2).nearest(objs["genAs_toE"], threshold=0.4).pt),
        ],
    ),
    "genA_muLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_muLj_ptRatio",
                   label="Muon Lepton Jet pT / (closest) Dark Photon pT"),
                   lambda objs, mask: derived_objs["mu_ljs"](objs).pt
                       / derived_objs["mu_ljs"](objs).nearest(objs["genAs_toMu"], threshold=0.4).pt),
        ],
    ),
    "genA_dsaMuonLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_dsaMuonLj_ptRatio",
                   label="DSA Muon Lepton Jet pT / (closest) Dark Photon pT"),
                   lambda objs, mask: (objs["dsaMuons"][mask]).nearest(objs["ljs"][mask], threshold=0.4).pt
                       / (objs["dsaMuons"][mask]).nearest(objs["genAs_toMu"][mask], threshold=0.4).pt),
        ],
    ),
    "genA_pfMuonLj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_pfMuonLj_ptRatio",
                   label="PF Muon Lepton Jet pT / (closest) Dark Photon pT"),
                   lambda objs, mask: (objs["muons"][mask]).nearest(objs["ljs"][mask], threshold=0.4).pt
                       / (objs["muons"][mask]).nearest(objs["genAs_toMu"][mask], threshold=0.4).pt),
        ],
    ),
    "genA_dsaMuon0Lj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_dsaMuonLj_ptRatio",
                   label="Lead DSA Muon Lepton Jet / (closest) Dark Photon pT"),
                   lambda objs, mask: ((objs["dsaMuons"][mask, 0:1]).nearest(objs["ljs"][mask], threshold=0.4).pt
                       / (objs["dsaMuons"][mask, 0:1]).nearest(objs["genAs_toMu"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["dsaMuons"], objs["genAs_toMu"], 0.4)) > 0)
                               & (ak.num(matched(objs["dsaMuons"], objs["ljs"], 0.4)) > 0)),
    ),
    "genA_pfMuon0Lj_ptRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_pfMuonLj_ptRatio",
                   label="Lead PF Muon Lepton Jet / (closest) Dark Photon pT"),
                   lambda objs, mask: ((objs["muons"][mask, 0:1]).nearest(objs["ljs"][mask], threshold=0.4).pt
                       / (objs["muons"][mask, 0:1]).nearest(objs["genAs_toMu"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["muons"], objs["genAs_toMu"], 0.4)) > 0)
                               & (ak.num(matched(objs["muons"], objs["ljs"], 0.4)) > 0)),
    ),
    # genA - LJ 0.4 matching radius, LJ Reco Lxy / True Lxy
    "genA_muLj_lxyRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_muLj_lxyRatio",
                                    label="Muon Lepton Jet Reco L$_{xy}$ / (closest) Dark Photon L$_{xy}$"),
                   lambda objs, mask: derived_objs["mu_ljs"](objs).kinvtx.lxy
                       / lxy(derived_objs["mu_ljs"](objs).nearest(objs["genAs"], threshold=0.4))),
        ],
    ),
    "genA_egmLj_lxyRatio": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="genA_egmLj_lxyRatio",
                                    label="EGM Lepton Jet Reco L$_{xy}$ / (closest) Dark Photon L$_{xy}$"),
                   lambda objs, mask: derived_objs["egm_ljs"](objs).kinvtx.lxy
                       / lxy(derived_objs["egm_ljs"](objs).nearest(objs["genAs"], threshold=0.4))),
        ],
    ),
    # LJ Res vs Reco Lxy
    "mu_lj_genA_ptRatio_vs_recolxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="mu_lj_genA_ptRatio"),
                   lambda objs, mask: derived_objs["mu_ljs"](objs).pt
                       / derived_objs["mu_ljs"](objs).nearest(objs["genAs"]).pt),
            h.Axis(hist.axis.Regular(100, 0, 300, name="mu_lj_recolxy"),
                   lambda objs, mask: derived_objs["mu_ljs"](objs).kinvtx.lxy),
        ],
    ),
    "egm_lj_genA_ptRatio_vs_recolxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="egm_lj_genA_ptRatio"),
                   lambda objs, mask: derived_objs["egm_ljs"](objs).pt
                       / derived_objs["egm_ljs"](objs).nearest(objs["genAs"]).pt),
            h.Axis(hist.axis.Regular(100, 0, 300, name="egm_lj_recolxy"),
                   lambda objs, mask: derived_objs["egm_ljs"](objs).kinvtx.lxy),
        ],
    ),
    # LJ Res vs True Lxy, 0.4 thresholds on dR matching
    "egm_lj_genA_ptRatio_vs_truelxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="egm_lj_genA_ptRatio"),
                   lambda objs, mask: derived_objs["egm_ljs"](objs).pt
                       / derived_objs["egm_ljs"](objs).nearest(objs["genAs"], threshold=0.4).pt),
            h.Axis(hist.axis.Regular(100, 0, 300, name="egm_lj_truelxy"),
                   lambda objs, mask: lxy(derived_objs["egm_ljs"](objs).nearest(objs["genAs"], threshold=0.4))),
        ],
    ),
    "mu_lj_genA_ptRatio_vs_truelxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 2.0, name="mu_lj_genA_ptRatio"),
                   lambda objs, mask: derived_objs["mu_ljs"](objs).pt
                       / derived_objs["mu_ljs"](objs).nearest(objs["genAs"], threshold=0.4).pt),
            h.Axis(hist.axis.Regular(100, 0, 300, name="mu_lj_truelxy"),
                   lambda objs, mask: lxy(derived_objs["mu_ljs"](objs).nearest(objs["genAs"], threshold=0.4))),
        ],
    ),
    "dsaMuon0_genMu0_ptRatio_vs_truelxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0., 2.0, name="dsaMuon0_genMu0_ptRatio"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].pt
                       / objs["dsaMuons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(100, 0, 300, name="dsaMuon0_lj_truelxy"),
                   lambda objs, mask: lxy(objs["dsaMuons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4))),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["genMus"], objs["dsaMuons"], 0.4)) > 0)
                               & (ak.num(matched(objs["genAs"], objs["dsaMuons"], 0.4)) > 0)),
    ),
    "muon0_genMu0_ptRatio_vs_truelxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0., 2.0, name="muon0_genMu0_ptRatio"),
                   lambda objs, mask: (objs["muons"][mask,0:1].pt
                       / objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(100, 0, 300, name="pfMuon0_lj_truelxy"),
                   lambda objs, mask: lxy(objs["muons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4))),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["genMus"], objs["muons"], 0.4)) > 0)
                               & (ak.num(matched(objs["genAs"], objs["muons"], 0.4)) > 0)),
    ),
    # LJ Res vs True pT, dR 0.4 matching window
    "dsaMuon0_genMu0_ptRatio_vs_truept": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2.0, name="dsaMuon0_genMu0_ptRatio"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].pt
                       / objs["dsaMuons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: (ak.num(objs["dsaMuons"]) > 0)
    ),
    "muon0_genMu0_ptRatio_vs_truept": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2.0, name="muon0_genMu0_ptRatio"),
                   lambda objs, mask: (objs["muons"][mask,0:1].pt
                       / objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: (ak.num(objs["muons"]) > 0)
    ),
    "dsaMuon0_muLj_ptRatio_vs_truept": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0, 2.0, name="dsaMuon0_genMu0_ptRatio"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].nearest(objs["ljs"][mask], threshold=0.4).pt
                       / objs["dsaMuons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: (ak.num(objs["dsaMuons"]) > 0)
    ),
    "muon0_muLj_ptRatio_vs_truept": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0., 2.0, name="muon0_muLj_ptRatio"),
                   lambda objs, mask: (objs["muons"][mask,0:1].nearest(objs["ljs"][mask], threshold=0.4).pt
                       / objs["muons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["ljs"], objs["muons"], 0.4)) > 0)
                               & (ak.num(matched(objs["genMus"], objs["muons"], 0.4)) > 0)),
    ),
    "egmLj_ptRatio_vs_egm_truept": h.Histogram(
        [
            h.Axis(hist.axis.Regular(100, 0., 2.0, name="egm_lj_genA_ptRatio"),
                   lambda objs, mask: derived_objs["egm_ljs"](objs).pt
                       / derived_objs["egm_ljs"](objs).nearest(objs["genAs"], threshold=0.4).pt),
            h.Axis(hist.axis.Regular(100, 0, 1000, name="genE_pt"),
                   lambda objs, mask: (derived_objs["egm_ljs"](objs).nearest(objs["genEs"], threshold=0.4).pt)[mask,0]),
        ],
    ),
    # LJ True pT vs True Lxy, dR 0.4 matching window
    "genMu0_truept_vs_dsaMuon0_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["dsaMuons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(100, 0, 300, name="dsaMuon0_lj_truelxy"),
                   lambda objs, mask: lxy(objs["dsaMuons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4))),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["genMus"], objs["dsaMuons"], 0.4)) > 0)
                               & (ak.num(matched(objs["genAs"], objs["dsaMuons"], 0.4)) > 0)),
    ),
    "genMu0_truept_vs_muon0_lxy": h.Histogram(
        [
            h.Axis(hist.axis.Regular(200, 0, 1000, name="genMu0_pt"),
                   lambda objs, mask: (objs["muons"][mask,0:1].nearest(objs["genMus"][mask], threshold=0.4).pt)),
            h.Axis(hist.axis.Regular(100, 0, 300, name="pfMuon0_lj_truelxy"),
                   lambda objs, mask: lxy(objs["muons"][mask,0:1].nearest(objs["genAs"][mask], threshold=0.4))),
        ],
        evt_mask=lambda objs: ((ak.num(matched(objs["genMus"], objs["muons"], 0.4)) > 0)
                               & (ak.num(matched(objs["genAs"], objs["muons"], 0.4)) > 0)),
    ),
}
