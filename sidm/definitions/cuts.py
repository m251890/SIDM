"""Define all available cuts"""

# columnar analysis
import awkward as ak
# local
from sidm.definitions.objects import derived_objs
from sidm.tools.utilities import dR, lxy, check_bit


obj_cut_defs = {
    "pvs": {
        "ndof > 4": lambda objs: objs["pvs"].ndof > 4,
        "|z| < 24 cm": lambda objs: abs(objs["pvs"].z) < 24,
        "|rho| < 0.2 mm": lambda objs: abs(objs["pvs"].rho) < 0.2,
    },
    "ljs": {
        "pT > 30 GeV": lambda objs: objs["ljs"].pt > 30,
        "|eta| < 2.4": lambda objs: abs(objs["ljs"].eta) < 2.4,
        "dR(LJ, A) < 0.2": lambda objs: dR(objs["ljs"], objs["genAs"]) < 0.2,
    },
    "genAs": {
        "dR(A, LJ) < 0.2": lambda objs: dR(objs["genAs"], objs["ljs"]) < 0.2,
        "lxy < 10 cm": lambda objs: lxy(objs["genAs"]) < 10,
        "10 cm <= lxy < 100 cm": lambda objs: ((lxy(objs["genAs"]) >= 10)
                                               & (lxy(objs["genAs"]) < 100)),
        "lxy >= 100 cm": lambda objs: lxy(objs["genAs"]) >= 100,
    },
    "electrons": {
        "pT > 10 GeV": lambda objs: objs["electrons"].pt > 10,
        "|eta| < 1.479": lambda objs: abs(objs["electrons"].eta) < 1.479,
        "1.479 < |eta| < 2.4": lambda objs: ((abs(objs["electrons"].eta) > 1.479)&
                                             (abs(objs["electrons"].eta) < 2.4)),
        "|eta| < 2.4": lambda objs: abs(objs["electrons"].eta) < 2.4,
        "dR(e, A) < 0.5": lambda objs: dR(objs["electrons"], objs["genAs_toE"]) < 0.5,
        #Loose ID = bit 1
        "looseID": lambda objs: check_bit(objs["electrons"].idResults,1),
        "barrel SigmaIEtaIEtaCut": lambda objs: (objs["electrons"].GsfEleFull5x5SigmaIEtaIEtaCut_0) < .0112,
        "barrel DEtaInSeedCut": lambda objs: (abs(objs["electrons"].GsfEleDEtaInSeedCut_0) < .00377),
        "barrel DPhiInCut": lambda objs: (abs(objs["electrons"].GsfEleDPhiInCut_0) < .0884),
        "barrel InverseCut": lambda objs: (objs["electrons"].GsfEleEInverseMinusPInverseCut_0) < .193,
        "barrel Iso": lambda objs: (objs["electrons"].GsfEleRelPFIsoScaledCut_0) < (.112+.506/(objs["electrons"].pt)),
        "barrel ConversionVeto": lambda objs: (abs(objs["electrons"].GsfEleConversionVetoCut_0) == 1),
        "barrel H/E": lambda objs: (objs["electrons"].GsfEleHadronicOverEMEnergyScaledCut_0) < .05,
        "barrel MissingHits": lambda objs: (abs(objs["electrons"].GsfEleMissingHitsCut_0) < 1),
    },
    "muons": {
        #Loose ID = bit 0
        #See https://gitlab.cern.ch/areinsvo/Firefighter/-/blob/master/ffNtuple/plugins/ffNtupleMuon.cc
        "looseID": lambda objs: check_bit(objs["muons"].selectors,0),
        "pT > 5 GeV": lambda objs: objs["muons"].pt > 5,
        "|eta| < 2.4": lambda objs: abs(objs["muons"].eta) < 2.4,
    },
    "photons":{
        "pT > 20 GeV": lambda objs: objs["photons"].pt > 20,
        "|eta| < 2.5": lambda objs: abs(objs["photons"].scEta) < 2.5, # fixme: do we want eta or scEta (which is only available in v4)?
        #Loose ID = bit 0
        "looseID": lambda objs: check_bit(objs["photons"].idResults,0),
    },
    "dsaMuons": {
        "pT > 10 GeV": lambda objs: objs["dsaMuons"].pt > 10,
        "|eta| < 2.4": lambda objs: abs(objs["dsaMuons"].eta) < 2.4,
        "ifcsczero": lambda objs: ak.where(((objs["dsaMuons"].CSCHits==0) 
                                           & (objs["dsaMuons"].DTHits<=18)), False, True),
        "segOverlap < 0.66": lambda objs: objs["dsaMuons"].segOverlapRatio < 0.66,
        "extrapolatedDr > 0.2": lambda objs: objs["dsaMuons"].extrapolatedDr > 0.2,
        "isSubsetAnyPFMuon False": lambda objs: objs["dsaMuons"].isSubsetAnyPFMuon == 0,
        "normChi2 < 4": lambda objs: objs["dsaMuons"].normChi2 < 4,
        "DT + CSC hits > 12": lambda objs: (objs["dsaMuons"].DTHits
                                            + objs["dsaMuons"].CSCHits) > 12,
        "DT + CSC stations >= 2": lambda objs: (objs["dsaMuons"].DTStations
                                                + objs["dsaMuons"].CSCStations) >= 2,
        "ptErrorOverPT < 1": lambda objs: objs["dsaMuons"].ptErrorOverPt <1.0,
    }
}

evt_cut_defs = {
    # This following will be True for every event. There's probably a more intuitive way to do this
    "Keep all evts": lambda objs: ak.num(objs["pvs"]) >= 0,
    "PV filter": lambda objs: ak.num(objs["pvs"]) >= 1,
    "Cosmic veto": lambda objs: objs["cosmicveto"].result,
    ">=2 LJs": lambda objs: ak.num(objs["ljs"]) >= 2,
    ">=2 matched As": lambda objs: ak.num(derived_objs["genAs_matched_lj"](objs, 0.2)) >= 2,
    # 4mu: leading two LJs are both mu-type
    "4mu": lambda objs: ak.count_nonzero(objs["ljs"][:, :2].muon_n >= 2, axis=-1) == 2,
    # 2mu2e: leading two LJs contain exactly 1 mu-type and exactly 1 egm-type
    "2mu2e": lambda objs: ((ak.count_nonzero(objs["ljs"][:, :2].muon_n >= 2, axis=-1) == 1)
                           & (ak.count_nonzero(objs["ljs"][:, :2].muon_n == 0, axis=-1) == 1)),
}
