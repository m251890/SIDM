"""Define all commonly used objects"""

from analysis.tools.utilities import dR


obj_defs = {
    # 3=pfmuon, 8=dsamuon, 2=pfelectron, 4=pfphoton
    "mu_ljs": lambda objs: objs["ljs"][(objs["ljs"].muon_n >= 2)],
    "egm_ljs": lambda objs: objs["ljs"][(objs["ljs"].muon_n == 0)],
    "matched_genAs": lambda objs, r: objs["genAs"][dR(objs["genAs"].p4, objs["ljs"].p4) < r]
}
