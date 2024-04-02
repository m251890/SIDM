"""Module to define the Selection and JaggedSelection classes"""

# columnar analysis
from coffea.analysis_tools import PackedSelection
# local
from sidm.definitions.cuts import evt_cut_defs, obj_cut_defs


class Selection:
    """Class to represent the collection of cuts that define a Selection

    A selection consists of event-level cuts which reject whole events.
    Cuts are stored as a PackedSelection.

    All available cuts are defined in sidm.definitions.cuts. The specific cuts that define each
    selection are accepted by Selection() as lists of strings.
    """

    def __init__(self, cuts):
        self.evt_cuts = cuts # list of names of cuts to be applied
        self.all_evt_cuts = PackedSelection() # will be filled later when cuts are evaluated

    def apply_evt_cuts(self, objs, verbose=False):
        """Evaluate all event cuts and apply results to object collections"""

        # evaluate all selected cuts
        for cut in self.evt_cuts:
            if verbose:
                print("Applying cut: ", cut)
            self.all_evt_cuts.add(cut, evt_cut_defs[cut](objs))

        # apply event cuts to object collections
        sel_objs = {}
        for name, obj in objs.items():
            sel_objs[name] = obj[self.all_evt_cuts.all(*self.evt_cuts)]
        return sel_objs


class JaggedSelection:
    """Class to represent the collection of cuts that define a JaggedSelection

    A JaggedSelection consists of object-level cuts (for example, electron or lepton-jet-level cuts).
    Object-level cuts slim object collections and are stored as a dictionary of masks.

    All available cuts are defined in sidm.definitions.cuts. The specific cuts that define each
    selection are accepted by JaggedSelection() as lists of strings.
    """

    def __init__(self, cuts):
        self.obj_cuts = cuts # dict of names of cuts to be applied
        self.evaluated_obj_cuts = {}

    def evaluate_obj_cuts(self, objs, verbose=False):
        """Evaluate all relevant object-level cuts that have not already been evaluated"""
        for obj, cuts in self.obj_cuts.items():
            if obj not in objs:
                print(f"Warning: {obj} not found in sample. "
                      f"The following cuts will not be applied: {cuts}")
                continue
            if obj not in self.evaluated_obj_cuts:
                self.evaluated_obj_cuts[obj] = {}
            for cut in cuts:
                if cut not in self.evaluated_obj_cuts[obj]:
                    if verbose:
                        print("Evaluating ", obj," ",cut)
                    try:
                        self.evaluated_obj_cuts[obj][cut] = obj_cut_defs[obj][cut](objs)
                    except:
                        print(f"Warning: Unable to apply {cut} for {obj}. Skipping.")

    def make_obj_masks(self, channel_cut_list, verbose=False):
        """Create one mask per object, using the subset of cuts specified in channel_cut_list"""
        obj_masks = {}
        for obj, cuts in channel_cut_list.items():
            if obj not in self.evaluated_obj_cuts:
                print(f"Warning: {obj} not found in sample. "
                      f"The following cuts will not be applied: {cuts}")
                continue
            for cut in cuts:
                if cut not in self.evaluated_obj_cuts[obj]:
                    print("Uh oh, haven't evaluated this cut yet! Make sure it was included in the list of cuts you used to initialize this JaggedSelection.  ", obj, ": ",cut)
                else:
                    if verbose:
                        print("Adding the following cut on ",obj, "to the mask: ",cut)
                    if obj not in obj_masks:
                        obj_masks[obj] = self.evaluated_obj_cuts[obj][cut]
                    else:
                        obj_masks[obj] = obj_masks[obj] & self.evaluated_obj_cuts[obj][cut]
        return obj_masks

    def apply_obj_masks(self, objs, obj_masks, verbose=False):
        """Filter object collections based on object masks """
        sel_objs = {}
        for name, obj in objs.items():
            # filter objects if mask exists, return collection unfiltered if mask does not exist
            sel_objs[name] = obj[obj_masks[name]] if name in obj_masks else obj

            if verbose:
                if name in obj_masks:
                    print("Applying mask to collection: ", name)
                else:
                    print("No mask available for collection; returning unfiltered: ",name)
        for collection_to_cut in obj_masks:
            if collection_to_cut not in objs.keys():
                print("WARNING! Trying to apply a cut to ",collection_to_cut," but that's not a valid object")
        return sel_objs
    
    def make_and_apply_obj_masks(self, objs, channel_cut_list, verbose=False):
        return self.apply_obj_masks(objs, self.make_obj_masks(channel_cut_list, verbose), verbose)
