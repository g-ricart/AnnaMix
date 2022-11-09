"""@package AnnaMix
A candidate mixing package for invariant mass spectrum using Ostap.

@author Gabriel RICART gabriel.ricart@etu.sorbonne-universite.fr
@date 2020-07-08
@version 1.0
"""

import ROOT
import warnings
import os.path
import sys
import array
import numbers
from collections import OrderedDict

# Always display warnings.
warnings.simplefilter("always")

class _AnnaOrderedDict(OrderedDict):
    """An OrderedDict with a prepend method.

    Methods
    -------
    prepend(key, value)
        Append given key with given value to the left of the dictionnary.
    """

    def prepend(self, key, value, dict_setitem=dict.__setitem__):
        """Append given key with given value to the left of the dictionnary.

        If given key is already a key of the dictionnary, it is moved to the
        left, and it's value is not updated.

        Parameters
        ----------
        key
            New key to add to the left of the dictionnary.
        value
            Value associated with the new key.
        """

        root = self._OrderedDict__root
        first = root[1]

        if key in self:
            link = self._OrderedDict__map[key]
            link_prev, link_next, _ = link
            link_prev[1] = link_next
            link_next[0] = link_prev
            link[0] = root
            link[1] = first
            root[1] = first[0] = link
        else:
            root[1] = first[0] = self._OrderedDict__map[key] = [root, first, key]
            dict_setitem(self, key, value)

class AnnaMixEvent:
    """Class to perform candidate mixing on invariant mass, pT and rapidity.

    Mixing is done by ordering given TTree by runNumber and eventNumber, mixing
    daughter particles from a dynamic pool of events of size 'train_length'.

    Typical use:
        >>> import ROOT
        >>> import AnnaMix
        >>> raw_tree = ROOT.TTree(...) # TTree containing data to be mixed.
        >>> mix = AnnaMix.AnnaMixEvent(train_length = 50,
                                       tree         = raw_tree,
                                       outfile_path = "output.root")
        >>> mix.addMixCombination("J_psi_1S", ["muplus", "muminus"])
        >>> mix.runMixing()

    The example above will create a "output.root" file containing a TNtuple
    with the following branches :
    J_psi_1S_M, J_psi_1S_PT, J_psi_1S_Y, muplus_M, muplus_PT, muplus_Y,
    muminus_M, muminus_PY, muminus_Y.

    Here, the given TTree must contain the following branches :
    runNumber, eventNumber, muplus_PX, muplus_PY, muplus_PZ, muplus_PZ,
    muplus_PE, muplus_M, muplus_PT, muplus_Y, muminus_PX, muminus_PY,
    muminus_PZ, muminus_PZ, muminus_PE, muminus_M, muminus_PT, muminus_Y.

    A progress indicator can be displayed:
    >>> mix.runMixing(progress = True)
    Processing entry 4583/78596 (6%)

    A verbose mode is available:
    >>> mix.runMixing(verbose = True)
    This mode displays every particles being currently mixed with their index
    in the given TTree, as well as their run and event number.

    Methods
    -------
    addMixCombination(mixed_cdt_name, stems)
        Add a combination of branches to be mixed.
    runMixing(progress=False, verbose=False)
        Method used to run mixing.
    """

    def __init__(self, train_length, tree, outfile_path):
        """Construct all the necessary attributes for the object 'AnnaMixEvent'.

        Parameters
        ----------
        train_length : int
            Number of events used for mixing.
        tree : instance of ROOT.TTree
            TTree to be mixed. Must have 'runNumber' and 'eventNumber'
            as branches.
        outfile_path : str
            Path for the .root file to be saved to.
        """

        self._train_length   = train_length
        self._tree           = tree
        self._tree_index     = self._buildTreeIndex()
        self._outfile_path   = outfile_path
        self._train          = _AnnaOrderedDict() # {runNumber_eventNumber:[indexes]}
        self._branchRegistry = [] # List[str], branches of final tuple.
        self._stems          = [] # List[str], daughter particules names.

    def _buildTreeIndex(self):
        """Build ordered index of self._tree.

        Returns
        -------
        ROOT.TTreeIndex
            Index of self._tree ordered in 'runNumber' and 'eventNumber'.

        Warns
        -----
        UserWarning
            If given TTree doesn't have 'runNumber' and 'eventNumber'
            as branches.
        """

        # Check if runNumber and eventNumber are branches of the given tree.
        for branch in ["runNumber", "eventNumber"]:
            if not branch in self._tree.GetListOfBranches():
                warnings.warn("'{}' is not a branch of the given tree, "
                              "expect problems!".format(branch))

        self._tree.BuildIndex("runNumber", "eventNumber")
        return self._tree.GetTreeIndex()

    def _buildVarlist(self):
        """Build varlist for final tuple.

        Returns
        -------
        str
            The varlist.
        """

        varlist = ""
        for var in self._branchRegistry[:-1]:
            varlist += var
            varlist += ":"
        varlist += self._branchRegistry[-1]
        return varlist

    def _mixAndFill(self, tuple, entry_index):
        """Mix entry corresponding to the given index with the train.

        Paramaters
        ----------
        tuple : ROOT.TNtuple
            TNtuple to fill with mixed candidates and daughters.
        entry_index : long
            Index of the entry to be mixed with the train.
        """

        if self._verbose is True:
            print self._train
            print ("----------------------------------------------------")

        # List of all the indexes in the train.
        mix_cand = []
        for key in self._train.keys():
            mix_cand.extend(self._train[key])

        # print len(self._train.keys())

        # print len(mix_cand)

        # Stems taken from train.
        train_stems = self._stems[1:]

        final_vector = ROOT.TLorentzVector()
        # Need to have a list to keep the info up to saving to the tree.
        mix_vectors = []
        for j in range(len(train_stems)):
            mix_vectors.append(ROOT.TLorentzVector())

        for n in range(len(mix_cand)):
            # Try to get last stem.
            try:
                mix_cand[n+len(train_stems)-1]
            except IndexError:
                continue

            # Get entry associated with given index.
            self._tree.GetEntry(entry_index)

            # Filler used to copy M, PT and Y of daughter particles.
            daughter_filler = [getattr(self._tree,
                                       self._stems[0] + "_M")]
            daughter_filler.append(getattr(self._tree,
                                           self._stems[0] + "_PT"))
            daughter_filler.append(getattr(self._tree,
                                           self._stems[0] + "_Y"))

            # Set current candidate 4-vector.
            final_vector.SetPxPyPzE(getattr(self._tree,
                                            self._stems[0] + "_PX"),
                                    getattr(self._tree,
                                            self._stems[0] + "_PY"),
                                    getattr(self._tree,
                                            self._stems[0] + "_PZ"),
                                    getattr(self._tree,
                                            self._stems[0] + "_PE"))

            if self._verbose is True:
                print self._stems[0], "from", entry_index, getattr(self._tree, "runNumber"), getattr(self._tree, "eventNumber")

            for k, stem in enumerate(train_stems):

                # Get index of train event.
                train_index = mix_cand[n+k]
                self._tree.GetEntry(train_index)

                # Set train candidate 4-vector.
                mix_vectors[k].SetPxPyPzE(getattr(self._tree, stem + "_PX"),
                                          getattr(self._tree, stem + "_PY"),
                                          getattr(self._tree, stem + "_PZ"),
                                          getattr(self._tree, stem + "_PE"))

                # Add candidate to mix 4-vector to final candidate 4-vector.
                final_vector += mix_vectors[k]

                daughter_filler.append(getattr(self._tree, stem + "_M" ))
                daughter_filler.append(getattr(self._tree, stem + "_PT"))
                daughter_filler.append(getattr(self._tree, stem + "_Y" ))

                if self._verbose is True:
                    print stem, "from", train_index, getattr(self._tree, "runNumber"), getattr(self._tree, "eventNumber")

            # Fill tuple:
            # Mixed candidate.
            filler = [final_vector.M(),
                      final_vector.Pt(),
                      final_vector.Rapidity()]
            # Daughter particules.
            filler.extend(daughter_filler)
            # Weight for first stem.
            w = len(mix_cand) - len(train_stems) + 1
            if self._verbose is True:
                print "Weight :", w
                print("----------------------------------------------------")
            filler.append(w)
            tuple.Fill(array.array('f', filler))


    def addMixCombination(self, mixed_cdt_name, stems):
        """Add a combination of branches to be mixed.

        A stem is the name of one of the particles used for mixing. If one want
        to mix a dimuon pair to form J/psi, it is likely that the given TTree
        contains the branches muplus_M, muplus_PZ, ... and muminus_M,
        muminus_PZ, ...
        In this case, the given stems must be ["muplus", "muminus"].

        Branches of the output TNtuple will be named accordingly to the given
        mixed candidate name and stems: <mixed_cdt_name>_<variable>
        and <stem>_<variable>, where <variable> can be Y, PT, M.

        Parameters
        ----------
        mixed_cdt_name : str
            Name of mixed candidate.
        stems : List[str]
            List of stems.

        Warns
        -----
        UserWarning
            If one of the necessary branches for mixing is missing.
        """

        # Check if every element of stems can be used to access PE, PX, PY
        # and PZ branches in the tree.
        branches = [stem + var for stem in stems for var in ["_PE", "_PX",
                                                             "_PY", "_PZ"]]
        for branch in branches :
            if branch not in self._tree.GetListOfBranches():
                warnings.warn("'{}' is not a branch of the given tree, "
                              "expect problems!".format(branch))

        self._stems = stems

        # Fill branch registry.
        vars              = ["_M", "_PT", "_Y"]
        mixed_branches    = [mixed_cdt_name + var for var in vars]
        daughter_branches = [stem + var for stem in stems for var in vars]
        self._branchRegistry.extend(mixed_branches)
        self._branchRegistry.extend(daughter_branches)
        self._branchRegistry.append("w_{}".format(stems[0]))

        daughters_string = stems[0]
        for stem in stems[1:]:
            daughters_string += ", {}".format(stem)

        print("Will mix {} to form {}.".format(daughters_string,
                                               mixed_cdt_name))

    def runMixing(self, progress=False, verbose=False):
        """Method used to run mixing.

        Parameters
        ----------
        progress : bool
            Display a progress indicator.
            Minimal performance impact.
            Default is False.
        verbose : bool
            Display runNumber, eventNumber and index in the given TTree of
            every particles being mixed, as well as the weight value.
            HUGE perfomance impact.
            Default is False.
        """

        self._progress = progress
        self._verbose  = verbose

        # Build varlist for saved tuple.
        varlist = self._buildVarlist()

        # Test
        # print self._branchRegistry
        # print varlist

        # Open file to save tuple.
        with ROOT.TFile(self._outfile_path, "RECREATE") as outfile:
            # Create tuple to be saved.
            new_tuple = ROOT.TNtuple("mt", "mixing tuple", varlist)

            # Initialize counters.
            wagon_counter = 0 # Number of different events processed.

            # Get size of index.
            index_size = self._tree_index.GetN()

            # Wagons.
            evt_wagon  = [] # Wagon of entries to be mixed with train.
            strd_wagon = [] # Stored wagon to be used later.

            # Add first ordered entry to train.
            self._tree.GetEntry(self._tree_index.GetIndex()[0])
            evt_wagon.append((self._tree_index.GetIndex()[0],
                              getattr(self._tree, "runNumber"),
                              getattr(self._tree, "eventNumber")))

            # Print progress if asked.
            if self._progress is True:
                sys.stdout.write("Processing entry 1/{} ({}%)".format(index_size, 100/index_size))

            # Run over ordered entries except the first one.
            for i in range(index_size)[1:]:
                # Print progress if asked.
                if self._progress is True:
                    sys.stdout.write("\rProcessing entry {}/{} ({}%)".format(i+1, index_size, (i+1)*100/index_size))

                cur_tree_index = self._tree_index.GetIndex()[i]
                self._tree.GetEntry(cur_tree_index)

                # Get current entry runNumber and eventNumber.
                cur_entry = (cur_tree_index,
                             getattr(self._tree, "runNumber"),
                             getattr(self._tree, "eventNumber"))

                # Get last entry stored in train.
                prev_entry = evt_wagon[-1]

                # Fill wagon until a new event is met.
                if cur_entry[1:] == prev_entry[1:]:
                    evt_wagon.append(cur_entry)
                    continue

                # When a new event is met, if train is empty.
                elif len(self._train.keys()) < 1:
                    # Counter for number of reversed event.
                    reverse_counter = 1
                    # Fill train with last events from the given tree.
                    # Get last ordered entry of the tree.
                    cur_rev_index = self._tree_index.GetIndex()[index_size-1]
                    self._tree.GetEntry(cur_rev_index)
                    cur_nrun   = getattr(self._tree, "runNumber")
                    cur_nevent = getattr(self._tree, "eventNumber")
                    cur_train_key = "{}_{}".format(cur_nrun, cur_nevent)

                    # Create dictionnary key and add current index.
                    self._train[cur_train_key] = [cur_rev_index]

                    for j in reversed(range(index_size - 1)):
                        # Get index, runNumber and eventNumber.
                        cur_rev_index  = self._tree_index.GetIndex()[j]
                        self._tree.GetEntry(cur_rev_index)
                        cur_nrun   = getattr(self._tree, "runNumber")
                        cur_nevent = getattr(self._tree, "eventNumber")
                        cur_train_key = "{}_{}".format(cur_nrun, cur_nevent)

                        if cur_train_key == self._train.keys()[-1]:
                            # Add current index to the dictionnary.
                            self._train[cur_train_key].append(cur_rev_index)

                            continue

                        elif reverse_counter < self._train_length:
                            # If train is not complete, start another wagon.
                            self._train[cur_train_key] = [cur_rev_index]
                            reverse_counter += 1

                            continue

                        else:
                            # When train is full, mix event wagon with it.
                            for wagon_entry in evt_wagon:
                                self._mixAndFill(new_tuple, wagon_entry[0])

                            # Store and replace event wagon content with
                            # current entry.
                            strd_wagon = evt_wagon
                            evt_wagon  = [cur_entry]

                            wagon_counter += 1

                            break

                    continue

                # When a new event is met, but train is neither empty nor full.
                else : # elif wagon_counter < self._train_length:
                    # Update train:
                    # Dump last wagon of train.
                    del self._train[self._train.keys()[-1]]
                    # Create new key.
                    new_key = "{}_{}".format(strd_wagon[0][1], strd_wagon[0][2])
                    # Add stored wagon to the train.
                    self._train.prepend(new_key, [a[0] for a in strd_wagon])

                    # Mix event wagon with train.
                    for wagon_entry in evt_wagon:
                        self._mixAndFill(new_tuple, wagon_entry[0])

                    # Store and replace event wagon content with current entry.
                    strd_wagon = evt_wagon
                    evt_wagon  = [cur_entry]

                    wagon_counter += 1

                    continue

            # At the end of the tree, mix last event.

            # Update train:
            # Dump last wagon of train.
            del self._train[self._train.keys()[-1]]
            # Create new key.
            new_key = "{}_{}".format(strd_wagon[0][1], strd_wagon[0][2])
            # Add stored wagon to the train.
            self._train.prepend(new_key, [a[0] for a in strd_wagon])

            # Mix event wagon with train.
            for wagon_entry in evt_wagon:
                self._mixAndFill(new_tuple, wagon_entry[0])

            wagon_counter += 1

            try:
                assert new_tuple.GetEntries() > 0
            except AssertionError:
                warnings.warn("Saved tuple is empty !")

            outfile[new_tuple.GetName()] = new_tuple

            if self._progress is True:
                sys.stdout.write("\n")

            print("Mixing done on {} events.".format(wagon_counter))
