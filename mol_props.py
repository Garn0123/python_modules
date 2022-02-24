# Import modules for dealing with chemical information
from rdkit import Chem as Chem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem as Chem2
from rdkit.Chem import Descriptors as Desc
from rdkit.Chem import Draw
from rdkit.Chem import rdmolfiles as RDFile

#Stuff for calculation of SASA by RDKit
import pickle
import os
import os.path as op
import math
from collections import defaultdict

#itemgetter for sorting purposes
from operator import itemgetter

#reads in a data file of single values - must be fed a datatype
#and will return correctly types values. 
def read_data_file(filename,datatype):
    line_list=[]
    with open(filename,"r") as f:
        for line in f:
            if datatype=="float":
                line_list.append(float(line.strip()))
            if datatype=="int":
                line_list.append(int(line.strip()))
            if datatype=="str":
                line_list.append(str(line.strip()))
    return(line_list)
    
#writes a data file - essentially just takes a list and writes 
#each entry of the list to an individual line.
def write_data_file(filename,data_list):
    with open(filename,"w") as f:
        for entry in data_list:
            f.write(str(entry)+"\n")

#Counts the number of mols in a given file based on the @<TRIPOS>MOLECULE line
def countMolsInFile(filename):
    mol_counter=0
    with open(filename,"r") as f:
        for line in f:
            if "@<TRIPOS>MOLECULE" in line:
                mol_counter+=1
    return(mol_counter)

#Takes in a list of list of fragments (mols) extracted from a file and
#the index of the line containing the fragname (lnk.0, scf.1, etc)
#and combines them into a sortable list [[frag, index]]
def combine_fragments_by_name(all_fragment_lists,name_line_id):
    combined_fragments_list=[]
    for frag_list in all_fragment_lists:
        for entry in frag_list:
            freq_pos=int(entry[name_line_id].split(".")[name_line_id].strip())
            combined_fragments_list.append([entry,freq_pos])
    return(combined_fragments_list)

#Takes fragments indexed by combine_fragments functions and writes them
#to a set of files with a given name scheme, starting from 0
def write_fragments(fragments_list, num_write, file_prefix):
    fragments_list=sorted(fragments_list, key=itemgetter(1), reverse=False)

    for i in range(0, num_write):
        with open("%s%s.mol2" % (file_prefix, i),"w") as f:
            for line in fragments_list[i][0]:
                f.write(line)


#Takes in a filename and the first entry of a header and returns a list
#of molecules with headers.
def extract_molecule_list(filename,header_start):
    line_list=[]
    molecule_list=[]
    with open(filename,"r") as f:
        append_val=0
        for line in f:
            #TYPE: is the first line in each header
            if header_start in line:
                append_val+=1 #shoddy way of tracking if we're at the second molecule
                if append_val>=2: #if we are, start appending things
                    append_next=True
                    molecule_list.append(line_list)
                    
                line_list=[]
                line_list.append(line)
            else:
                line_list.append(line)
        molecule_list.append(line_list)
    return(molecule_list)

#Takes in a filename and the score as listed in the header and returns
#a list of scores.
def extract_scores_from_file(filename,score_type):
    score_list=[]
    with open(filename,"r") as f:
        for line in f:
            if score_type in line:
                score_list.append(float(line.split(":")[1].strip()))

    return(score_list)

#Takes a list of molecules extracted from a file and returns
#a list of scores pulled from the header.
def extract_scores_from_list(list_of_molecules,score_type):
    score_list=[]
    for molecule in list_of_molecules:
        for line in molecule:
            if score_type in line:
                score_list.append(float(line.split(":")[1].strip()))
    return(score_list)



def separate_top_molecules(list_of_molecules,score_type,num_sep):
    sep_molecules=[]
    mol_scores=[]
    for i,entry in enumerate(list_of_molecules):
        for line in entry:
            if score_type in line:
                mol_scores.append((i,float(line.split(":")[1].strip())))
    mol_scores.sort(key = lambda x: x[1])
    for i in range(0,num_sep):
        sep_molecules.append(list_of_molecules[mol_scores[i][0]])
    return(sep_molecules)




# Takes a list of rdkitmols and returns the scores for graphing.
# Current scores included in this function are:
# Molecular Weight, LogP, #Aromatic Rings, #HDonors,
# #HAcceptors, Synth Accessiblity, 
def calc_rdkit_scores_for_graphing(rdkit_mols):
#########################
#Calculated by QED code
    mw_list=[]
    logp_list=[]
    aro_ring_list=[]
    hbd_list=[]
    hba_list=[]
    rot_bond_list=[]
    pains_alerts_list=[]
    qed_list=[]
#########################
    sasa_list=[]
    spiro_list=[]
    stereo_list=[]

    for mol in rdkit_mols:
        #calculates QED properties and then farms out those numbers
        #to the various lists needed for graphing
        QED_props_list=Chem.QED.properties(mol)
        mw_list.append(QED_props_list[0])
        logp_list.append(QED_props_list[1])
        hba_list.append(QED_props_list[2])
        hbd_list.append(QED_props_list[3])
        rot_bond_list.append(QED_props_list[5])
        aro_ring_list.append(QED_props_list[6])
        pains_alerts_list.append(QED_props_list[7])
        qed_list.append(Chem.QED.qed(mol,QED_props_list))
        sasa_list.append(calculateSAScore(mol))
        spiro_list.append(rdMolDescriptors.CalcNumSpiroAtoms(mol))
        stereo_list.append(rdMolDescriptors.CalcNumAtomStereoCenters(mol))
    return [mw_list, logp_list, hba_list, hbd_list, rot_bond_list, \
            aro_ring_list, pains_alerts_list, qed_list, sasa_list, \
            stereo_list, spiro_list]

#This one is for a list molecules, returns by molecule
def calc_rdkit_scores_list_mol_for_graphing(rdkit_mols):
    stats_by_mol=[]
    for mol in rdkit_mols:
        #calculates QED properties and then appends them as a group so they
        # are all associated with the properties for each molecule
        QED_props_list=Chem.QED.properties(mol)
        mw = QED_props_list[0]
        logp = QED_props_list[1]
        hba = QED_props_list[2]
        hbd = QED_props_list[3]
        rot_bond = QED_props_list[5]
        aro_ring = QED_props_list[6]
        struct_alerts = QED_props_list[7]
        #qed = Chem.QED.qed(mol,QED_props_list)
        qed = Chem.QED.qed(mol)
        sasa = calculateSAScore(mol)
        spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        stereo = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        stats_by_mol.append([mw, logp, hba, hbd, rot_bond, aro_ring,\
                            struct_alerts, qed, sasa, spiro, stereo])
    return stats_by_mol

######### This one is for a single molecule calculation.
def calc_rdkit_scores_single_mol_for_graphing(mol):
    stats_by_mol=[]
    #calculates QED properties and then appends them as a group so they
    # are all associated with the properties for each molecule
    QED_props_list=Chem.QED.properties(mol)
    mw = QED_props_list[0]
    logp = QED_props_list[1]
    hba = QED_props_list[2]
    hbd = QED_props_list[3]
    rot_bond = QED_props_list[5]
    aro_ring = QED_props_list[6]
    struct_alerts = QED_props_list[7]
    psa = QED_props_list[4]
    #qed = Chem.QED.qed(mol,QED_props_list)
    qed = Chem.QED.qed(mol)
    sasa = calculateSAScore(mol)
    spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    stereo = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    return([mw, logp, hba, hbd, rot_bond, aro_ring,\
            struct_alerts, qed, sasa, spiro, stereo, psa])

def calc_spiro_stereo(rdkit_mols):
    spiro_list=[]
    stereo_list=[]
    for mol in rdkit_mols:
        spiro_list.append(rdMolDescriptors.CalcNumSpiroAtoms(mol))
        stereo_list.append(rdMolDescriptors.CalcNumAtomStereoCenters(mol))
    return((spiro_list,stereo_list))


###TANIMOTO COMPARISON OF FRAGMENTS
def pw_tan_return_index(master_mol,compare_list):
    for i in range(0,len(compare_list)):
        if (DataStructs.FingerprintSimilarity(master_mol,compare_list[i]) == 1):
            return(i)


def pw_tan_print_tan(master_mol,compare_list):
    for i in range(0,len(compare_list)):
        if (DataStructs.FingerprintSimilarity(master_mol,compare_list[i]) > 0.6):
            print(DataStructs.FingerprintSimilarity(master_mol,compare_list[i]))
            print(i)


def du_to_h(mol_list):
    for i in range(0,len(mol_list)):
        for x in range(0,len(mol_list[i])):
            if "Du" in mol_list[i][x]:
                line_change=mol_list[i][x][:47] + "H " + mol_list[i][x][49:]
                mol_list[i][x] = line_change
    return(mol_list)

#This function is currently built to graph the parameters extracted 
#by the calc_rdkit_scores functions elsewhere in this library
#SPLITS UP GRAPHING FOR ALL INTEGER VALUES OF SASA
#Saves file in experiment_list[0] directory
def graph_9panel_mol_props_by_sasa(experiment_list,workdir):
    graph_values=["MW","LOGP","HBA","HBD","ROT","ARO","PAINS","QED","SASA","STER","SPIR"]
    graph_titles=["Molecular Weight",\
                "ALogP",\
                "HBond Acceptors",\
                "HBond Donors",\
                "# Rotatable Bonds",\
                "# Aromatic Rings",\
                "# Structural Alerts",\
                "QED",\
                "Novartis Synthetic Accessibility",\
                "Stereo Centers",\
                "Spiro Atoms"]
    for experiment in experiment_list:
        data_list=[]
        #data extraction from previously saved data
        with open(f"{workdir}/{experiment}/{experiment}_props.dat") as f:
            for line in f:
                data_list.append(line.strip().split(","))
        #converts all values to floats
        for i,entry in enumerate(data_list):
            data_list[i] = [float(x) for x in entry]
        data_collected.append(data_list)

    #this function graphs by SASA, which here is 0 to 10
    for i in range(0,10):
        fig,axs = plt.subplots(4,3,figsize=(20,15),dpi=500,sharey=True)
        pos_x=0
        pos_y=0
        graph_val_position=0
        pos = np.arange(len(experiment_list)+1,1,-1)
        for val_ind, value in enumerate(graph_values):
            data_list=[]
            #separates out all molecule data by the SynthA
            for exp_ind, experiment in enumerate(experiment_list):
                #data_list.append([x[val_ind] for x in data_collected[exp_ind]])
                data_list.append([x[val_ind] for x in data_collected[exp_ind] \
                                if x[8] > i and x[8] < i+1])
            for index, entry in enumerate(data_list):
                #Just in the case there are NO MOLECULES MADE for this value, 
                if len(entry) == 0:
                    data_list[index]=[float('nan'), float('nan')]
            if value == "MW":
                plt.setp(axs[pos_x,pos_y], xlim=(0,850))
                axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
            if value == "LOGP":
                plt.setp(axs[pos_x,pos_y], xlim=(-8,10))
                axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
            if value == "HBA" or value == "ROT":
                plt.setp(axs[pos_x,pos_y], xlim=(0,16))
                axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
            if (value == "HBD") or (value == "ARO"):
                plt.setp(axs[pos_x,pos_y], xlim=(0,10))
                axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
            if value == "SASA":
                plt.setp(axs[pos_x,pos_y], xlim=(0,10))
                axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
            if (value == "PAINS") or (value == "STER") or (value =="SPIR"):
                plt.setp(axs[pos_x,pos_y], xlim=(0,5))
                axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
            if value == "QED":
                plt.setp(axs[pos_x,pos_y], xlim=(0,1))
                axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
            axs[pos_x,pos_y].set_title(graph_titles[graph_val_position])
            axs[pos_x,pos_y].set_yticks(pos)
            axs[pos_x,pos_y].set_yticklabels(experiment_list)
            if (pos_y > 0) and (pos_y % 2 == 0):
                pos_x+=1
                pos_y=0
            else:
                pos_y+=1
            graph_val_position+=1
            fig.subplots_adjust(top=0.94,hspace = 0.2)
        st = fig.suptitle(f"{experiment[0]}_{i}_to_{i+1}",fontsize=20)
        st.set_y(1.00)
        plt.savefig(f"{workdir}/{experiment[0]}/zzz.{experiment[0]}_sasa_{i}.png")
        plt.close()

# Function that reads multi-molecule MOL2 files. Adapted from:
# https://chem-workflows.com/articles/2019/07/18/building-a-multi-molecule-mol2-reader-for-rdkit/
# further adapted from function written by Guilherme Duarte, Rizzo Lab
# by John Bickel, Rizzo Lab

#faster than reading in the file and doing it in-place. 
#extract and use this function instead.
def convert_list_to_rdkitmol(molecule_list):
    #  This function extracts all the molecules in the multi-molecule
    #  MOL2 file `file` and returns a list of rdkit.Chem.rdchem.mol 
    #  object.
    #  
    #  Variable         I/O          dtype           default value?
    #  ------------------------------------------------------------
    #  file              I           string                  None
    #  mols              O           list                    N/A
    #  mols[i]           O           rdkit.Chem.rdchem.mol   N/A


    # does some cleaning in case there's any header information - 
    # This will pull everything from @TRIPOS<MOLECULE> to ROOT
    conversion_list=[]

    for molecule in molecule_list:
        temp_mol=[]
        record=False
        #runs through all the lines per molecule
        for line in molecule:
            #checks for first line
            if ("@<TRIPOS>MOLECULE") in line:
                record=True
                temp_mol.append(line)
            #will record if @<TRIPOS>MOLECULE found
            elif record==True:
                temp_mol.append(line)
            #if root is in line, we're done with this molecule for most small molecules
            if ("ROOT") in line:
                #block formatting
                block = ",".join(temp_mol).replace(',','')
                #conversion step
                m=Chem.MolFromMol2Block(block,
                                        sanitize=False,
                                        cleanupSubstructures=False)
                conversion_list.append(m)
    for mol in conversion_list:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    return(conversion_list)

def convert_mol_to_rdkitmol(mol):
    #  This function extracts all the molecules in the multi-molecule
    #  MOL2 file `file` and returns a list of rdkit.Chem.rdchem.mol 
    #  object.
    #  
    #  Variable         I/O          dtype           default value?
    #  ------------------------------------------------------------
    #  file              I           string                  None
    #  mols              O           list                    N/A
    #  mols[i]           O           rdkit.Chem.rdchem.mol   N/A


    # does some cleaning in case there's any header information - 
    # This will pull everything from @TRIPOS<MOLECULE> to ROOT
    conversion_list=[]

    temp_mol=[]
    record=False
    #runs through all the lines per molecule
    for line in mol:
        #checks for first line
        if ("@<TRIPOS>MOLECULE") in line:
            record=True
            temp_mol.append(line)
        #will record if @<TRIPOS>MOLECULE found
        elif record==True:
            temp_mol.append(line)
        #if root is in line, we're done with this molecule for most small molecules
        if ("ROOT") in line:
            #block formatting
            block = ",".join(temp_mol).replace(',','')
            #conversion step
            m=Chem.MolFromMol2Block(block,
                                    sanitize=False,
                                    cleanupSubstructures=False)
            m.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(m,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    
            return(m)



# Function that reads multi-molecule MOL2 files. Adapted from:
# https://chem-workflows.com/articles/2019/07/18/building-a-multi-molecule-mol2-reader-for-rdkit/
# further adapted from function written by Guilherme Duarte, Rizzo Lab
# by John Bickel, Rizzo Lab

#THIS STEP IS SLOWER TO DO IN-PLACE THAN EXTRACTING AND THEN USING
#CONVERT_LIST_TO_RDKITMOL
def mol2_file_to_rdkitmol(filename):
    #  This function extracts all the molecules in the multi-molecule
    #  MOL2 file `file` and returns a list of rdkit.Chem.rdchem.mol 
    #  object.
    #  
    #  Variable         I/O          dtype           default value?
    #  ------------------------------------------------------------
    #  file              I           string                  None
    #  mols              O           list                    N/A
    #  mols[i]           O           rdkit.Chem.rdchem.mol   N/A

    conversion_list=[]
    with open(filename,"r") as f:
        temp_mol=[]
        record=False
        for line in f:
            #checks for first line
            if ("@<TRIPOS>MOLECULE") in line:
                record=True
                temp_mol.append(line)
            #will record if @<TRIPOS>MOLECULE found
            elif record==True:
                temp_mol.append(line)
            #if root is in line, we're done with this molecule for most small molecules
            if ("ROOT") in line:
                #block formatting
                block = ",".join(temp_mol).replace(',','')
                #conversion step
                m=Chem.MolFromMol2Block(block,
                                        sanitize=False,
                                        cleanupSubstructures=False)
                conversion_list.append(m)
    for mol in conversion_list:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    return(conversion_list)






###### FUNCTIONS FOR CALCULATION OF SYNTHETIC ACCESSIBLITY ########
_fscores = None


def readFragmentScores(name='fpscores'):
    import gzip
    global _fscores
    # generate the full path filename:
    if name == "fpscores":
        print(os.getcwd())
        name = op.join(os.getcwd(), name)
    data = pickle.load(gzip.open('%s.pkl.gz' % name))
    outDict = {}
    for i in data:
        for j in range(1, len(i)):
            outDict[i[j]] = float(i[0])
    _fscores = outDict


def numBridgeheadsAndSpiro(mol, ri=None):
    nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    return nBridgehead, nSpiro


def calculateSAScore(m):
    if _fscores is None:
        readFragmentScores()
    m.UpdatePropertyCache(strict=False)
    # fragment score
    fp = rdMolDescriptors.GetMorganFingerprint(m,
                                               2)  # <- 2 is the *radius* of the circular fingerprint
    fps = fp.GetNonzeroElements()
    score1 = 0.
    nf = 0
    for bitId, v in fps.items():
        nf += v
        sfp = bitId
        score1 += _fscores.get(sfp, -4) * v
    score1 /= nf

    # features score
    nAtoms = m.GetNumAtoms()
    nChiralCenters = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
    ri = m.GetRingInfo()
    nBridgeheads, nSpiro = numBridgeheadsAndSpiro(m, ri)
    nMacrocycles = 0
    for x in ri.AtomRings():
        if len(x) > 8:
            nMacrocycles += 1

    sizePenalty = nAtoms**1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1)
    spiroPenalty = math.log10(nSpiro + 1)
    bridgePenalty = math.log10(nBridgeheads + 1)
    macrocyclePenalty = 0.
    # ---------------------------------------
    # This differs from the paper, which defines:
    #  macrocyclePenalty = math.log10(nMacrocycles+1)
    # This form generates better results when 2 or more macrocycles are present
    if nMacrocycles > 0:
        macrocyclePenalty = math.log10(2)

    score2 = 0. - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

    # correction for the fingerprint density
    # not in the original publication, added in version 1.1
    # to make highly symmetrical molecules easier to synthetise
    score3 = 0.
    if nAtoms > len(fps):
        score3 = math.log(float(nAtoms) / len(fps)) * .5

    sascore = score1 + score2 + score3

    # need to transform "raw" value into scale between 1 and 10
    min = -4.0
    max = 2.5
    sascore = 11. - (sascore - min + 1) / (max - min) * 9.
    # smooth the 10-end
    if sascore > 8.:
        sascore = 8. + math.log(sascore + 1. - 9.)
    if sascore > 10.:
        sascore = 10.0
    elif sascore < 1.:
        sascore = 1.0
    return sascore
