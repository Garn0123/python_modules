#Author: John Bickel
#Create Date: 09 June 2022
#Purpose: Using the RDKIT_SMILES header, find all unique molecules based on smiles string
# from an input multimol2. 
#Last Edit: 09 June 2022

import argparse, os, sys

#Extracts molecules paired with their smiles strings - also adds a unique index to each
#molecule
def extract_molecule_list_pair_smiles(filename,verbosity):
    molecule_list=[]
    with open(filename,"r") as f:
        found_header=False
        header_start=""
        smiles_string=""
        append_val=0
        line_list=[]
        for line in f:
            #automatically takes the first ### line and will identify all new molecules
            #by that line
            if "###" in line and found_header==False:
                found_header=True
                header_start=line.split(":")[0].strip()+":"

            if header_start in line and found_header == True:
                append_val+=1 #shoddy way of tracking if we're at the second molecule
                if append_val>=2: #if we are, start appending things
                    molecule_list.append([smiles_string,line_list,append_val-2])
                line_list=[]
                line_list.append(line)
            else:
                line_list.append(line)
            if "RD_SMILES:" in line:
                smiles_string=line.split(":")[1].strip()

        molecule_list.append([smiles_string,line_list,append_val])

        if (verbosity):
            print(f"Number of molecules extracted: {len(molecule_list)}.")
    return(molecule_list)

def unique_by_smiles(molecule_list,verbosity):

    #get a unique list of smiles to double check numbers
    uni_smi_list=[x[0] for x in molecule_list]
    unique_numbering=list(set(uni_smi_list))

    list_of_seen=[]
    uni_mols=[]
    dupe_mols=[]
    pos=0
    #checks all molecule SMILES against the unique - makes secondary
    #list showing if they've been seen.
    for entry in molecule_list:
        if entry[0] not in list_of_seen:
            list_of_seen.append(entry[0])
            uni_mols.append(entry[1])
        else:
            dupe_mols.append(entry[1])


    if (verbosity):
        print(f"Number of unique SMILES strings: {len(unique_numbering)}")
        print(f"Number of unique molecules extracted: {len(uni_mols)}")
    return(uni_mols,dupe_mols)

def write_out_mols(molecule_list,duplicate_list,outfile,outfile_duplicates,verbosity):
    with open(outfile,"w") as of:
        for entry in molecule_list:
            for line in entry:
                of.write(line)
    with open(outfile_duplicates,"w") as of:
        for entry in duplicate_list:
            for line in entry:
                of.write(line)
    if(verbosity):
        print(f"Wrote {len(molecule_list)} unique molecules to {outfile}.")
        print(f"Wrote {len(duplicate_list)} duplicate molecules to {outfile_duplicates}.")
        print(f"This means {len(molecule_list) + len(duplicate_list)} molecules were written to file.")
if __name__ == "__main__":
    #setup all the parameters
    program_desc="This script is used to plot a set of footprints with varying types out of output."
    parser=argparse.ArgumentParser(program_desc)
    parser.add_argument("-fi","--input_file",
                        help="Input MOL2 with RDKIT_Smiles in the header.",
                        type=str)
    parser.add_argument("-fo", "--fileout", help="Specifies the name of the unique output file. Defaults to 'unique_out.mol2'",
                         default="unique_out.mol2", type=str)
    parser.add_argument("-v","--verbosity",
                        help="Turns verbose output on.",
                        action="store_true")
    parser.add_argument("-fd","--file_dupes",
                        help="Specifies the name of the duplicate output file. Defaults to 'dupe_mols.mol2'",
                        default="dupe_mols.mol2",type=str)

    args = parser.parse_args()
    if (args.verbosity):
        print("Verbose output turned on. Printing to console.")
    #run the functions
    mol_list=extract_molecule_list_pair_smiles(args.input_file,args.verbosity)
    mol_list,dupe_list=unique_by_smiles(mol_list,args.verbosity)
    write_out_mols(mol_list,dupe_list,args.fileout,args.file_dupes,args.verbosity)