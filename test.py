import math
import sys
import numpy as np

#DICTIONNARIES FOR HBONDS CALCULATIONS 

side_comb_don = { #combinaison for donnors
                  "ARG": ["NE","NE1","NE2"], 
                  "ASN": ["ND"],
                  "CYS": ["SG"], 
                  "GLN": ["NE"],
                  "HIS": ["ND", "NE"], 
                  "LYS": ["NZ"], 
                  "SER": ["OG"], 
                  "THR": ["OG"], 
                  "TRP": ["NE"], 
                  "TYR": ["OE"], 
}

side_comb_accept = { #combinaison for acceptors
                     "ASN": ["OD"], 
                     "ASP": ["OD1","OD2"], 
                     "CYS": ["SG"], 
                     "GLN": ["OE"], 
                     "GLU": ["OE1", "OE2"], 
                     "HIS": ["ND","NE"], 
                     "MET": ["SD"], 
                     "SER": ["OG"], 
                     "THR": ["OG"], 
                     "TYR": ["OE"]
}

side_accept = { #previous carbon before the donnor
                "CD": ["NE","NE2", "OE1","OE2"],
                "CZ": ["NH1","NH2","OH"], 
                "CG": ["ND1","ND2","OD1","OD2"], 
                "CB": ["OG","OG1"]
}

side_don = { #Hydrogens associates with the donnors
             "NE": ["HE"],
             "NH1": ["HH11","HH12"],
             "NH2": ["HH21","HH22"],
             "ND2": ["HD21","HD22"],
             "NE2": ["HE21","HE22"],
             "NZ": ["HZ2"],
             "OG": ["HG"],
             "NE1": ["HE1"],
             "OH": ["HH"],
             "OG1": ["HG1"]
}

#Parsing the pdb file    

def parse_pdb(filename):
    """Parse the PDB file

    Parameters
    ----------
    filename: str
        PDB filename.

    Returns
    -------
    atom_info: dict
        keys: tuple with atom_name, res_name, res_seqnum, chain_id
        values: x, y, z coordinates
    """
    atom_info = {}

    with open(filename, "r") as pdb_input:
        for line in pdb_input:
            if line.startswith("ATOM"):
                # Extract pattern from PDB
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                res_seqnum = int(line[22:26].strip())
                chain_id = line[21:22].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                atom_info[(atom_name, res_name, res_seqnum, chain_id)] = np.array([x, y, z])

    return atom_info

#Calculating the distance

def calc_distance(atom1, atom2):
    """
    Calculate the Euclidean distance between two 3D points.
    
    Args:
    atom1 (list): List containing the x, y, and z coordinates of the first atom.
    atom2 (list): List containing the x, y, and z coordinates of the second atom.
    
    Returns:
    float: The Euclidean distance between the two atoms.
    """
    distance = np.sqrt((atom1[0] - atom2[0]) ** 2 +
                       (atom1[1] - atom2[1]) ** 2 +
                       (atom1[2] - atom2[2]) ** 2)
    
    return distance

#Matrix for the distance

def distance_matrix(atom_info):
    """
    Create a distance matrix for a set of atoms.
    
    Args:
    atom_info (dict): A dictionary where keys are atom names and values are lists containing
                      the x, y, and z coordinates of the atoms.
    
    Returns:
    dict: A nested dictionary representing the distance matrix between atoms.
    """
    dist_matrix = {}

    for atom1, coord1 in atom_info.items():
        dist_matrix[atom1] = {}
        for atom2, coord2 in atom_info.items():
            dist_matrix[atom1][atom2] = calc_distance(coord1, coord2)
    
    return dist_matrix

#Angle calculation

def calc_angle(atom1, atom2, atom3):

    vector1 = np.array([atom1[0]-atom2[0],atom1[1]-atom2[1],atom1[2]-atom2[2]])
    vector2 = np.array([atom3[0]-atom2[0],atom3[1]-atom2[1],atom3[2]-atom2[2]])

    cosine_angle = np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
    angle = np.arccos(cosine_angle)
    degrees = np.degrees(angle)
    return degrees
    
def hydrophobic_interactions(dist_matrix, atom_info):
    interacting_res = []

    atom_sides = ["CB", "CD", "CD1", "CD2","CE","CE1","CE2","CE3", "CZ", "CZ1","CZ2", "CZ3", "CG", "CG1", "CG2", "CH2", "SD"]  # Atoms in the side chains
    apolar_residues = ["ALA", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL", "PRO"]

     # Get the number of atoms

    for atom1, _ in dist_matrix.items():
        for atom2, _ in dist_matrix.items():
            distance = dist_matrix[atom1][atom2]  # Calculate the distance within the loop
        
            if (
                atom1[1] in apolar_residues
                and atom1[0] in atom_sides
                and atom2[1] in apolar_residues
                and atom2[0] in atom_sides
                and atom1[2] != atom2[2]
                and distance < 5.0
            ):  
                interacting_res_key = (atom1[1], atom2[1], atom1[3])
                if interacting_res_key not in interacting_res:
                   interacting_res.append(interacting_res_key)
                   print(
                         f"POSITION: {atom1[2]} RESIDUE: {atom1[1]}, CHAIN: {atom1[3]}, "
                         f"POSITION: {atom2[2]} RESIDUE: {atom2[1]}, CHAIN: {atom2[3]}"
                   )



def disulfite_bridges(dist_matrix, atom_info):
    
    interacting_res = []

    for atom1, _ in dist_matrix.items():   
        for atom2, _ in dist_matrix.items():
            distance = dist_matrix[atom1][atom2]
     
            if (
                atom1[0] == "SG"
                and atom1[1] == "CYS"
                and atom1[2]<= atom2[2]
                and distance <= 2.2
            ):  
                interacting_res_key = (atom1[2], atom2[2], atom1[3])
                if interacting_res_key not in interacting_res:
                   interacting_res.append(interacting_res_key)
                   print(
                         f"POSITION: {atom1[2]} RESIDUE: {atom1[1]}, CHAIN: {atom1[3]}, "
                         f"POSITION: {atom2[2]} RESIDUE: {atom2[1]}, CHAIN: {atom2[3]}, "
                         f"DISTANCE: {dist_matrix[atom1][atom2]}"
                   )
                   

                   
def ionic_interactions(dist_matrix, atom_info):
    cationic_res = ["ARG", "LYS", "HIS"]
    anionic_res = ["ASP", "GLU"]
    ion_atoms = ["ND1", "NZ", "CZ", "OD2", "OE2"]
    interacting_residues = []

    for atom1, _ in dist_matrix.items():   
        for atom2, _ in dist_matrix.items():
            distance = dist_matrix[atom1][atom2]

            if (atom1[1] in cationic_res and atom1[0] in ion_atoms
                and atom2[1] in anionic_res and atom2[0] in ion_atoms) or \
               (atom2[1] in cationic_res and atom2[0] in ion_atoms
                and atom1[1] in anionic_res and atom1[0] in ion_atoms):
                if distance <= 6 and atom2[2]>atom1[2]:
                   interacting_res_key = (atom1[2], atom2[2], atom1[3])
                   if interacting_res_key not in interacting_residues:
                      interacting_residues.append(interacting_res_key)
                      print(
                            f"POSITION: {atom1[2]} RESIDUE: {atom1[1]}, CHAIN: {atom1[3]}, "
                            f"POSITION: {atom2[2]} RESIDUE: {atom2[1]}, CHAIN: {atom2[3]}"
                      )

def hbond_main_main(dist_matrix, atom_info):
    
    for atom1, _ in dist_matrix.items():
        for atom2, _ in dist_matrix.items():
            distance = dist_matrix[atom1][atom2]
            if (
                atom1[0] == "N"
                and atom2[0] in {"O","OXT"}
                and distance < 3.5
                and abs(atom2[2]-atom1[2])>=2
                and atom2[1] != "PRO"
                ):
                 # Define the relevant atoms for clarity
                 donor = atom_info[atom1]
                 acceptor = atom_info[atom2]
                 c_acceptor = atom_info["C", atom2[1], atom2[2], atom2[3]]
                 h_donor = atom_info["H", atom1[1], atom1[2], atom1[3]]

                 # Calculate distances and angles
                 dist_h_acceptor = dist_matrix["H", atom1[1], atom1[2], atom1[3]][atom2[0], atom2[1], atom2[2], atom2[3]]
                 angle_dha = calc_angle(donor, h_donor, acceptor)
                 angle_adc = calc_angle(acceptor, donor, c_acceptor)
                    
                 print(
                      f"POSITION: {atom1[2]} RESIDUE: {atom1[1]}, CHAIN: {atom1[3]}, ATOM : {atom1[0]}, "
                      f"POSITION: {atom2[2]} RESIDUE: {atom2[1]}, CHAIN: {atom2[3]}, ATOM : {atom2[0]},  "
                      f"Dd-a: {distance}, Dh-a {dist_h_acceptor}, "
                      f"A(d-H-N) : {angle_dha}, A(a-O=C) : {angle_adc}"
                      ) 
                    # Further processing or printing can be added here

def hbond_side_side(dist_matrix, atom_info):
    for atom1, _ in dist_matrix.items():
        for atom2, _ in dist_matrix.items():
            distance = dist_matrix[atom1][atom2]

            if (
                atom1[0] in side_comb_don.keys()
                and atom1[1] in side_comb_don.values()
                and atom2[0] in side_comb_accept.keys()
                and atom2[1] in side_comb_accept.values()
            ):
                donor = atom_info[atom1]
                acceptor = atom_info[atom2]

                if donor is not None and acceptor is not None:
                    if (
                        (abs(donor[2] - acceptor[2]) >= 2
                         and (donor[0] == "SG" or acceptor[0] == "SD")
                         and distance <= 4
                         )
                         or distance <= 3.5
                    ):
                        c_accept_name = side_accept.get(acceptor[0], [])
                        c_acceptor = atom_info.get((c_accept_name, acceptor[1], acceptor[2], acceptor[3]))
                        h_name = side_don.get(donor[0], [])
                        h = atom_info.get((h_name, donor[1], donor[2], donor[3]))

                        if h is not None and c_acceptor is not None:
                            dist_don_accept = dist_matrix.get((donor, acceptor), None)
                            dist_h_accept = dist_matrix.get((h, acceptor), None)
                            angle_accept = calc_angle(acceptor, donor, c_acceptor)
                            angle_don = calc_angle(donor, h, acceptor)

                            if dist_don_accept is not None and dist_h_accept is not None:
                                print(
                                    f"POSITION: {donor[2]} RESIDUE: {donor[1]}, CHAIN: {donor[3]}, ATOM : {donor[0]}, "
                                    f"POSITION: {acceptor[2]} RESIDUE: {acceptor[1]}, CHAIN: {acceptor[3]}, ATOM : {acceptor[0]},  "
                                    f"Dd-a: {dist_don_accept}, Dh-a {dist_h_accept}, "
                                    f"A(d-H-N) : {angle_don}, A(a-O=C) : {angle_accept}"
                                )
							 

def hbond_main_side(dist_matrix, atom_info):

    for (atom1_key, atom2_key), distance in dist_matrix.items():
        atom1 = atom_info.get(atom1_key)  # Retrieve atom information
        atom2 = atom_info.get(atom2_key)
        donor = None
        acceptor = None
        h = None
        c_acceptor = None  # Initialize c_acceptor
        
        # Check the donor-acceptor pairs
        if atom1[0] == "N" and (atom2[0] in side_comb_accept.keys() and atom2[1] in side_comb_accept[atom2[0]]):
            donor = atom1
            acceptor = atom2
        elif atom1[0] in ["O", "OXT"] and (atom2[0] in side_comb_donn.keys() and atom2[1] in side_comb_accept[atom2[0]]):
            donor = atom2
            acceptor = atom1
        elif atom2[0] == "N" and (atom1[0] in side_comb_accept.keys() and atom1[1] in side_comb_accept[atom1[0]]):
            donor = atom2
            acceptor = atom1
        elif atom2[0] in ["O", "OXT"] and (atom1[0] in side_comb_donn.keys() and atom1[1] in side_comb_accept[atom1[0]]):
            donor = atom1
            acceptor = atom2
        else:
            continue
        if (
            donor is not None
            and acceptor is not None
            and donor[1] != "PRO"
            and abs(donor[2] - acceptor[2]) != 1
        ):
            if (donor[0] == "SG" or acceptor[0] == "SD") and distance <= 4.0 or distance <= 3.9:
                if acceptor[0] in side_comb_accept.keys():
                    c_accept_name = side_accept.get(acceptor[0], [])
                    c_acceptor = atom_info.get(
                        (c_accept_name, acceptor[1], acceptor[2], acceptor[3])
                    )
                else:
                    c_acceptor = atom_info.get(("C", acceptor[1], acceptor[2], acceptor[3]))
                if donor[0] in side_comb_accept.keys():
                    h_name = side_don.get(donor[0], [])
                    h = atom_info.get(
                        (h_name, donor[1], donor[2], donor[3])
                    )
                else:
                    h = atom_info.get(("H", donor[1], donor[2], donor[3]))
                 
                if h is not None and c_acceptor is not None:
                    dist_don_accept = dist_matrix[donor][acceptor]
                    dist_h_accept = dist_matrix[h][acceptor]
                    angle_accept = calc_angle(acceptor, donor, c_acceptor)
                    angle_don = calc_angle(donor, h, acceptor)
                    
                    print(
                        f"POSITION: {donor[2]} RESIDUE: {donor[1]}, CHAIN: {donor[3]}, ATOM: {donor[0]}, "
                        f"POSITION: {acceptor[2]} RESIDUE: {acceptor[1]}, CHAIN: {acceptor[3]}, ATOM: {acceptor[0]},  "
                        f"Dd-a: {dist_don_accept}, Dh-a: {dist_h_accept}, "
                        f"A(d-H-N): {angle_don}, A(a-O=C): {angle_accept}"
                    )


def main():

    filename = input("Enter the file path = ")
    pdb_name = filename[:-6]
    atom_info = parse_pdb(filename)

    # Calculate the distance matrix
    dist_matrix = distance_matrix(atom_info)
    
    # Analyze various interactions
    
    print("\nHydrophobic Interactions:")
    hydrophobic_interactions(dist_matrix, atom_info)
    
    print("\nDisulfite bridges:")
    disulfite bridges(dist_matrix, atom_info)
    
    print("\nIonic Interactions:")
    ionic_interactions(dist_matrix, atom_info)
    
    print("\nHydrogen Bonds (Main Chain - Main Chain):")
    hbond_main_main(dist_matrix, atom_info)
    

    

if __name__ == "__main__":
    main() 
































    
        
        
 
                           


    


