import csv
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

# List of common inorganic ions and counterions to remove
INORGANIC_IONS = {"[Na+]", "[Cl-]", "[K+]", "[Ca+2]", "[Mg+2]", "[Fe+2]", "[Fe+3]", "[SO4-2]", "[NO3-]", "[PO4-3]", "[H+]", "[Ga+3]"}

# List of metal atomic numbers to remove (common metals in salts)
METAL_ATOMIC_NUMS = {3, 11, 12, 13, 19, 20, 26, 30, 31}  # Li, Na, Mg, Al, K, Ca, Fe, Zn, Ga

def is_organic(mol):
    """Check if a molecule is organic (contains Carbon and is not an inorganic ion)."""
    if not mol:
        return False
    has_carbon = any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    no_metals = all(atom.GetAtomicNum() not in METAL_ATOMIC_NUMS for atom in mol.GetAtoms())
    return has_carbon and no_metals

def desalt_smiles(smiles: str) -> str:
    """Remove salts, counterions, and inorganic fragments from a given SMILES string."""
    if not smiles:
        return "Invalid_SMILES"
    
    fragments = smiles.split('.')
    mols = [Chem.MolFromSmiles(frag) for frag in fragments if Chem.MolFromSmiles(frag)]
    
    # Remove known inorganic ions explicitly
    mols = [mol for mol in mols if mol and Chem.MolToSmiles(mol) not in INORGANIC_IONS]
    
    # Keep only organic fragments
    organic_mols = [mol for mol in mols if is_organic(mol)]
    
    if not organic_mols:
        return "Invalid_SMILES"  # No valid organic fragments found
    
    # Select the largest organic molecule based on heavy atom count
    largest_mol = max(organic_mols, key=lambda m: rdMolDescriptors.CalcNumHeavyAtoms(m))
    return Chem.MolToSmiles(largest_mol)

# File paths
input_file = "input_smiles.csv"
desalted_file = "desalted_smiles.csv"
not_desalted_file = "not_desalted_smiles.csv"
invalid_file = "invalid_smiles.csv"

# Counters
total_entries = 0
perfectly_desalted = 0
not_desalted = 0
invalid_entries = 0

desalted_data = []
not_desalted_data = []
invalid_data = []

# Read and process SMILES
with open(input_file, "r") as infile:
    reader = csv.reader(infile)
    header = next(reader, None)  # Skip header if present
    
    for row_num, row in enumerate(reader, start=2):  # Start from row 2 (after header)
        if not row or len(row) < 1 or not row[0].strip():  # Skip empty or invalid rows
            continue
        
        total_entries += 1
        original_smiles = row[0].strip()
        cleaned_smiles = desalt_smiles(original_smiles)
        
        if cleaned_smiles == "Invalid_SMILES":
            invalid_entries += 1
            invalid_data.append([row_num, original_smiles])
        elif cleaned_smiles == original_smiles:
            not_desalted += 1
            not_desalted_data.append([row_num, original_smiles])
        else:
            perfectly_desalted += 1
            desalted_data.append([row_num, original_smiles, cleaned_smiles])

# Calculate percentage of valid inputs
valid_percentage = ((perfectly_desalted + not_desalted) / total_entries) * 100 if total_entries > 0 else 0

# Write sorted outputs
def write_csv(filename, data, headers):
    with open(filename, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(headers)  # Header
        writer.writerows(data)

# Save outputs
write_csv(desalted_file, desalted_data, ["Row_Number", "Original_SMILES", "Cleaned_SMILES"])
write_csv(not_desalted_file, not_desalted_data, ["Row_Number", "Original_SMILES"])
write_csv(invalid_file, invalid_data, ["Row_Number", "Original_SMILES"])

# Print analysis results
print(f"Total SMILES processed: {total_entries}")
print(f"Perfectly desalted molecules: {perfectly_desalted}")
print(f"Not desalted (unchanged molecules): {not_desalted}")
print(f"Invalid entries: {invalid_entries}")
print(f"✅ Percentage of valid inputs: {valid_percentage:.2f}%")
print("Sorted results saved in:")
print(f"- {desalted_file} (desalted molecules with row numbers)")
print(f"- {not_desalted_file} (not desalted molecules with row numbers)")
print(f"- {invalid_file} (invalid molecules with row numbers)")

