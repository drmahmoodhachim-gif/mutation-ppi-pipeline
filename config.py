"""Configuration and tissue-specific data for the mutation pipeline."""

# Gene symbol -> UniProt ID mapping (common cardiac/channel genes)
GENE_UNIPROT = {
    "SCN5A": "Q14524",
    "KCNQ1": "P51787",
    "KCNH2": "Q12809",
    "MYH7": "P12883",
    "TNNT2": "P45379",
    "TNNI3": "P19429",
    "MYBPC3": "Q14896",
    "LMNA": "P02545",
    "RYR2": "Q92736",
    "CACNA1C": "Q13936",
    "FGF13": "Q92913",
    "CALM1": "P0DP23",
    "VCL": "P18206",  # Vinculin
}

# Tissue-specific protein interactors for cardiac myocytes
CARDIAC_MYOCYTE_INTERACTORS = {
    "SCN5A": [
        {"partner": "FGF13", "uniprot": "Q92913", "role": "Late Na+ current regulation, FHF binding"},
        {"partner": "Calmodulin", "uniprot": "P0DP23", "role": "IQ domain binding, inactivation"},
        {"partner": "SCN1B", "uniprot": "Q07699", "role": "β1 subunit, trafficking"},
        {"partner": "SCN2B", "uniprot": "Q9UQ35", "role": "β2 subunit"},
        {"partner": "SCN3B", "uniprot": "Q9NY72", "role": "β3 subunit"},
        {"partner": "SCN4B", "uniprot": "Q9H0B8", "role": "β4 subunit"},
        {"partner": "CSNK2A2", "uniprot": "P19784", "role": "Casein kinase 2, predicted phosphorylation (Group-Based Prediction)"},
    ],
    "KCNQ1": [
        {"partner": "KCNE1", "uniprot": "P15382", "role": "MinK, IKs modulation"},
        {"partner": "Calmodulin", "uniprot": "P0DP23", "role": "Calcium-dependent regulation"},
    ],
    "VCL": [
        {"partner": "Talin", "uniprot": "Q9Y490", "role": "Focal adhesion, actin linkage"},
        {"partner": "Actin", "uniprot": "P68133", "role": "F-actin binding, cytoskeleton"},
        {"partner": "α-Catenin", "uniprot": "P35221", "role": "Cadherin adhesion complex"},
        {"partner": "VASP", "uniprot": "P50552", "role": "Vasodilator-stimulated phosphoprotein, actin dynamics"},
    ],
}

# Default interactors when gene not in tissue-specific map
DEFAULT_INTERACTORS = {
    "Calmodulin": "P0DP23",
    "FGF13": "Q92913",
}

# PDB complexes for PPI structure-based predictions (gene_pair -> pdb_id)
PPI_PDB_COMPLEXES = {}

# Known PDB structures for key proteins
PROTEIN_PDB = {
    "SCN5A": ["6UZ3", "6UZ0", "8BEI"],  # Nav1.5; AlphaFold for disordered regions (e.g. ~1054)
    "KCNQ1": ["6UZZ"],
    "MYH7": ["4DB1"],
    "VCL": ["1TR2", "1SDR", "1SYQ"],  # 1TR2 preferred for linker region
}

# AlphaMissense API
ALPHAMISSENSE_API = "https://alphamissense.hegelab.org/hotspotapi"
UNIPROT_MAPPING_API = "https://rest.uniprot.org/uniprotkb/search"
RCSB_PDB_API = "https://data.rcsb.org/rest/v1/core/uniprot"
