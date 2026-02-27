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
# partner = symbol for keys (e.g. SCN5A_CALM1); label = display name in UI
CARDIAC_MYOCYTE_INTERACTORS = {
    "SCN5A": [
        {"partner": "FGF13", "label": "FGF13", "uniprot": "Q92913", "role": "Late Na+ current regulation, FHF binding"},
        {"partner": "CALM1", "label": "Calmodulin", "uniprot": "P0DP23", "role": "IQ domain binding, inactivation"},
        {"partner": "SCN1B", "label": "SCN1B", "uniprot": "Q07699", "role": "β1 subunit, trafficking"},
        {"partner": "SCN2B", "label": "SCN2B", "uniprot": "Q9UQ35", "role": "β2 subunit"},
        {"partner": "SCN3B", "label": "SCN3B", "uniprot": "Q9NY72", "role": "β3 subunit"},
        {"partner": "SCN4B", "label": "SCN4B", "uniprot": "Q9H0B8", "role": "β4 subunit"},
    ],
    "KCNQ1": [
        {"partner": "KCNE1", "label": "KCNE1", "uniprot": "P15382", "role": "MinK, IKs modulation"},
        {"partner": "CALM1", "label": "Calmodulin", "uniprot": "P0DP23", "role": "Calcium-dependent regulation"},
    ],
    "VCL": [
        {"partner": "TLN1", "label": "Talin", "uniprot": "Q9Y490", "role": "Focal adhesion, actin linkage"},
        {"partner": "ACTB", "label": "Actin", "uniprot": "P68133", "role": "F-actin binding, cytoskeleton"},
        {"partner": "CTNNA1", "label": "α-Catenin", "uniprot": "P35221", "role": "Cadherin adhesion complex"},
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
    "SCN5A": ["6UZ3", "6UZ0", "8BEI"],  # Nav1.5 cardiac sodium channel
    "KCNQ1": ["6UZZ"],
    "MYH7": ["4DB1"],
    "VCL": ["1TR2", "1SDR", "1SYQ"],
}

EXAMPLE_VARIANTS = [
    {"gene": "SCN5A", "mutation": "c.1577G>A, p.R526H", "label": "SCN5A p.R526H"},
    {"gene": "SCN5A", "mutation": "c.3160T>G, p.Ser1054Ala", "label": "SCN5A p.Ser1054Ala"},
    {"gene": "VCL", "mutation": "c.2507A>G, p.Gln836Arg", "label": "VCL p.Gln836Arg"},
]

# AlphaMissense API
ALPHAMISSENSE_API = "https://alphamissense.hegelab.org/hotspotapi"
UNIPROT_MAPPING_API = "https://rest.uniprot.org/uniprotkb/search"
RCSB_PDB_API = "https://data.rcsb.org/rest/v1/core/uniprot"
