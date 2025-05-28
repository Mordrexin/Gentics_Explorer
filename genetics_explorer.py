import gradio as gr
import pandas as pd
import matplotlib.pyplot as plt
import os
import io
import time # For timestamped filenames
import yaml # For loading YAML configuration
# import shutil # Could be used for more aggressive cleanup if needed

# --- 0. Configuration Loading ---
CONFIG_FILE = "genetics_explorer_config.yaml"
DEFAULT_CONFIG = {
    "app_title": "🌱 Plant Genetics Explorer (Default Config)",
    "app_description": "Default description. Create `genetics_explorer_config.yaml`.",
    "data_directory": "data_default",
    "database_options": ["NCBI (SRA, GenBank, RefSeq)", "Ensembl Plants"],
    "default_database": "NCBI (SRA, GenBank, RefSeq)",
    "mock_gene_data": { # Minimal fallback mock data
        "AT_DEFAULT": {
            "name": "Default Gene", "description": "Default description.",
            "traits": ["Default Trait"], "expression_levels": {"tissue1": 50},
            "go_terms": ["GO:0000000"], "ncbi_gene_id": "000000", "ensembl_plants_id": "AT_DEFAULT"
        }
    }
}

def load_config(config_path):
    """Loads configuration from a YAML file."""
    try:
        with open(config_path, 'r') as f:
            config_data = yaml.safe_load(f)
        print(f"Successfully loaded configuration from '{config_path}'")
        # Provide defaults for keys that might be missing to avoid KeyErrors later
        # This makes the config file more flexible if a user omits non-critical parts
        return {
            "app_title": config_data.get("app_title", DEFAULT_CONFIG["app_title"]),
            "app_description": config_data.get("app_description", DEFAULT_CONFIG["app_description"]),
            "data_directory": config_data.get("data_directory", DEFAULT_CONFIG["data_directory"]),
            "database_options": config_data.get("database_options", DEFAULT_CONFIG["database_options"]),
            "default_database": config_data.get("default_database", DEFAULT_CONFIG["default_database"]),
            "mock_gene_data": config_data.get("mock_gene_data", DEFAULT_CONFIG["mock_gene_data"])
        }
    except FileNotFoundError:
        print(f"Warning: Configuration file '{config_path}' not found. Using default configuration.")
        return DEFAULT_CONFIG
    except yaml.YAMLError as e:
        print(f"Error parsing YAML file '{config_path}': {e}. Using default configuration.")
        return DEFAULT_CONFIG
    except Exception as e:
        print(f"An unexpected error occurred while loading config: {e}. Using default configuration.")
        return DEFAULT_CONFIG

# Load configuration
config = load_config(CONFIG_FILE)

DATA_DIR = config["data_directory"]
MOCK_GENE_DATA = config["mock_gene_data"]

# Create a reverse mapping for traits to genes (dynamically from loaded MOCK_GENE_DATA)
MOCK_TRAIT_TO_GENES = {}
if MOCK_GENE_DATA: # Check if mock_gene_data is not empty
    for gene_id, data in MOCK_GENE_DATA.items():
        if "traits" in data and data["traits"]: # Check if 'traits' key exists and is not empty
            for trait in data["traits"]:
                if trait not in MOCK_TRAIT_TO_GENES:
                    MOCK_TRAIT_TO_GENES[trait] = []
                MOCK_TRAIT_TO_GENES[trait].append(gene_id)
else:
    print("Warning: MOCK_GENE_DATA is empty or not loaded correctly. Trait explorer might not function as expected.")


def setup_data_directory():
    """Creates the data directory if it doesn't exist."""
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Data directory '{DATA_DIR}' ensured at {os.path.abspath(DATA_DIR)}")

def list_data_files_md():
    """Returns a markdown string listing files in the DATA_DIR."""
    if not os.path.exists(DATA_DIR) or not os.listdir(DATA_DIR):
        return f"No files in data directory (`{DATA_DIR}`) yet."
    
    files = sorted(os.listdir(DATA_DIR)) # Sort for consistent order
    md_list = f"### Files in `./{DATA_DIR}` directory:\n"
    for f_name in files:
        md_list += f"- `{f_name}`\n"
    return md_list

# --- 1. Data Acquisition Guide Functions ---

def get_database_info(db_name, species, search_term):
    info = f"### Information for {db_name}\n"
    info += f"**Species:** {species if species else 'Not specified'}\n"
    info += f"**Search Term:** {search_term if search_term else 'Not specified'}\n\n"

    # Database specific information
    # Note: The actual logic here doesn't change much, but the available db_name options
    # can now be controlled by the YAML.
    if db_name == "NCBI (SRA, GenBank, RefSeq)":
        info += "**Access Points:**\n"
        info += "- **NCBI Search:** [https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/)\n"
        info += "- **SRA (Sequence Read Archive):** For raw sequencing data. Search by project, sample, or run accession.\n"
        info += "  - *Example Tool:* `sra-tools` (e.g., `prefetch SRRXXXXXX && fastq-dump SRRXXXXXX`)\n"
        info += "- **GenBank/RefSeq:** For assembled genomes, genes, proteins. Search by gene name, accession, or organism.\n"
        info += "  - *Example Tool:* `ncbi-datasets-cli` (e.g., `ncbi-datasets-cli download genome accession GCF_000001735.4` for Arabidopsis)\n"
        if species:
            info += f"  - *NCBI Taxonomy ID for {species} (example):* You'd look this up (e.g., Arabidopsis thaliana is 3702)\n"
        if search_term:
             info += f"  - *Direct NCBI Gene Search for '{search_term}':* [https://www.ncbi.nlm.nih.gov/gene/?term={search_term.replace(' ', '+')}%5Bgene%5D+AND+{species.replace(' ', '+')}%5Borgn%5D](https://www.ncbi.nlm.nih.gov/gene/?term={search_term.replace(' ', '+')}%5Bgene%5D+AND+{species.replace(' ', '+')}%5Borgn%5D)\n"

    elif db_name == "Ensembl Plants":
        info += "**Access Points:**\n"
        info += "- **Website:** [https://plants.ensembl.org/](https://plants.ensembl.org/)\n"
        info += "- **BioMart:** Data mining tool for custom dataset downloads.\n"
        info += "- **FTP Server:** For bulk downloads of genomes, annotations (GFF3), sequences (FASTA).\n"
        if species:
            species_formatted = species.lower().replace(" ", "_")
            info += f"  - *Example Species Page:* [https://plants.ensembl.org/{species_formatted}/Info/Index](https://plants.ensembl.org/{species_formatted}/Info/Index)\n"
        if search_term:
            info += f"  - *Example Search for '{search_term}' in {species if species else 'all species'}:* [https://plants.ensembl.org/Multi/Search/Results?species={species_formatted if species else 'all'};idx=;q={search_term}](https://plants.ensembl.org/Multi/Search/Results?species={species_formatted if species else 'all'};idx=;q={search_term})\n"

    elif db_name == "Phytozome (JGI)":
        info += "**Access Points:**\n"
        info += "- **Website:** [https://phytozome-next.jgi.doe.gov/](https://phytozome-next.jgi.doe.gov/)\n"
        info += "- **JGI Data Portal:** For browsing and downloading JGI-sequenced plant genomes.\n"
        if species:
             info += f"  - *Typically, you select species from their browser, then find download sections for FASTA, GFF3, etc.*\n"
    elif db_name == "The Arabidopsis Information Resource (TAIR)": # Example for new DB
        info += "**Access Points:**\n"
        info += "- **Website:** [https://www.arabidopsis.org/](https://www.arabidopsis.org/)\n"
        info += "- **Search:** Use the search bar for genes, loci, keywords.\n"
        info += "- **Bulk Data:** Often found in FTP or download sections.\n"
        if search_term:
            info += f"  - *Example Search for '{search_term}':* [https://www.arabidopsis.org/servlets/Search?type=general&search_action=detail&method=1&name={search_term.replace(' ', '+')}&sub_type=gene](https://www.arabidopsis.org/servlets/Search?type=general&search_action=detail&method=1&name={search_term.replace(' ', '+')}&sub_type=gene)\n"
    else:
        info += "Database not recognized or information not yet available for this configured option."

    # --- Create and save the sample FASTA file ---
    fasta_content = ""
    if MOCK_GENE_DATA:
        for gene_id, data in MOCK_GENE_DATA.items():
            fasta_content += f">{gene_id} {data.get('name', 'Unknown Name')}\n" # Use .get for safety
            # Dummy sequence
            fasta_content += "ATGCGTAGCATCGATCGATCGATCGTAGCATGCTAGCATCGATCGATCGTAGCTAGCATCG\n"
            fasta_content += "GCTAGCATCGATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCATCGATCGATCGTAG\n"
    else:
        fasta_content = ">NO_MOCK_DATA_AVAILABLE\nNNNNNNNNNN\n"
    
    safe_db_name = db_name.split(' ')[0].lower().replace('(', '').replace(')', '')
    timestamp = int(time.time())
    sample_filename = f"{safe_db_name}_sample_genes_{timestamp}.fasta"
    sample_file_path = os.path.join(DATA_DIR, sample_filename)
    
    file_for_download = None
    try:
        with open(sample_file_path, "w") as f:
            f.write(fasta_content)
        info += f"\n**Mock data file generated:** `{sample_filename}` has been saved to the `{DATA_DIR}` directory and is available for download."
        file_for_download = sample_file_path
    except Exception as e:
        info += f"\n**Error saving file:** Could not write to `{sample_file_path}`. Error: {e}"
    
    data_dir_listing_md = list_data_files_md()
        
    return gr.Markdown(info), file_for_download, gr.Markdown(data_dir_listing_md)

# --- 2. Gene-Trait Explorer Functions ---
def get_gene_info(gene_id):
    if not MOCK_GENE_DATA:
        return "Mock gene data not loaded. Cannot retrieve gene info.", None
    if not gene_id or gene_id not in MOCK_GENE_DATA:
        return "Please select a valid gene.", None

    data = MOCK_GENE_DATA[gene_id]
    info_md = f"### Gene: {gene_id} ({data.get('name', 'N/A')})\n"
    info_md += f"**Description:** {data.get('description', 'N/A')}\n"
    info_md += f"**Associated Traits:** {', '.join(data.get('traits', []))}\n"
    info_md += f"**NCBI Gene ID:** [{data.get('ncbi_gene_id', 'N/A')}](https://www.ncbi.nlm.nih.gov/gene/{data.get('ncbi_gene_id', '')})\n"
    info_md += f"**EnsemblPlants ID:** [{data.get('ensembl_plants_id', 'N/A')}](https://plants.ensembl.org/Arabidopsis_thaliana/Gene/Summary?g={data.get('ensembl_plants_id', '')}) (Assuming Arabidopsis for link)\n"
    info_md += "**Gene Ontology Terms:**\n"
    for go_term in data.get('go_terms', []):
        info_md += f"- [{go_term}](http://amigo.geneontology.org/amigo/term/{go_term})\n"

    plot_output = None
    if data.get('expression_levels'):
        tissues = list(data['expression_levels'].keys())
        levels = list(data['expression_levels'].values())
        
        fig, ax = plt.subplots()
        ax.bar(tissues, levels, color=['skyblue', 'lightgreen', 'lightcoral'])
        ax.set_ylabel('Mock Expression Level')
        ax.set_title(f'Mock Expression for {gene_id}')
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plot_output = fig
        
    return gr.Markdown(info_md), plot_output

def get_trait_info(trait_name):
    if not MOCK_TRAIT_TO_GENES:
         return "Trait data not initialized (likely due to missing mock gene data).", None
    if not trait_name or trait_name not in MOCK_TRAIT_TO_GENES:
        return "Please select a valid trait.", None

    associated_genes = MOCK_TRAIT_TO_GENES[trait_name]
    info_md = f"### Trait: {trait_name}\n"
    info_md += f"**Genes associated with this trait (mock data):**\n"
    gene_list_for_table = []
    for gene_id in associated_genes:
        gene_data = MOCK_GENE_DATA.get(gene_id, {}) # Use .get for safety
        info_md += f"- {gene_id} ({gene_data.get('name', 'N/A')})\n"
        gene_list_for_table.append({
            "Gene ID": gene_id,
            "Gene Name": gene_data.get('name', 'N/A'),
            "Description": gene_data.get('description', 'N/A')
        })
    df_genes = pd.DataFrame(gene_list_for_table)

    return gr.Markdown(info_md), df_genes


# --- 3. Gradio Interface ---
with gr.Blocks(theme=gr.themes.Soft()) as demo:
    gr.Markdown(config["app_title"]) # Use title from config
    gr.Markdown(config["app_description"]) # Use description from config

    with gr.Tab("Data Acquisition Guide"):
        gr.Markdown("## Select and 'Download' Open Genetic Repositories")
        gr.Markdown(
            "Choose a database and provide optional species/search terms to get information and links. "
            f"A sample FASTA file based on mock data will be **saved to the local `./{DATA_DIR}` directory** and offered for browser download."
        )
        with gr.Row():
            db_dropdown = gr.Dropdown(
                choices=config["database_options"], # Use choices from config
                label="Select Genetic Database",
                value=config["default_database"] # Use default from config
            )
        with gr.Row():
            species_input = gr.Textbox(label="Plant Species (e.g., Arabidopsis thaliana, Oryza sativa)", placeholder="Optional")
            search_term_input = gr.Textbox(label="Search Term (e.g., Gene ID, keyword like 'kinase')", placeholder="Optional")
        
        search_db_btn = gr.Button("Get Database Info & Generate Sample File", variant="primary")
        
        db_info_output = gr.Markdown(label="Database Information & Access Instructions")
        download_file_output = gr.File(label="Download Generated Sample File")
        
        gr.Markdown("---")
        data_dir_contents_output = gr.Markdown(value=list_data_files_md(), label=f"Files in ./{DATA_DIR}")
        
        search_db_btn.click(
            fn=get_database_info,
            inputs=[db_dropdown, species_input, search_term_input],
            outputs=[db_info_output, download_file_output, data_dir_contents_output]
        )

    with gr.Tab("Gene-Trait Explorer (Mock Data)"):
        gr.Markdown("## Evaluate How Genes Affect Plant Traits (Using Mock Data)")
        
        gene_ids_list = list(MOCK_GENE_DATA.keys()) if MOCK_GENE_DATA else []
        trait_list = list(MOCK_TRAIT_TO_GENES.keys()) if MOCK_TRAIT_TO_GENES else []

        with gr.Row():
            with gr.Column(scale=1):
                gr.Markdown("### Explore by Gene")
                gene_dropdown = gr.Dropdown(choices=gene_ids_list, label="Select Gene ID", 
                                            info="Populated from mock_gene_data in config.")
                gene_info_btn = gr.Button("Get Gene Info", variant="secondary")
            with gr.Column(scale=1):
                gr.Markdown("### Explore by Trait")
                trait_dropdown = gr.Dropdown(choices=trait_list, label="Select Plant Trait",
                                             info="Populated from mock_gene_data in config.")
                trait_info_btn = gr.Button("Get Trait Info", variant="secondary")

        with gr.Row():
            with gr.Column(scale=2):
                gene_details_output = gr.Markdown(label="Gene Details")
                trait_details_output = gr.Markdown(label="Trait Details")
                trait_genes_table_output = gr.DataFrame(label="Genes for Selected Trait")
            with gr.Column(scale=1):
                expression_plot_output = gr.Plot(label="Mock Gene Expression Plot")
        
        if gene_ids_list: # Only attach click handler if there are genes
            gene_info_btn.click(
                fn=get_gene_info,
                inputs=[gene_dropdown],
                outputs=[gene_details_output, expression_plot_output]
            ).then(lambda: (None, None), outputs=[trait_details_output, trait_genes_table_output])
        else:
            gr.Markdown("*Gene selection disabled: No mock gene data found or loaded.*")


        if trait_list: # Only attach click handler if there are traits
            trait_info_btn.click(
                fn=get_trait_info,
                inputs=[trait_dropdown],
                outputs=[trait_details_output, trait_genes_table_output]
            ).then(lambda: (None, None), outputs=[gene_details_output, expression_plot_output])
        else:
             gr.Markdown("*Trait selection disabled: No mock trait data found or loaded.*")


if __name__ == "__main__":
    setup_data_directory() # Ensure data directory (from config) exists on startup
    
    print(f"Application configured with DATA_DIR: '{DATA_DIR}'")
    print(f"Mock gene data keys loaded: {list(MOCK_GENE_DATA.keys()) if MOCK_GENE_DATA else 'None'}")
    print("Open application in a browser at http://127.0.0.1:7860 (or the address shown by Gradio)")
    demo.launch(debug=True) # debug=True is helpful for development
    
    print(f"Application closed. Files are stored in ./{DATA_DIR}")
    print(f"Full path to data directory: {os.path.abspath(DATA_DIR)}")