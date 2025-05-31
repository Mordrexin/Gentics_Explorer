

import gradio as gr
import pandas as pd 
import matplotlib.pyplot as plt 
import os
import io
import time 
import yaml 
import urllib.parse 

from Bio import Entrez, SeqIO
from xml.etree import ElementTree 

# --- 0. Configuration Loading ---
CONFIG_FILE = "genetics_explorer_config.yaml"
DEFAULT_CONFIG = {
    "app_title": "🌱 Plant Genetics Explorer (Live NCBI Data - Default)",
    "app_description": "Explore live gene data from NCBI. Configure `genetics_explorer_config.yaml`.",
    "data_directory": "genetic_data_live_default",
    "entrez_email": "your.email@example.com", 
    "ncbi_plant_species_options": ["Arabidopsis thaliana", "Oryza sativa Japonica Group"],
}

def load_config(config_path):
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config_data = yaml.safe_load(f)
        print(f"Successfully loaded configuration from '{config_path}'")
        return {
            "app_title": config_data.get("app_title", DEFAULT_CONFIG["app_title"]),
            "app_description": config_data.get("app_description", DEFAULT_CONFIG["app_description"]),
            "data_directory": config_data.get("data_directory", DEFAULT_CONFIG["data_directory"]),
            "entrez_email": config_data.get("entrez_email", DEFAULT_CONFIG["entrez_email"]),
            "ncbi_plant_species_options": config_data.get("ncbi_plant_species_options", DEFAULT_CONFIG["ncbi_plant_species_options"]),
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

config = load_config(CONFIG_FILE)

DATA_DIR = config["data_directory"]
ENTREZ_EMAIL = config["entrez_email"]
NCBI_PLANT_SPECIES = config["ncbi_plant_species_options"]

if ENTREZ_EMAIL == "your.email@example.com" or not ENTREZ_EMAIL:
    print("CRITICAL WARNING: NCBI Entrez email is not set or is default. Update in YAML for reliable NCBI access.")
Entrez.email = ENTREZ_EMAIL

app_state = {
    "last_fetched_ncbi_gene_details": None
}

def setup_data_directory():
    os.makedirs(DATA_DIR, exist_ok=True)
    print(f"Data directory '{DATA_DIR}' ensured at {os.path.abspath(DATA_DIR)}")

def list_data_files_md():
    if not os.path.exists(DATA_DIR) or not os.listdir(DATA_DIR):
        return f"No files in data directory (`{DATA_DIR}`) yet."
    files = sorted(os.listdir(DATA_DIR)) 
    md_list = f"### Files in `./{DATA_DIR}` directory:\n"
    for f_name in files:
        md_list += f"- `{f_name}`\n"
    return md_list

# --- 2. Gene Explorer Function (Solely for Live NCBI Data) ---
def display_fetched_ncbi_gene_in_explorer():
    global app_state 
    if app_state.get("last_fetched_ncbi_gene_details"):
        gene_data = app_state["last_fetched_ncbi_gene_details"]
        if not gene_data:
             return gr.Markdown("Error: Last fetched NCBI gene details are empty. Please re-fetch from NCBI tab."), gr.Markdown("")
    else:
        return gr.Markdown("No NCBI Gene data has been successfully fetched yet. Use the 'NCBI Gene/Protein Search' tab to fetch a **Gene record**."), gr.Markdown("")
    info_md = f"### Gene Information (Source: Live NCBI Data)\n"
    info_md += f"**Input ID/Symbol for Fetch:** {gene_data.get('input_id', 'N/A')}\n" # Display original input
    info_md += f"**Resolved NCBI GeneID:** {gene_data.get('ncbi_gene_id', 'N/A')}\n"
    info_md += f"**Official Symbol:** {gene_data.get('name', 'N/A')}\n"
    info_md += f"**Official Full Name:** {gene_data.get('full_name_from_summary', gene_data.get('description', 'N/A'))}\n" # Use full name, fallback to desc
    info_md += f"**Summary/Description:** {gene_data.get('description', 'N/A')}\n"

    if 'organism' in gene_data:
        info_md += f"**Organism:** {gene_data.get('organism', 'N/A')}\n"
    
    ncbi_gene_id_val = gene_data.get('ncbi_gene_id', None) 
    if ncbi_gene_id_val and ncbi_gene_id_val != "UNKNOWN": 
        info_md += f"**NCBI Gene Page Link:** [{ncbi_gene_id_val}](https://www.ncbi.nlm.nih.gov/gene/{ncbi_gene_id_val})\n"
    
    go_terms = gene_data.get('go_terms', [])
    if go_terms:
        info_md += "**Gene Ontology (GO) Terms:**\n"
        for go_term_entry in go_terms:
            if isinstance(go_term_entry, str) and ';' in go_term_entry:
                term_id, term_desc_cat = go_term_entry.split(';', 1)
                term_id = term_id.strip()
                term_desc_cat = term_desc_cat.strip()
                info_md += f"- [{term_id}](http://amigo.geneontology.org/amigo/term/{term_id}) {term_desc_cat}\n"
            elif isinstance(go_term_entry, str): 
                info_md += f"- [{go_term_entry}](http://amigo.geneontology.org/amigo/term/{go_term_entry})\n"
    else:
        info_md += "**Gene Ontology (GO) Terms:** Not available or not fetched.\n"
    expression_info_md = "*Expression plot not applicable for NCBI Gene records via this tool.*"
    return gr.Markdown(info_md), gr.Markdown(expression_info_md)

# --- 3. NCBI Data Retrieval Functions ---
def fetch_go_terms_for_gene(numeric_gene_id_str):
    go_terms_list = []
    if not numeric_gene_id_str or not numeric_gene_id_str.isdigit():
        print(f"fetch_go_terms_for_gene: Invalid or non-numeric Gene ID '{numeric_gene_id_str}'. Cannot fetch GO terms.")
        return go_terms_list          
    try:
        print(f"Fetching GO terms for numeric Gene ID: {numeric_gene_id_str}")
        handle = Entrez.esummary(db="gene", id=numeric_gene_id_str, retmode="xml")
        tree = ElementTree.parse(handle)
        handle.close()
        root = tree.getroot()
        for doc_summary in root.findall(".//DocumentSummary"): # Should be only one for a single ID
            for go_info in doc_summary.findall(".//GeneOntologyInfo/GeneOntology"):
                go_id_elem = go_info.find("GOTermID")
                go_term_name_elem = go_info.find("GOTermName")
                go_category_elem = go_info.find("GOCategory")
                go_id = go_id_elem.text if go_id_elem is not None else None
                go_term_name = go_term_name_elem.text if go_term_name_elem is not None else None
                go_category = go_category_elem.text if go_category_elem is not None else None
                if go_id and go_term_name and go_category:
                    go_terms_list.append(f"{go_id}; {go_term_name} ({go_category})")
    except Exception as e:
        print(f"Error fetching/parsing GO terms for Gene ID {numeric_gene_id_str}: {e}")
    return go_terms_list

def fetch_ncbi_data_by_id(id_to_use_for_fetch, db_type, original_input_id=None):
    """
    Fetches NCBI data using a NUMERIC GeneID (for Gene type) or Accession (for Protein type).
    id_to_use_for_fetch: The ID ESearch resolved to (numeric for Gene) or direct input.
    original_input_id: The very first ID the user typed (symbol, LOC, numeric etc.) for context.
    """
    global app_state 
    app_state["last_fetched_ncbi_gene_details"] = None 
    if not id_to_use_for_fetch:
        return "Error: No ID provided to fetch_ncbi_data_by_id.", "", None 
    
    if original_input_id is None: # If not passed, use the fetch ID as original
        original_input_id = id_to_use_for_fetch

    info_md_content = ""
    fasta_str_content = ""
    download_file_path = None
    
    try:
        if db_type == "Gene":
            # This function now primarily expects id_to_use_for_fetch to be the NUMERIC GeneID
            if not id_to_use_for_fetch.isdigit():
                msg = f"Internal Error: fetch_ncbi_data_by_id for Gene type expects a numeric ID, but received '{id_to_use_for_fetch}'. ESearch resolution might be needed."
                print(msg)
                return f"**{msg}**", "", None

            numeric_gene_id = id_to_use_for_fetch
            print(f"Fetching details for NUMERIC Gene ID: {numeric_gene_id}")

            # 1. Fetch Gene Summary using ESummary
            official_symbol_from_summary = original_input_id # Default to original input if symbol not found in summary
            full_name_from_summary = "N/A"
            description_from_summary = "N/A"
            organism_from_summary = "N/A"
            
            print(f"Fetching ESummary for Gene ID: {numeric_gene_id}")
            handle_summary = Entrez.esummary(db="gene", id=numeric_gene_id, retmode="xml")
            summary_tree = ElementTree.parse(handle_summary)
            handle_summary.close()
            summary_root = summary_tree.getroot()
            doc_sum = summary_root.find(".//DocumentSummary") # Should only be one for a single ID
            if doc_sum is not None:
                name_elem = doc_sum.find("Name") 
                desc_elem = doc_sum.find("Summary") 
                if desc_elem is None: desc_elem = doc_sum.find("Description")
                org_elem = doc_sum.find("Organism/ScientificName")
                off_sym_elem = doc_sum.find("NomenclatureSymbol")

                if name_elem is not None and name_elem.text: full_name_from_summary = name_elem.text
                if desc_elem is not None and desc_elem.text: description_from_summary = desc_elem.text
                if org_elem is not None and org_elem.text: organism_from_summary = org_elem.text
                if off_sym_elem is not None and off_sym_elem.text: official_symbol_from_summary = off_sym_elem.text
            else:
                print(f"Could not find DocumentSummary in ESummary XML for {numeric_gene_id}")

            # 2. Fetch Sequence using EFetch with rettype="fasta_cds_na" or "fasta" (genomic)
            sequence_length = 0
            print(f"Fetching FASTA sequence for Gene ID: {numeric_gene_id}")
            try:
                handle_fasta_seq = Entrez.efetch(db="gene", id=numeric_gene_id, rettype="fasta_cds_na", retmode="text")
                fasta_records = list(SeqIO.parse(handle_fasta_seq, "fasta"))
                handle_fasta_seq.close()
                if fasta_records:
                    seq_record = fasta_records[0] 
                    fasta_str_content = f">{seq_record.id} {seq_record.description}\n{str(seq_record.seq)}\n"
                    sequence_length = len(seq_record.seq)
                    print(f"Fetched CDS FASTA: {seq_record.id}")
                else: 
                    print(f"No CDS FASTA found for {numeric_gene_id}. Trying genomic FASTA.")
                    handle_fasta_genomic = Entrez.efetch(db="gene", id=numeric_gene_id, rettype="fasta", retmode="text")
                    fasta_records_genomic = list(SeqIO.parse(handle_fasta_genomic, "fasta"))
                    handle_fasta_genomic.close()
                    if fasta_records_genomic:
                        seq_record = fasta_records_genomic[0]
                        fasta_str_content = f">{seq_record.id} {seq_record.description}\n{str(seq_record.seq)}\n"
                        sequence_length = len(seq_record.seq)
                        print(f"Fetched genomic FASTA: {seq_record.id}")
                    else:
                        print(f"No FASTA sequence (CDS or genomic) found for Gene ID {numeric_gene_id}.")
            except Exception as e_fasta:
                print(f"Error during FASTA sequence retrieval for {numeric_gene_id}: {e_fasta}")
                # fasta_str_content remains ""

            # 3. Fetch GO Terms
            go_terms = fetch_go_terms_for_gene(numeric_gene_id)
            
            fetched_details_for_explorer = {
                "input_id": original_input_id, # Store the original user input for context
                "id": numeric_gene_id, 
                "name": official_symbol_from_summary, 
                "full_name_from_summary": full_name_from_summary,
                "description": description_from_summary, 
                "organism": organism_from_summary,
                "ncbi_gene_id": numeric_gene_id, 
                "go_terms": go_terms,
                "sequence_length": sequence_length
            }
            app_state["last_fetched_ncbi_gene_details"] = fetched_details_for_explorer
            
            info_md_content += f"### NCBI Gene: {official_symbol_from_summary} (GeneID: {numeric_gene_id})\n"
            info_md_content += f"*(Successfully resolved from input: '{original_input_id}')*\n" if original_input_id != numeric_gene_id else ""
            info_md_content += f"**Official Full Name:** {full_name_from_summary}\n"
            info_md_content += f"**Summary/Description:** {description_from_summary}\n"
            # ... (rest of info_md_content)
            info_md_content += f"**Organism:** {organism_from_summary}\n"
            info_md_content += f"**Sequence Length:** {sequence_length} bp (first CDS/genomic record)\n"
            if go_terms:
                 info_md_content += f"**Fetched GO Terms:** {len(go_terms)} (see Gene Explorer tab for list).\n"
            else:
                 info_md_content += f"**GO Terms:** None found or error during fetch.\n"
            info_md_content += f"**View on NCBI:** [https://www.ncbi.nlm.nih.gov/gene/{numeric_gene_id}](https://www.ncbi.nlm.nih.gov/gene/{numeric_gene_id})\n"
            filename_id_part = official_symbol_from_summary if official_symbol_from_summary != "N/A" else numeric_gene_id

        elif db_type == "Protein": 
            # id_to_use_for_fetch is the protein accession here
            handle_fasta_prot = Entrez.efetch(db="protein", id=id_to_use_for_fetch, rettype="fasta", retmode="text")
            record_prot = SeqIO.read(handle_fasta_prot, "fasta")
            handle_fasta_prot.close()
            app_state["last_fetched_ncbi_gene_details"] = None 
            info_md_content += f"### NCBI Protein: {record_prot.id}\n"
            info_md_content += f"**Description:** {record_prot.description}\n"
            info_md_content += f"**Sequence Length:** {len(record_prot.seq)} aa\n"
            info_md_content += f"**View on NCBI:** [https://www.ncbi.nlm.nih.gov/protein/{id_to_use_for_fetch}](https://www.ncbi.nlm.nih.gov/protein/{id_to_use_for_fetch})\n"
            fasta_str_content = f">{record_prot.id} {record_prot.description}\n{str(record_prot.seq)}\n"
            filename_id_part = record_prot.id.split('|')[-2] if '|' in record_prot.id and len(record_prot.id.split('|')) > 2 else record_prot.id
            info_md_content += "\n*Note: Protein data loaded. The 'Gene Explorer' tab is for 'Gene' record details.*"
        else:
            return f"Error: Invalid database type '{db_type}'.", "", None

        if fasta_str_content:
            timestamp = int(time.time())
            safe_filename_id = "".join(c for c in str(filename_id_part) if c.isalnum() or c in ('_', '-')).strip() or "record"
            fasta_filename = f"ncbi_{db_type.lower()}_{safe_filename_id}_{timestamp}.fasta"
            temp_file_path = os.path.join(DATA_DIR, fasta_filename)
            with open(temp_file_path, "w") as f:
                f.write(fasta_str_content)
            download_file_path = temp_file_path
            info_md_content += f"\n**File saved:** `{fasta_filename}` in `{DATA_DIR}` and available for download."
        else:
            info_md_content += "\n*Note: No sequence data was retrieved for this record or an error occurred during sequence fetch.*"

        if db_type == "Gene" and app_state["last_fetched_ncbi_gene_details"]:
             info_md_content += "\n**Data ready to load in 'Gene Explorer' tab.**"
        
        return info_md_content, fasta_str_content, download_file_path

    except Entrez.HTTPError as http_err: 
        error_code = http_err.code
        error_reason = str(http_err) # http_err itself can be printed for full details
        error_msg_detail = f"Error fetching NCBI record '{id_to_use_for_fetch}' ({db_type}): HTTP Error {error_code} ({error_reason})."
        if error_code == 400: # Bad Request
             error_msg_detail += " This often means the identifier is ambiguous, not found for direct EFetch/ESummary, or there's an issue with the request parameters."
        print(f"Detailed HTTP error in fetch_ncbi_data_by_id for '{id_to_use_for_fetch}': {error_msg_detail}")
        app_state["last_fetched_ncbi_gene_details"] = None 
        return f"**{error_msg_detail}**\n\n*Please verify the ID and try again. If the problem persists, NCBI services might be temporarily unavailable or the ID is not directly retrievable with EFetch/ESummary.*", "", None
    except Exception as e: 
        error_msg_detail = f"Unexpected error processing NCBI record '{id_to_use_for_fetch}' ({db_type}): {str(e)}"
        print(f"Detailed error in fetch_ncbi_data_by_id for '{id_to_use_for_fetch}': {error_msg_detail}") 
        app_state["last_fetched_ncbi_gene_details"] = None 
        return f"**{error_msg_detail}**", "", None


def perform_ncbi_search(selected_species, keyword_identifier, db_type):
    global app_state
    app_state["last_fetched_ncbi_gene_details"] = None 
    
    if not Entrez.email or Entrez.email == "your.email@example.com":
        error_md = "CRITICAL: Entrez email not configured in YAML. NCBI access may fail."
        return gr.Markdown(""), gr.Textbox(value="", visible=False), gr.File(value=None, visible=False), gr.Markdown(error_md), gr.Markdown(list_data_files_md())
    if not keyword_identifier:
        error_md = "Please enter a search term/ID."
        return gr.Markdown(""), gr.Textbox(value="", visible=False), gr.File(value=None, visible=False), gr.Markdown(error_md), gr.Markdown(list_data_files_md())

    original_input_id = keyword_identifier.strip() # Store the original input for context
    info_str_md_content = ""
    fasta_str_content = ""
    download_file_obj = None 
    error_md_content = ""
    
    is_numeric_gene_id = db_type == "Gene" and original_input_id.isdigit()
    is_loc_tag = db_type == "Gene" and original_input_id.upper().startswith("LOC") and original_input_id[3:].isdigit()
    is_protein_accession = db_type == "Protein" and \
                           (original_input_id.startswith(("NP_", "XP_", "WP_", "YP_", "AP_")) or \
                            (original_input_id.count('.') == 1 and len(original_input_id) > 5 and original_input_id.replace('.', '').isalnum()) or \
                            (len(original_input_id) >= 6 and original_input_id.isalnum() and not original_input_id.isdigit() and original_input_id.count('_') == 0) )
    is_potential_gene_symbol = db_type == "Gene" and not is_numeric_gene_id and not is_loc_tag and \
                               original_input_id.isalnum() and not (' ' in original_input_id or ',' in original_input_id or ';' in original_input_id)

    if is_numeric_gene_id or is_loc_tag or is_protein_accession:
        # For LOC tags or direct numeric GeneIDs, or protein accessions, pass them directly to fetch_ncbi_data_by_id
        # fetch_ncbi_data_by_id will handle them (expecting numeric for Gene or accession for Protein)
        id_for_fetch = original_input_id
        if is_loc_tag and not id_for_fetch.isdigit(): # For LOC tag, fetch_ncbi_data_by_id expects the full LOCxxxx string
             pass # No change needed, fetch_ncbi_data_by_id handles LOC tags
        elif is_numeric_gene_id:
             pass # Already numeric
        
        print(f"Attempting direct fetch for ID: {id_for_fetch}, DB: {db_type}")
        info_str_md_content, fasta_str_content, download_file_path = fetch_ncbi_data_by_id(id_for_fetch, db_type, original_input_id=original_input_id)
        if download_file_path:
             download_file_obj = download_file_path 
        if not fasta_str_content and not download_file_path and ("Error:" in info_str_md_content or "**Error" in info_str_md_content):
            error_md_content = info_str_md_content 
            info_str_md_content = "" 
    
    elif is_potential_gene_symbol:
        print(f"Input '{original_input_id}' identified as potential gene symbol. Using ESearch with species context '{selected_species}'.")
        search_term_for_esearch = f"{original_input_id}[Gene Name]" # Search by Gene Name field for symbols
        if selected_species and selected_species != "All Plant Species":
            search_term_for_esearch += f" AND \"{selected_species}\"[Organism]" # Enclose species in quotes for ESearch
        
        try:
            print(f"ESearch term: {search_term_for_esearch}")
            handle_esearch = Entrez.esearch(db="gene", term=search_term_for_esearch, retmax="5", sort="relevance")
            results = Entrez.read(handle_esearch)
            handle_esearch.close()
            id_list = results["IdList"]

            if len(id_list) == 1:
                numeric_gene_id_from_esearch = id_list[0]
                info_str_md_content = f"Found unique Gene ID '{numeric_gene_id_from_esearch}' for symbol '{original_input_id}' (context: {selected_species}). Fetching details...\n\n"
                details_md, fasta_str_content, download_file_path = fetch_ncbi_data_by_id(numeric_gene_id_from_esearch, "Gene", original_input_id=original_input_id)
                info_str_md_content += details_md # Append fetched details to the ESearch result message
                if download_file_path:
                    download_file_obj = download_file_path
                if not fasta_str_content and not download_file_path and ("Error:" in info_str_md_content or "**Error" in info_str_md_content):
                     # Error from fetching the resolved ID, info_str_md_content already contains the error from fetch_ncbi_data_by_id
                     # We might want to make the initial part of info_str_md_content conditional on success
                     if "Error:" in details_md or "**Error" in details_md: # Check if fetching the resolved ID failed
                          info_str_md_content = f"Found unique Gene ID '{numeric_gene_id_from_esearch}' for symbol '{original_input_id}' but an error occurred fetching its details. {details_md}"
                          error_md_content = "" # Error is already in info_str_md_content
                     else: # Success
                          pass # info_str_md_content is already good.

            elif len(id_list) > 1:
                summaries_md = ""
                # Optionally fetch summaries for ambiguous results
                # try:
                #     handle_sum_ambiguous = Entrez.esummary(db="gene", id=",".join(id_list[:3]), retmode="xml") # Summaries for top 3
                #     sum_tree_amb = ElementTree.parse(handle_sum_ambiguous)
                #     handle_sum_ambiguous.close()
                #     for doc_sum_amb in sum_tree_amb.findall(".//DocumentSummary"):
                #         sum_name = doc_sum_amb.find("Name").text if doc_sum_amb.find("Name") is not None else "N/A"
                #         sum_desc = doc_sum_amb.find("Summary").text if doc_sum_amb.find("Summary") is not None else "N/A"
                #         sum_id = doc_sum_amb.get("uid")
                #         summaries_md += f"- **ID:** {sum_id}, **Name:** {sum_name}, **Summary:** {sum_desc[:100]}...\n"
                # except Exception as e_sum_amb:
                #     summaries_md = f"\nCould not fetch summaries for ambiguous results: {e_sum_amb}"


                info_str_md_content = (f"Symbol '{original_input_id}' (context: {selected_species}) is ambiguous and matches multiple Gene IDs: **{', '.join(id_list)}**. "
                                       f"Please try one of these numeric IDs directly for a specific record.{summaries_md}")
                fasta_str_content = ""
                download_file_obj = None
            else: # len(id_list) == 0
                info_str_md_content = f"No direct NCBI Gene ID found for symbol '{original_input_id}' with species context '{selected_species}' via ESearch. It might be a broader keyword or the symbol/species combination is not found."
                # Fall through to generate keyword link
                fasta_str_content = ""
                download_file_obj = None
                # Generate keyword link because ESearch failed to find specific symbol
                query_parts_kw = [f"({original_input_id}[All Fields])"] 
                if selected_species and selected_species != "All Plant Species":
                    query_parts_kw.append(f"\"{selected_species}\"[Organism]")
                full_query_kw = " AND ".join(query_parts_kw)
                search_url_db_part_kw = "gene" 
                encoded_query_kw = urllib.parse.quote_plus(full_query_kw)
                search_url_kw = f"https://www.ncbi.nlm.nih.gov/{search_url_db_part_kw}/?term={encoded_query_kw}"
                info_str_md_content += (f"\n\nYou can try a broader search on NCBI: "
                                        f"**[Search NCBI for '{original_input_id}' in '{selected_species}']({search_url_kw})**")


        except Exception as e_esearch:
            error_md_content = f"Error during ESearch for symbol '{original_input_id}': {str(e_esearch)}"
            print(error_md_content)
            # Also provide a general keyword search link as a fallback
            query_parts_kw_fail = [f"({original_input_id}[All Fields])"]
            if selected_species and selected_species != "All Plant Species":
                query_parts_kw_fail.append(f"\"{selected_species}\"[Organism]")
            full_query_kw_fail = " AND ".join(query_parts_kw_fail)
            search_url_db_part_kw_fail = "gene"
            encoded_query_kw_fail = urllib.parse.quote_plus(full_query_kw_fail)
            search_url_kw_fail = f"https://www.ncbi.nlm.nih.gov/{search_url_db_part_kw_fail}/?term={encoded_query_kw_fail}"
            info_str_md_content = error_md_content + (f"\n\nAs ESearch failed, you can try a broader search on NCBI: "
                                        f"**[Search NCBI for '{original_input_id}' in '{selected_species}']({search_url_kw_fail})**")
            fasta_str_content = ""
            download_file_obj = None
    else: 
        # General keyword search if not a direct ID and not a potential symbol
        print(f"Treating input '{original_input_id}' as general keyword for DB: {db_type}, Species: {selected_species}")
        query_parts_gen = [f"({original_input_id}[All Fields])"] 
        if selected_species and selected_species != "All Plant Species":
            query_parts_gen.append(f"\"{selected_species}\"[Organism]")
        
        full_query_gen = " AND ".join(query_parts_gen)
        search_url_db_part_gen = "gene" if db_type == "Gene" else "protein"
        encoded_query_gen = urllib.parse.quote_plus(full_query_gen)
        search_url_gen = f"https://www.ncbi.nlm.nih.gov/{search_url_db_part_gen}/?term={encoded_query_gen}"
        info_str_md_content = (f"**Keyword-based Search Information:**\n\n"
                               f"Input '{original_input_id}' was treated as a general keyword. "
                               f"A link to the NCBI website is provided below.\n\n"
                               f"- **Database:** {db_type}\n"
                               f"- **Species Context:** {selected_species if selected_species and selected_species != 'All Plant Species' else 'Any/Not specified'}\n"
                               f"- **Search Term:** {original_input_id}\n"
                               f"- **Constructed NCBI Query String:** `{full_query_gen}`\n\n"
                               f"**➡️ [Click here to perform this search directly on the NCBI website]({search_url_gen})**")
        fasta_str_content = "" 
        download_file_obj = None 
            
    current_file_list_md = list_data_files_md()
    show_sequence_and_download = bool(fasta_str_content)

    # Consolidate error display
    if error_md_content and not info_str_md_content: # If only error, show it as primary info
        info_str_md_content = error_md_content
        error_md_content = "" # Clear separate error display if main info is error
    elif error_md_content: # Append to info if both exist
        info_str_md_content += f"\n\n**Additional Errors/Status:**\n{error_md_content}"


    return gr.Markdown(info_str_md_content), \
           gr.Textbox(value=fasta_str_content, label="Fetched Sequence (FASTA Format)", lines=10, interactive=False, visible=show_sequence_and_download), \
           gr.File(value=download_file_obj, label="Download Fetched FASTA File", visible=show_sequence_and_download), \
           gr.Markdown(""), \
           gr.Markdown(current_file_list_md) # Error display is now part of main info if needed


# --- 4. Gradio Interface ---
with gr.Blocks(theme=gr.themes.Soft()) as demo:
    gr.Markdown(f"# {config['app_title']}")
    gr.Markdown(config["app_description"])
    shared_data_dir_contents_output = gr.Markdown(value=list_data_files_md(), label=f"Files in ./{DATA_DIR}")

    # Tab for "Data Acquisition Guide" is REMOVED

    with gr.Tab("NCBI Gene/Protein Search (Live Data)"):
        gr.Markdown("## Fetch Gene or Protein Data Directly from NCBI")
        gr.Markdown("Enter an NCBI Gene ID (e.g., `816394`), Locus Tag (e.g. `LOC542363`), Gene Symbol (e.g. `PHYB`), or Protein Accession (e.g., `NP_171631.1`). If a **Gene ID/Symbol/Locus Tag** is provided and resolved to a unique Gene record, its details will be fetched for the 'Gene Explorer' tab. For broader keywords or ambiguous symbols, a search link to NCBI may be generated.")
        if not Entrez.email or Entrez.email == "your.email@example.com":
             gr.Markdown("<p style='color:red;font-weight:bold;'>WARNING: NCBI Entrez email is not configured! Update in `genetics_explorer_config.yaml`. NCBI lookups may fail.</p>")
        with gr.Row():
            species_dropdown_ncbi = gr.Dropdown(
                choices=["All Plant Species"] + NCBI_PLANT_SPECIES, value="All Plant Species",
                label="Select Plant Species (Context for Symbol/Keyword Search)",
                info="Helps disambiguate gene symbols. Less critical for direct numeric IDs/Accessions."
            )
            ncbi_db_type_dropdown = gr.Dropdown(choices=["Gene", "Protein"], value="Gene", 
                                                label="Select NCBI Database Type (Choose 'Gene' to populate Gene Explorer)")
        keyword_identifier_input_ncbi = gr.Textbox(
            label="Gene Symbol, NCBI Gene ID, Locus Tag, Protein Acc, or Keyword", 
            placeholder="e.g., PHYB, 816394, LOC542363, NP_171631.1"
        )
        fetch_ncbi_btn = gr.Button("Search/Fetch from NCBI", variant="primary")
        ncbi_info_output = gr.Markdown(label="Search Results / Fetched Information") # Main output for info and errors
        ncbi_sequence_output = gr.Textbox(label="Fetched Sequence (FASTA Format)", lines=10, interactive=False, visible=False)
        ncbi_download_file_output = gr.File(label="Download Fetched FASTA File", visible=False)
        # ncbi_error_output = gr.Markdown(label="Status/Errors") # Errors now part of ncbi_info_output
        gr.Markdown("---")
        gr.Markdown(render=shared_data_dir_contents_output)
        fetch_ncbi_btn.click(
            fn=perform_ncbi_search,
            inputs=[species_dropdown_ncbi, keyword_identifier_input_ncbi, ncbi_db_type_dropdown],
            outputs=[ncbi_info_output, ncbi_sequence_output, ncbi_download_file_output, gr.Markdown(""), shared_data_dir_contents_output] # Pass empty md to removed error output
        )
        
    with gr.Tab("Gene Explorer (Live NCBI Data)"):
        gr.Markdown("## Explore Details of Fetched NCBI Gene Record")
        gr.Markdown("After successfully fetching a **Gene record** from the 'NCBI Gene/Protein Search' tab (often by providing a numeric GeneID, or a symbol that resolves uniquely), click the button below to load and display its details, including GO terms.")
        with gr.Row():
            load_fetched_gene_btn = gr.Button("Load Last Fetched NCBI Gene Details", variant="primary", scale=1)
        with gr.Row():
            with gr.Column(scale=2):
                gene_details_output_explorer = gr.Markdown(label="Fetched Gene Details") 
            with gr.Column(scale=1):
                expression_info_explorer = gr.Markdown(label="Expression Info") 
        load_fetched_gene_btn.click(
            fn=display_fetched_ncbi_gene_in_explorer,
            inputs=[], 
            outputs=[gene_details_output_explorer, expression_info_explorer] 
        )
        if not app_state.get("last_fetched_ncbi_gene_details"):
            gr.Markdown("*No NCBI Gene data loaded yet. Use the 'NCBI Gene/Protein Search' tab first to fetch a 'Gene' record by its ID/Symbol.*", elem_id="initial_gene_explorer_message")

if __name__ == "__main__":
    setup_data_directory()
    print(f"Application configured with DATA_DIR: '{DATA_DIR}'")
    print(f"NCBI Entrez Email set to: '{Entrez.email}' (CRITICAL: Update in config if this is the default!)")
    print(f"Configured Plant species for NCBI search: {NCBI_PLANT_SPECIES}")
    print("Application is now focused on live NCBI data. Mock data and Guide tab have been removed.")
    print("Open application in a browser at http://127.0.0.1:7860 (or the address shown by Gradio)")
    demo.launch(debug=True, show_error=True) 
    print(f"Application closed. Files are stored in ./{DATA_DIR}")
    print(f"Full path to data directory: {os.path.abspath(DATA_DIR)}")

