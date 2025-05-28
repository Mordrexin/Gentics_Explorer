# 🌱 Plant Genetics Explorer

A Gradio-based web application for exploring plant genetic information. This tool currently serves as a demonstrator, providing a user-friendly interface to:

1.  **Data Acquisition Guide:** Offers information on accessing major plant genetic databases (NCBI, Ensembl Plants, Phytozome, TAIR - configurable) and generates a sample mock FASTA file.
2.  **Gene-Trait Explorer:** Allows exploration of mock gene data, displaying gene details, associated traits, and mock expression level plots.

The application is designed to be configurable via a `genetics_explorer_config.yaml` file.

## Features

*   **Interactive Web Interface:** Built with Gradio for ease of use.
*   **Configurable:**
    *   Application title and description.
    *   Data output directory.
    *   List of genetic databases and the default selection.
    *   Mock gene data for demonstration.
*   **Data Acquisition Guidance:** Provides links and example usage for selected genetic databases.
*   **Sample Data Generation:** Creates a mock FASTA file based on configured gene data and saves it to a local directory, also offering it for download.
*   **Gene Information Display:** Shows details for selected mock genes, including description, associated traits, external DB links, and GO terms.
*   **Trait Information Display:** Lists mock genes associated with selected plant traits.
*   **Expression Visualization:** Generates a simple bar plot for mock gene expression levels.
*   **Local File Management:** Lists files created in the configured data directory.

## Setup and Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/plant-genetics-explorer.git
    cd plant-genetics-explorer
    ```

2.  **Create a Python virtual environment (recommended):**
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows: venv\Scripts\activate
    ```

3.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```
    *(`requirements.txt` should initially contain:*
    ```
    gradio
    pandas
    matplotlib
    PyYAML
    biopython
    ```
    *)*

4.  **Configure the application:**
    *   Copy or rename `genetics_explorer_config.example.yaml` to `genetics_explorer_config.yaml`.
    *   Edit `genetics_explorer_config.yaml` to customize settings.

5.  **Run the application:**
    ```bash
    python app.py
    ```
    *(Assuming your main script is named `app.py`)*

6.  Open your web browser and navigate to the URL provided by Gradio (usually `http://127.0.0.1:7860`).

## Configuration (`genetics_explorer_config.yaml`)

The application's behavior is controlled by `genetics_explorer_config.yaml`. Key sections include:

*   `app_title`, `app_description`: Customize the UI text.
*   `data_directory`: Specifies where generated files (like sample FASTA) are stored.
*   `database_options`, `default_database`: Define the list of genetic databases shown in the "Data Acquisition Guide" tab.
*   `mock_gene_data`: A dictionary defining the mock genes used for demonstration. This will be phased out as real data retrieval is implemented.


## Initial Implementation Roadmap (TODO)

The following steps outline the initial development plan to integrate real data capabilities using Biopython and prepare for AlphaFold integration.

1.  **Integrate Biopython for NCBI Data Retrieval:**
    *   Modify the "Data Acquisition Guide" or create a new tab/section.
    *   Implement functionality to allow users to input an NCBI Gene ID or Protein Accession.
    *   Use `Bio.Entrez` to fetch gene/protein records (e.g., GenBank format).
    *   Parse the fetched record using `Bio.SeqIO` to extract relevant information (sequence, description, features, organism).
    *   Display the extracted information in the Gradio interface (e.g., sequence in a textbox, metadata as Markdown).
    *   Allow downloading the fetched sequence as a FASTA file.

2.  **Basic Sequence Analysis with Biopython:**
    *   Once a sequence is fetched (from step 1) or uploaded by the user (see step 3):
        *   Display basic sequence properties (length, GC content) using Biopython.
        *   Implement simple translation (DNA to Protein) if applicable.

3.  **User FASTA File Upload and Parsing:**
    *   Add a Gradio `gr.File` component to allow users to upload their own FASTA files.
    *   Use `Bio.SeqIO` to parse the uploaded FASTA file.
    *   If multiple sequences are present, allow the user to select one for display/analysis or display information for all.
    *   Connect this to the basic sequence analysis functionalities (step 2).

4.  **Integrate AlphaFold Database API for Structure Retrieval:**
    *   Add a section (likely near gene/protein info display) for protein structure.
    *   Allow users to input a UniProt ID (which can be obtained from NCBI records or user input).
    *   Implement a function to query the AlphaFold Protein Structure Database API (e.g., `https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}`) to get the PDB/mmCIF file URL.
    *   Provide a direct download link for the predicted structure file.

5.  **Basic 3D Structure Visualization:**
    *   Integrate a simple client-side 3D molecule viewer (like `py3Dmol` if it can be embedded well, or provide instructions/links for users to view the downloaded PDB file in external viewers like ChimeraX, PyMOL, or web-based viewers like Mol\*).
    *   If using an embedded viewer, load the PDB/mmCIF file fetched in step 4.
    *   *Note: Full local AlphaFold prediction is a more advanced step and deferred for now.*

6.  **Refactor Mock Data Usage:**
    *   Gradually replace reliance on `MOCK_GENE_DATA` in the "Gene-Trait Explorer" with data fetched via Biopython or from user uploads.
    *   The mock data can remain as a fallback or for demonstration if external services are unavailable.

7.  **Modularize Code:**
    *   As new functionalities are added, start organizing the Python code into separate modules (e.g., `ncbi_utils.py`, `alphafold_utils.py`, `sequence_analysis.py`) to keep `app.py` cleaner and more manageable.

