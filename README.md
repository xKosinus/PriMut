# PriMut - Primer Design, Protocol Generation, and Protein Variant Management
This repository provides a complete workflow for Site‑Directed Mutagenesis experiment planning — from primer design, to mutation step planning, to generating printable protocols and tube labels. It is tailored for molecular biology workflows involving multi‑site mutagenesis with automated primer design and PDF output.

✨ Key Features
- Automated Primer Design

  - Validates wildtype sequence start/stop codons and flanking regions.
  - Groups nearby mutations into single primer sets (adjustable distance threshold).
  - Uses primer3‑py to optimize melting temperatures and primer lengths.
  - Outputs TXT, CSV, and JSON primer files.

- Smart Stepwise Mutagenesis Planning
  - Calculates minimal mutation pathways, reusing intermediate variants where possible.
  - Supports configurable maximum mutations per step.
  - Maintains a persistent variant_databank.json to avoid duplicate work.

- Protocol PDF Generation
  - Rich formatting with overview tables, per-step details, and all variants generated.
  - Includes notice of pre-existing variants already in databank.

- Label Sheet PDF Generation
  - Automatically prints labels for all variants from the current run (final + intermediates).
  - Customizable start column/row to reuse partially used sheets.
  - Uses standard Hema 10900 label layout.

- Interactive Jupyter/Colab UI (optional)
  - Paste sequence and list of variants in a form.
  - Choose experiment parameters and run everything with a single button.

📂 Outputs
- primer_list.txt / .csv / .json – primer sequences and statistics.
- mutagenesis_protocol.pdf – full mutagenesis plan, with steps and overview tables.
- herma_10900_labels.pdf – printable tube labels.
- variant_databank.json – persistent variant registry.

🛠 Requirements
- Python 3.x
  - primer3-py (primer design)
  - reportlab (PDF generation)
  - (optional) ipywidgets, IPython.display for interactive use in Jupyter/Google Colab
