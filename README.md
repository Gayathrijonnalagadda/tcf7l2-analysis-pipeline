# TCF7L2 Analysis Pipeline

A simple Python project to explore the **TCF7L2** gene, which plays a big role in how our body handles blood sugar and insulin.

## Why I Made This Project

I got interested in this because of continuous glucose monitoring (CGM) data and diabetes management. I read some papers about how genes can affect blood sugar responses after eating glucose. One key gene is TCF7L2 – certain changes (variants) in it make people more likely to get type 2 diabetes by reducing insulin release when glucose levels rise.

I wanted to learn how to use bioinformatics tools to look at real gene data myself. This project helped me practice Python, Biopython, and working with public databases like NCBI.

## What the Project Does

This is a small pipeline (script) that:
- Fetches the DNA sequence of TCF7L2 (using accession CR536574.1, a full coding sequence).
- Calculates basic info: sequence length, GC content (percentage of G and C bases), and shows the first 200 bases.
- Makes simple plots: GC content across the sequence and base composition (A, T, G, C counts).
- Searches for important genetic variants (SNPs) linked to disease, like rs7903146.
- Saves everything in an "output" folder: FASTA file, text report, plots, and variant list.

Run `python main.py` to see it in action!

## Scientific Question We Are Answering

**How do genetic variants in the TCF7L2 gene affect glucose metabolism and risk of type 2 diabetes?**

From the data and papers:
- TCF7L2 helps control insulin production in the pancreas.
- The common variant rs7903146 (T allele instead of C) is an "intron variant" – it doesn't change the protein directly but reduces how much TCF7L2 is made.
- People with the T allele have lower insulin release after eating glucose, leading to higher blood sugar and higher risk of type 2 diabetes.
- This shows why some people respond differently to the same diet – genetics matter for precision nutrition!

We learned real variants from dbSNP, including rs7903146 listed as a "risk factor".

## Tools Used

- **Python** with **Biopython** (great library for biology data from NCBI).
  - Biopython documentation: https://biopython.org/
- Matplotlib for plots.
- NCBI databases: Nucleotide and dbSNP.

## Key References (Papers That Inspired This)

These papers explain how genes like TCF7L2 and others affect glucose and insulin:

1. Lyssenko V, et al. (2008). Genetic variation in GIPR influences the glucose and insulin responses to an oral glucose challenge. *Nature Genetics*. https://pubmed.ncbi.nlm.nih.gov/20081857/

2. Mirzaei K, et al. (2014). Variants in glucose- and circadian rhythm-related genes affect the response of energy expenditure to weight-loss diets: the POUNDS LOST Trial. *American Journal of Clinical Nutrition*. https://pubmed.ncbi.nlm.nih.gov/24335056/

3. Villareal DT, et al. (2010). TCF7L2 variant rs7903146 affects the risk of type 2 diabetes by modulating incretin action. *Diabetes*. https://pubmed.ncbi.nlm.nih.gov/19934000/

More on TCF7L2 and rs7903146: Search "TCF7L2 rs7903146 diabetes" on PubMed for many studies.

## How to Run

1. Install requirements: `pip install biopython matplotlib`
2. Put your email in `main.py` (NCBI needs it).
3. Run `python main.py`

Outputs go to the "output" folder.

Feel free to fork and improve! This is my first bioinformatics pipeline – feedback welcome.

Created by Gayathri Jonnalagadda – December 2025

