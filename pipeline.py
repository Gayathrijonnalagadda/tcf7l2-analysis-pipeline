# pipeline.py
from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
import os
import xml.etree.ElementTree as ET


def setup_output_directory(output_dir):
    """Create output directory if it doesn't exist."""
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory ready: {output_dir}")
    return output_dir


def fetch_sequence(accession, email, output_dir):
    """Fetch sequence from NCBI and save as FASTA."""
    Entrez.email = email

    print(f"Fetching sequence for {accession}...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()

    fasta_path = os.path.join(output_dir, "tcf7l2_sequence.fasta")
    with open(fasta_path, "w") as f:
        f.write(f">{record.description}\n")
        for i in range(0, len(record.seq), 80):
            f.write(str(record.seq[i:i + 80]) + "\n")

    print(f"Sequence saved to {fasta_path}")
    return record


def calculate_sequence_stats(record):
    """Calculate length, GC content, and first 200 bases."""
    sequence = record.seq
    length = len(sequence)
    gc_content = (sequence.count("G") + sequence.count("C")) / length * 100
    first_200 = str(sequence[:200])

    stats = {
        "accession": record.id,
        "description": record.description,
        "length": length,
        "gc_content": gc_content,
        "first_200": first_200
    }
    return stats


def save_summary_report(stats, output_dir):
    """Save a nice summary text report."""
    report_path = os.path.join(output_dir, "tcf7l2_summary.txt")
    with open(report_path, "w") as f:
        f.write("TCF7L2 Sequence Analysis Report\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Accession: {stats['accession']}\n")
        f.write(f"Description: {stats['description']}\n\n")
        f.write(f"Sequence Length: {stats['length']} bases\n")
        f.write(f"GC Content: {stats['gc_content']:.2f}%\n\n")
        f.write("First 200 bases:\n")
        f.write(stats['first_200'] + "\n\n")
        f.write("Biological Note:\n")
        f.write(
            "TCF7L2 is a key transcription factor in the Wnt pathway and strongly associated with type 2 diabetes risk.\n")
        f.write("Common variant rs7903146 (T allele) reduces insulin secretion in response to glucose.\n")

    print(f"Summary report saved to {report_path}")


def plot_gc_content(record, output_dir):
    """Plot sliding window GC content."""
    sequence = record.seq
    window_size = 200
    gc_values = []
    positions = []

    for i in range(0, len(sequence) - window_size + 1, window_size):
        window = sequence[i:i + window_size]
        gc = (window.count("G") + window.count("C")) / window_size * 100
        gc_values.append(gc)
        positions.append(i + window_size / 2)

    plt.figure(figsize=(10, 6))
    plt.plot(positions, gc_values, label="GC Content", color="blue", marker="o")
    plt.title("Sliding Window GC Content in TCF7L2 (CR536574.1)")
    plt.xlabel("Position (bases)")
    plt.ylabel("GC %")
    plt.legend()
    plt.grid(True)
    plot_path = os.path.join(output_dir, "tcf7l2_gc_plot.png")
    plt.savefig(plot_path)
    plt.close()
    print(f"GC plot saved to {plot_path}")


def plot_base_composition(record, output_dir):
    """Plot A/T/G/C base counts."""
    sequence = record.seq
    counts = {'A': sequence.count('A'), 'T': sequence.count('T'),
              'G': sequence.count('G'), 'C': sequence.count('C')}

    plt.figure(figsize=(6, 4))
    plt.bar(counts.keys(), counts.values(), color=['green', 'red', 'blue', 'purple'])
    plt.title("Base Composition in TCF7L2 Sequence")
    plt.xlabel("Base")
    plt.ylabel("Count")
    plot_path = os.path.join(output_dir, "tcf7l2_base_composition.png")
    plt.savefig(plot_path)
    plt.close()
    print(f"Base composition plot saved to {plot_path}")


def fetch_variants(gene_name="TCF7L2", output_dir="."):
    """Fetch clinically significant variants for TCF7L2 from dbSNP (with namespace handling)."""
    print(f"Searching dbSNP for clinically significant variants in {gene_name}...")

    search_term = f"{gene_name}[Gene Name] AND human[Organism] AND (pathogenic[Clinical Significance] OR risk factor[Clinical Significance])"
    handle = Entrez.esearch(db="snp", term=search_term, retmax=20)
    search_results = Entrez.read(handle)
    handle.close()

    variant_ids = search_results["IdList"]
    print("Found variant IDs:", variant_ids)
    print("Including key rs7903146?", "7903146" in variant_ids)

    if "7903146" not in variant_ids:
        variant_ids.append("7903146")

    variants = []

    # Namespace for dbSNP XML
    ns = {'ds': 'https://www.ncbi.nlm.nih.gov/SNP/docsum'}

    for rs_id in variant_ids:
        print(f"Fetching data for rs{rs_id}...")
        fetch_handle = Entrez.efetch(db="snp", id=rs_id, retmode="xml")
        xml_data = fetch_handle.read().decode('utf-8')
        fetch_handle.close()

        try:
            root = ET.fromstring(xml_data)

            for doc in root.findall(".//ds:DocumentSummary", namespaces=ns):
                rsid_elem = doc.find("ds:SNP_ID", ns)
                rsid = rsid_elem.text if rsid_elem is not None else "N/A"

                alleles_set = set()
                for freq in doc.findall(".//ds:FREQ", ns):
                    text = freq.text or ""
                    if "=" in text:
                        base = text.split("=")[0].strip().upper()
                        if base in "ACGT":
                            alleles_set.add(base)
                alleles = "/".join(sorted(alleles_set)) if alleles_set else "N/A"

                # Fallback for known rs7903146 (reference C, variant T)
                if rsid == "7903146" and len(alleles_set) == 1:
                    alleles = "C/T"

                fxn = doc.find("ds:FXN_CLASS", ns)
                fxn_class = fxn.text.strip() if fxn is not None and fxn.text else "N/A"

                clin = doc.find("ds:CLINICAL_SIGNIFICANCE", ns)
                clin_sig = clin.text.strip() if clin is not None and clin.text else "N/A"

                variant_dict = {
                    "rsid": f"rs{rsid}",
                    "alleles": alleles,
                    "function": fxn_class,
                    "clinical_significance": clin_sig
                }
                variants.append(variant_dict)

                print(f"Parsed: rs{rsid} | Alleles: {alleles} | Function: {fxn_class} | Significance: {clin_sig}")

        except ET.ParseError as e:
            print(f"Parse error for rs{rs_id}: {e}")

    if variants:
        variant_path = os.path.join(output_dir, "tcf7l2_variants.txt")
        with open(variant_path, "w") as f:
            f.write("TCF7L2 Clinically Significant Variants (from dbSNP)\n")
            f.write("=" * 70 + "\n")
            f.write(f"{'RSID':<12} {'Alleles':<15} {'Function Class':<25} {'Clinical Significance'}\n")
            f.write("-" * 100 + "\n")
            for v in variants:
                f.write(f"{v['rsid']:<12} {v['alleles']:<15} {v['function']:<25} {v['clinical_significance']}\n")

        print(f"\nVariant report saved to {variant_path} with {len(variants)} entries!")
    else:
        print("\nNo variants parsed — check if tags exist in full XML.")

    return variants