# main.py
from pipeline import (setup_output_directory, fetch_sequence, calculate_sequence_stats,
                      save_summary_report, plot_gc_content, plot_base_composition, fetch_variants)

# === CONFIGURATION ===
ACCESSION = "CR536574.1"
EMAIL = "Gayathri.jonnalagaddaa@gmail.com"  # CHANGE THIS!
OUTPUT_DIR = r"C:\Users\jgaya\PyCharmMiscProject\tcf7l2_analysis_pipeline\output"


def main():
    print("Starting TCF7L2 Analysis Pipeline\n")

    # Use your fixed absolute path
    output_dir = setup_output_directory(OUTPUT_DIR)

    record = fetch_sequence(ACCESSION, EMAIL, output_dir)
    stats = calculate_sequence_stats(record)

    save_summary_report(stats, output_dir)
    plot_gc_content(record, output_dir)
    plot_base_composition(record, output_dir)

    print("\nPipeline completed successfully!")
    print(f"All results saved in: {output_dir}")

# Add variant analysis
    variants = fetch_variants("TCF7L2", output_dir)

if __name__ == "__main__":
    main()