#!/usr/bin/env python3
import os
import subprocess
import argparse

def run_command(command, description):
    """Run a shell command and print its description."""
    print(f"\n=== {description} ===")
    print(f"Executing: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        print("Command executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        exit(1)

def process_genome(fasta_url, gff_url, output_dir, strain, date):
    """Download, process, and load genome data for JBrowse."""
    # Step 1: Download and process reference genome
    fasta_file = f"GCA_{date}_genomic.fna.gz"
    fasta_unzipped = f"{strain}_{date}.fna"

    run_command(
        f"wget {fasta_url} -O {fasta_file}",
        "Downloading reference genome"
    )
    run_command(
        f"gunzip {fasta_file}",
        "Unzipping reference genome"
    )
    run_command(
        f"mv {fasta_file[:-3]} {fasta_unzipped}",
        "Renaming reference genome"
    )
    run_command(
        f"samtools faidx {fasta_unzipped}",
        "Indexing reference genome"
    )

    # Step 2: Load genome into JBrowse
    jbrowse_out = os.path.join(output_dir, "jbrowse2")
    run_command(
        f"jbrowse add-assembly {fasta_unzipped} --out {jbrowse_out} --load copy",
        "Loading genome into JBrowse"
    )

    # Step 3: Download and process genome annotations
    gff_file = f"GCA_{date}_genomic.gff.gz"
    gff_processed = f"{date}_genes.gff"
    gff_bgzip = f"{gff_processed}.gz"

    run_command(
        f"wget {gff_url} -O {gff_file}",
        "Downloading genome annotations"
    )
    run_command(
        f"gunzip {gff_file}",
        "Unzipping genome annotations"
    )
    run_command(
        f"jbrowse sort-gff {gff_file[:-3]} > {gff_processed}",
        "Sorting genome annotations"
    )
    run_command(
        f"bgzip {gff_processed}",
        "Compressing sorted annotations with bgzip"
    )
    run_command(
        f"tabix {gff_bgzip}",
        "Indexing compressed annotations with tabix"
    )

    # Step 4: Load annotation track into JBrowse
    run_command(
        f"jbrowse add-track {gff_bgzip} --out {jbrowse_out} --load copy --assemblyNames H3N2_{date}",
        "Loading annotation track into JBrowse"
    )

    # Step 5: Index for search-by-gene
    run_command(
        f"jbrowse text-index --out {jbrowse_out}",
        "Indexing for search-by-gene"
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and load genome data for JBrowse.")
    parser.add_argument('--fasta_url', required=True, help="URL of the reference genome (FASTA format).")
    parser.add_argument('--gff_url', required=True, help="URL of the genome annotations (GFF3 format).")
    parser.add_argument('--output_dir', default="/var/www/html", help="Output directory for JBrowse data.")
    parser.add_argument('--strain', required=True, help="strain of Influenza A virus (H1N1, H3N2).")
    parser.add_argument('--date', required=True, help="date identifier for the genome (e.g., 2015).")
    args = parser.parse_args()

    process_genome(args.fasta_url, args.gff_url, args.output_dir, args.strain, args.date)
