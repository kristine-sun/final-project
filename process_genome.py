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
        print(f"Error while running: {command}")
        print(f"Return code: {e.returncode}")
        print(f"Output: {e.output}")
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

    # Step 2: Check if assembly already exists
    jbrowse_out = os.path.join(output_dir, "jbrowse2")
    assembly_name = f"{strain}_{date}"

    # Check if the assembly already exists
    check_assembly_command = f"jbrowse list-assemblies --out {jbrowse_out}"
    result = subprocess.run(check_assembly_command, shell=True, capture_output=True, text=True)

    if assembly_name in result.stdout:
        print(f"Assembly {assembly_name} already exists. Skipping addition...")
    else:
        # Step 3: Load genome into JBrowse if the assembly does not exist
        run_command(
            f"jbrowse add-assembly {fasta_unzipped} --out {jbrowse_out} --load copy",
            "Loading genome into JBrowse"
        )

    # Step 4: Download and process genome annotations
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

    # Step 5: Load annotation track into JBrowse
    run_command(
        f"jbrowse add-track {gff_bgzip} --out {jbrowse_out} --load copy --assemblyNames {assembly_name}",
        "Loading annotation track into JBrowse"
    )

    # Step 6: Index for search-by-gene
    run_command(
        f"jbrowse text-index --out {jbrowse_out}",
        "Indexing for search-by-gene"
    )

# def process_genome2(fasta_url, gff_url, output_dir, strain, date):
#     """Download, process, and load genome data for JBrowse."""
#     # Step 1: Download and process reference genome
#     fasta_file = f"GCA_{date}_genomic.fna.gz"
#     fasta_unzipped = f"{strain}_{date}.fna"

#     run_command(
#         f"wget {fasta_url} -O {fasta_file}",
#         "Downloading reference genome"
#     )
#     run_command(
#         f"gunzip {fasta_file}",
#         "Unzipping reference genome"
#     )
#     run_command(
#         f"mv {fasta_file[:-3]} {fasta_unzipped}",
#         "Renaming reference genome"
#     )
#     run_command(
#         f"samtools faidx {fasta_unzipped}",
#         "Indexing reference genome"
#     )

#     # Step 2: Load genome into the same assembly in JBrowse
#     jbrowse_out = os.path.join(output_dir, "jbrowse2")
    
#     # Use the same assembly name for all genomes (e.g., 'all_genomes_assembly')
#     run_command(
#         f"jbrowse add-assembly {fasta_unzipped} --out {jbrowse_out} --load copy --name all_genomes_assembly",
#         "Loading genome into the same assembly"
#     )

#     # Step 3: Download and process genome annotations
#     gff_file = f"GCA_{date}_genomic.gff.gz"
#     gff_processed = f"{date}_genes.gff"
#     gff_bgzip = f"{gff_processed}.gz"

#     run_command(
#         f"wget {gff_url} -O {gff_file}",
#         "Downloading genome annotations"
#     )
#     run_command(
#         f"gunzip {gff_file}",
#         "Unzipping genome annotations"
#     )
#     run_command(
#         f"jbrowse sort-gff {gff_file[:-3]} > {gff_processed}",
#         "Sorting genome annotations"
#     )
#     run_command(
#         f"bgzip {gff_processed}",
#         "Compressing sorted annotations with bgzip"
#     )
#     run_command(
#         f"tabix {gff_bgzip}",
#         "Indexing compressed annotations with tabix"
#     )

#     # Step 4: Load annotation track into the same assembly in JBrowse
#     run_command(
#         f"jbrowse add-track {gff_bgzip} --out {jbrowse_out} --load copy --assemblyNames all_genomes_assembly",
#         "Loading annotation track into the same assembly"
#     )

#     # Step 5: Index for search-by-gene
#     run_command(
#         f"jbrowse text-index --out {jbrowse_out}",
#         "Indexing for search-by-gene"
#     )

def parse_args():
    parser = argparse.ArgumentParser(description="Process genomic data.")
    parser.add_argument("--fasta_url", required=True, help="URL to the FASTA file")
    parser.add_argument("--gff_url", required=True, help="URL to the GFF file")
    parser.add_argument("--output_dir", default="/var/www/html/jbrowse2", help="Directory for output files")
    parser.add_argument("--date", required=True, help="Date of the genome data (e.g., 2019_04_30)")
    parser.add_argument("--strain", required=True, help="Strain of the virus (e.g., H1N1, H3N2)")
    return parser.parse_args()

def main():
    args = parse_args()
    print(f"Processing genome data for strain {args.strain} from date {args.date}.")
    
    # Call the process_genome2 function with the parsed arguments
    process_genome2(
        fasta_url=args.fasta_url,
        gff_url=args.gff_url,
        output_dir=args.output_dir,
        strain=args.strain,
        date=args.date
    )
    
if __name__ == "__main__":
     main()
