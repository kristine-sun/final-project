import requests
from Bio import SeqIO
from io import StringIO
import sys

# Get URL from the command-line arguments
url = sys.argv[1]

# Fetch the file
response = requests.get(url)
response.raise_for_status()

# Read the content into a StringIO object
gbff_content = StringIO(response.text)

# Convert GBFF to GFF3
output_file = "output.gff3"
with open(output_file, "w") as gff_file:
    for record in SeqIO.parse(gbff_content, "genbank"):
        for feature in record.features:
            if feature.type != "source":
                qualifiers = ";".join([f"{key}={','.join(value)}" for key, value in feature.qualifiers.items()])
                start = feature.location.start.position + 1
                end = feature.location.end.position
                strand = "+" if feature.strand == 1 else "-"
                gff_line = f"{record.id}\t.\t{feature.type}\t{start}\t{end}\t.\t{strand}\t.\t{qualifiers}\n"
                gff_file.write(gff_line)

print(f"GFF3 file saved to {output_file}")
