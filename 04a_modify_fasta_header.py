import re
import sys
from pathlib import Path

def modify_fasta_headers(input_file, output_file=None):
    """
    Modify FASTA headers by extracting content between 'cmp_ref=' and ';'
    and placing it at the start of the header with a space.
    
    Example: >A375_batch_1.1.1 gene_name=MSTRG.5;...;cmp_ref=MSTRG.5.1;...
    becomes: >MSTRG.5.1 A375_batch_1.1.1 gene_name=MSTRG.5;...;cmp_ref=MSTRG.5.1;...
    
    Args:
        input_file: Path to input FASTA file
        output_file: Path to output file (if None, modifies in place)
    """
    if output_file is None:
        output_file = input_file
    
    # Pattern to extract content between cmp_ref= and ;
    pattern = r'cmp_ref=([^;]+);'
    
    with open(input_file, 'r', encoding='utf-8') as infile:
        lines = infile.readlines()
    
    modified_lines = []
    for line in lines:
        if line.startswith('>'):  # Header line in FASTA
            # Search for cmp_ref= pattern
            match = re.search(pattern, line)
            if match:
                cmp_ref_value = match.group(1).strip()
                if cmp_ref_value:  # Only modify if there's actual content
                    # Remove '>' and add it back with the new format
                    header_content = line[1:].rstrip('\n')
                    modified_line = f">{cmp_ref_value} {header_content}\n"
                    modified_lines.append(modified_line)
                else:
                    modified_lines.append(line)
            else:
                modified_lines.append(line)
        else:
            modified_lines.append(line)
    
    with open(output_file, 'w', encoding='utf-8') as outfile:
        outfile.writelines(modified_lines)
    
    print(f"Modified FASTA file written to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_fasta> [output_fasta]")
        print("If output_fasta is not specified, the input file will be modified in place.")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    modify_fasta_headers(input_file, output_file)