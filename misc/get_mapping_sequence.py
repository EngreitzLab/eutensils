import pandas as pd

def variant_id_to_mapping_sequence(variant_id, ref_sequence, ref_start_1based, ref_end_1based, padding=7):
    """
    Convert a VariantID to a MappingSequence (edited sequence with padding).
    
    Parameters:
    -----------
    variant_id : str
        Format: chrN:position:reference>alternate
        Position is 0-BASED genomic coordinate
    ref_sequence : str
        Reference sequence for the region
    ref_start_1based : int
        1-based start coordinate of the reference sequence (e.g., 10,673,890)
    ref_end_1based : int
        1-based end coordinate of the reference sequence (e.g., 10,674,170), inclusive
    padding : int
        Number of base pairs to include on each side of the edit
        
    Returns:
    --------
    str
        The edited sequence centered with padding
    """
    # Convert 1-based reference coordinates to 0-based for internal calculations
    ref_start_0based = ref_start_1based - 1
    
    # Parse the VariantID
    parts = variant_id.split(':')
    chrom = parts[0]
    position = int(parts[1])  # 0-based genomic coordinate
    edit_parts = parts[2].split('>')
    reference_allele = edit_parts[0]
    alternate_allele = edit_parts[1]
    
    # Convert genomic position to position within our reference sequence
    local_position = position - ref_start_0based
    
    # Validate that the reference allele matches
    ref_allele_in_seq = ref_sequence[local_position:local_position + len(reference_allele)]
    if ref_allele_in_seq != reference_allele:
        raise ValueError(f"Reference allele mismatch at position {position}. "
                        f"Expected '{reference_allele}', found '{ref_allele_in_seq}'")
    
    # Create the edited sequence
    edited_sequence = (ref_sequence[:local_position] + 
                       alternate_allele + 
                       ref_sequence[local_position + len(reference_allele):])
    
    # Calculate the position of the edit in the edited sequence
    edit_start = local_position
    edit_end = local_position + len(alternate_allele)
    
    # Extract mapping sequence with padding
    mapping_start = max(0, edit_start - padding)
    mapping_end = min(len(edited_sequence), edit_end + padding)
    
    mapping_sequence = edited_sequence[mapping_start:mapping_end]
    
    return mapping_sequence


def process_variant_file(input_file, ref_sequence, ref_start_1based, ref_end_1based, 
                         padding=7, output_file=None):
    """
    Process a file of variant IDs and return a table with mapping sequences.
    
    Parameters:
    -----------
    input_file : str
        Path to input file with one column containing variant IDs (with header)
    ref_sequence : str
        Reference sequence for the region
    ref_start_1based : int
        1-based start coordinate of the reference sequence
    ref_end_1based : int
        1-based end coordinate of the reference sequence
    padding : int
        Number of base pairs to include on each side of the edit
    output_file : str, optional
        Path to save output TSV file. If None, doesn't save.
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: VariantID, MappingSequence, Status
    """
    # Read the variant IDs from file
    df = pd.read_csv(input_file, sep='\t')  # Adjust separator if needed
    
    # Get the column name (first column)
    variant_col = df.columns[0]
    
    # Initialize results
    results = []
    
    # Process each variant
    for idx, variant_id in enumerate(df[variant_col]):
        try:
            mapping_seq = variant_id_to_mapping_sequence(
                variant_id, 
                ref_sequence, 
                ref_start_1based, 
                ref_end_1based, 
                padding=padding
            )
            results.append({
                'VariantID': variant_id,
                'MappingSequence': mapping_seq,
                'Status': 'Success'
            })
        except Exception as e:
            results.append({
                'VariantID': variant_id,
                'MappingSequence': None,
                'Status': f'Error: {str(e)}'
            })
            print(f"Warning: Failed to process variant {variant_id}: {e}")
    
    # Create output DataFrame
    output_df = pd.DataFrame(results)
    
    # Save to file if requested
    if output_file:
        output_df.to_csv(output_file, index=False, sep="\t")
        print(f"Results saved to {output_file}")
    
    return output_df


# Example usage:
if __name__ == "__main__":
    # Your reference sequence
    ref_seq = ("GGGGAGGGAGGGGAGAAAAAAAAAAACCAGCCTAGCTCGCGGGCCGGCCG"
               "CAGGTAACACAATGACGCGTGCCCGCCCGGCTCTCGGAGAAGGACCCGGA"
               "GAGCCCGTCTGGCAGCAGCGGCCGGGGCTGGCCACCTCTACCCAGCACGC"
               "CGGGCAGGGCGCATGCGCGCTTATTAATATTCATGAGAGGGCGTGCTCAC"
               "CCTGGGCACGCCCCTCCCCTTCACGTTGCTGGGGAGGGGGTAGTGCGAGG"
               "AGGAACTTGGAAGGGGTTGGGGGCAGCGGGA")
    
    ref_start = 10673890
    ref_end = 10674170
    
    # Process the file
    results_df = process_variant_file(
        input_file='JAG1_VariantIDs_hg38.txt',  # Replace with your file path
        ref_sequence=ref_seq,
        ref_start_1based=ref_start,
        ref_end_1based=ref_end,
        padding=7,  # Adjust padding as needed
        output_file='JAG1_VariantIDs_hg38_with_mapping.tsv'  # Output file path
    )
    
    # Display results
    print("\nResults:")
    print(results_df.to_string())
    
    # Show summary
    print(f"\n\nSummary:")
    print(f"Total variants processed: {len(results_df)}")
    print(f"Successful: {(results_df['Status'] == 'Success').sum()}")
    print(f"Failed: {(results_df['Status'] != 'Success').sum()}")