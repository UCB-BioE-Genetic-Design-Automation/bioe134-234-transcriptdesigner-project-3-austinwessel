"""
test_transcript_designer.py:

Purpose: Tests the TranscriptDesigner class to ensure mRNA sequences are correctly generated from protein sequences.
Coverage:
Test the reverse-translation of protein sequences into optimized DNA.
Ensure proper RBS assignment for each mRNA.
Validate handling of codon optimization and RNA folding requirements.
"""

import os
import traceback
import csv
import time
from statistics import mean
from genedesign.seq_utils.Translate import Translate
from genedesign.transcript_designer import TranscriptDesigner
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.codon_checker import CodonChecker


## helper functions
def extract_fasta_subset(input_file, output_file, num_sequences=10):
    """
    Extracts a subset of sequences from a large FASTA file and writes them to a new file.

    Parameters:
        input_file (str): Path to the input FASTA file.
        output_file (str): Path to the output FASTA file to save the subset.
        num_sequences (int): Number of sequences to extract.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        count = 0
        write_sequence = False
        for line in infile:
            if line.startswith(">"):  # Header line
                if count >= num_sequences:
                    break
                count += 1
                write_sequence = True
                outfile.write(line)
            elif write_sequence:  # Sequence line
                outfile.write(line)

    print(f"Subset of {num_sequences} sequences written to {output_file}")

def parse_fasta(fasta_file):
    """
    Parses the FASTA file to extract gene names and protein sequences.
    """
    sequences = {}
    current_gene = None
    current_sequence = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_gene:
                    sequences[current_gene] = ''.join(current_sequence)
                gene_name = None
                parts = line.split()
                for part in parts:
                    if part.startswith("GN="):
                        gene_name = part.split("=")[1]
                        break
                if not gene_name:
                    gene_name = line.split('|')[2].split(' ')[0]
                current_gene = gene_name
                current_sequence = []
            else:
                current_sequence.append(line)
        if current_gene:
            sequences[current_gene] = ''.join(current_sequence)
    
    return sequences

# Example usage
#input_fasta = "large_input.fasta"  # Replace with your input FASTA file path
#output_fasta = "subset_output.fasta"  # Replace with your desired output file path
#extract_fasta_subset(input_fasta, output_fasta, num_sequences=100)







## TranscriptDesigner test for subset of .fasta file 

def benchmark_proteome(fasta_file):
    """
    Benchmarks the proteome using TranscriptDesigner.
    """
    designer = TranscriptDesigner()
    designer.initiate()
    extract_fasta_subset(fasta_file,"tests/benchmarking/first_100.fasta", 10)
    proteome = parse_fasta("tests/benchmarking/first_100.fasta")
    successful_results = []
    error_results = []

    for gene, protein in proteome.items():
        try:
            print(f"Processing gene: {gene} with protein sequence: {protein[:30]}...")
            ignores = set()
            transcript = designer.run(protein, ignores)
            successful_results.append({
                'gene': gene,
                'protein': protein,
                'transcript': transcript
            })
        except Exception as e:
            error_results.append({
                'gene': gene,
                'protein': protein,
                'error': f"Error: {str(e)}\nTraceback: {traceback.format_exc()}"
            })
    
    return successful_results, error_results

def analyze_errors(error_results):
    """
    Write the error analysis to a text file.
    """
    error_summary = {}
    with open('first_100_error_summary.txt', 'w') as f:
        for error in error_results:
            error_message = error['error'].split("\n")[0]
            error_summary[error_message] = error_summary.get(error_message, 0) + 1
            f.write(f"Gene: {error['gene']}\n{error['error']}\n\n")
    
    return error_summary

def validate_transcripts(successful_results):
    """
    Validate the successful transcripts using various checkers, now including CodonChecker.
    """
    forbidden_checker = ForbiddenSequenceChecker()
    forbidden_checker.initiate()
    promoter_checker = PromoterChecker()
    promoter_checker.initiate()
    translator = Translate()
    translator.initiate()
    codon_checker = CodonChecker()  # Initialize CodonChecker
    codon_checker.initiate()  # Load the codon usage data

    validation_failures = []
    for result in successful_results:
        cds = ''.join(result['transcript'].codons)
        try:
            # Check if CDS length is a multiple of 3
            if len(cds) % 3 != 0:
                raise ValueError("CDS length is not a multiple of 3.")

            # Verify that the translated protein matches the original protein
            original_protein = result['protein']
            translated_protein = translator.run(cds)
            if original_protein != translated_protein:
                raise ValueError(f"Translation mismatch: Original {original_protein}, Translated {translated_protein}")

            # Ensure CDS starts with valid start codon and ends with stop codon
            if not (cds.startswith(("ATG", "GTG", "TTG")) and cds.endswith(("TAA", "TGA", "TAG"))):
                raise ValueError("CDS does not start with a valid start codon or end with a valid stop codon.")
        except ValueError as e:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': cds,
                'site': f"Translation or completeness error: {str(e)}"
            })
            continue

        # Validate against hairpins, forbidden sequences, and internal promoters
        transcript_dna = result['transcript'].rbs.utr.upper() + cds
        passed_hairpin, hairpin_string = hairpin_checker(transcript_dna)
        if not passed_hairpin:
            formatted_hairpin = hairpin_string.replace('\n', ' ').replace('"', "'")
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Hairpin detected: {formatted_hairpin}"
            })

        passed_forbidden, forbidden_site = forbidden_checker.run(transcript_dna)
        if not passed_forbidden:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Forbidden sequence: {forbidden_site}"
            })

        passed_promoter, found_promoter = promoter_checker.run(transcript_dna)
        if not passed_promoter:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': transcript_dna,
                'site': f"Constitutive promoter detected: {found_promoter}" if found_promoter else "Constitutive promoter detected"
            })

        codons_above_board, codon_diversity, rare_codon_count, cai_value = codon_checker.run(result['transcript'].codons)
        if not codons_above_board:
            validation_failures.append({
                'gene': result['gene'],
                'protein': result['protein'],
                'cds': cds,
                'site': f"Codon usage check failed: Diversity={codon_diversity}, Rare Codons={rare_codon_count}, CAI={cai_value}"
            })
    
    return validation_failures

def write_validation_report(validation_failures):
    """
    Writes validation results to a TSV file.
    """
    with open('first_10_validation_failures.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['gene', 'protein', 'cds', 'site'])
        for failure in validation_failures:
            writer.writerow([failure['gene'], failure['protein'], failure['cds'], failure['site']])

def generate_summary(total_genes, parsing_time, execution_time, errors_summary, validation_failures):
    """
    Generates a streamlined summary report categorizing validation failures by checker.
    """
    total_validation_failures = len(validation_failures)
    
    # Categorize failures by checker type
    checker_failures = {
        'Forbidden Sequence Checker': 0,
        'Hairpin Checker': 0,
        'Codon Usage Checker': 0,
        'Promoter Checker': 0,
        'Translation/Completeness Checker': 0
    }

    # Increment the appropriate checker category based on the failure site
    for failure in validation_failures:
        site = failure['site']
        if "Forbidden sequence" in site:
            checker_failures['Forbidden Sequence Checker'] += 1
        elif "Hairpin detected" in site:
            checker_failures['Hairpin Checker'] += 1
        elif "Codon usage check failed" in site:
            checker_failures['Codon Usage Checker'] += 1
        elif "Constitutive promoter detected" in site:
            checker_failures['Promoter Checker'] += 1
        elif "Translation or completeness error" in site:
            checker_failures['Translation/Completeness Checker'] += 1

    # Generate the summary report
    with open('first_10_summary_report.txt', 'w') as f:
        f.write(f"Total genes processed: {total_genes}\n")
        f.write(f"Parsing runtime: {parsing_time:.2f} seconds\n")
        f.write(f"Execution runtime: {execution_time:.2f} seconds\n")
        f.write(f"Total exceptions: {sum(errors_summary.values())}\n")

        if errors_summary:
            f.write(f"\nTop 3 most common exceptions:\n")
            for error, count in sorted(errors_summary.items(), key=lambda x: x[1], reverse=True)[:3]:
                f.write(f"- {error}: {count} occurrences\n")
        else:
            f.write("No exceptions encountered.\n")

        f.write(f"\nTotal validation failures: {total_validation_failures}\n")

        # Categorize validation failures by checker
        f.write("\nValidation Failures by Checker:\n")
        for checker, count in checker_failures.items():
            f.write(f"- {checker}: {count} occurrences\n")

def run_benchmark(fasta_file):
    """
    Runs the complete benchmark process: parsing, running TranscriptDesigner, validating, and generating reports.
    """
    start_time = time.time()
    
    # Benchmark the proteome
    parsing_start = time.time()
    successful_results, error_results = benchmark_proteome(fasta_file)
    parsing_time = time.time() - parsing_start
    
    # Analyze and log errors
    errors_summary = analyze_errors(error_results)
    
    # Validate the successful transcripts
    validation_start = time.time()
    validation_failures = validate_transcripts(successful_results)
    execution_time = time.time() - validation_start

    # Write validation and error reports
    write_validation_report(validation_failures)
    
    # Generate the summary report
    total_genes = len(successful_results) + len(error_results)
    generate_summary(total_genes, parsing_time, execution_time, errors_summary, validation_failures)

if __name__ == "__main__":
    fasta_file = "tests/benchmarking/uniprotkb_proteome_UP000054015_2024_09_24.fasta"
    run_benchmark(fasta_file)