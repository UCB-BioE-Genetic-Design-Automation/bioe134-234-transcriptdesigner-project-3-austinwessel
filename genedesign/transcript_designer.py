import random
import csv
import sys
from genedesign.rbs_chooser import RBSChooser
from genedesign.models.transcript import Transcript
from genedesign.checkers.hairpin_checker import hairpin_checker
from genedesign.checkers.forbidden_sequence_checker import ForbiddenSequenceChecker
from genedesign.checkers.internal_promoter_checker import PromoterChecker
from genedesign.checkers.codon_checker import CodonChecker
from genedesign.checkers.hairpin_post_checker import HairpinPostChecker
random.seed(42)

class TranscriptDesigner:
    """
    Designs optimized DNA sequences for a given peptide using sliding windows
    and guided Monte Carlo sampling, ensuring that the sequences pass all the specified checks.
    """

    def __init__(self):
        self.amino_acid_to_codons = {}
        self.codon_to_amino_acid = {}
        self.codon_usage_freq = {}
        self.highest_cai_codons = {}
        self.codon_checker = CodonChecker()
        self.forbidden_checker = ForbiddenSequenceChecker()
        self.promoter_checker = PromoterChecker()
        self.rbsChooser = RBSChooser()
        self.hpost_checker = HairpinPostChecker()

    def initiate(self):
        # Initialize checkers
        self.codon_checker.initiate()
        self.forbidden_checker.initiate()
        self.promoter_checker.initiate()
        self.rbsChooser.initiate()
        
        # Load codon usage frequencies
        self.load_codon_usage()
        self.hpost_checker.initiate(c2aa=self.codon_to_amino_acid, aa2c=self.amino_acid_to_codons)

    def load_codon_usage(self):
        # Read and parse codon usage file
        with open('genedesign/data/codon_usage.txt', 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if len(row) < 3:
                    continue
                codon, aa, freq = row[0], row[1], float(row[2])
                self.codon_usage_freq.setdefault(aa, []).append((codon, freq))
                self.codon_to_amino_acid[codon] = aa

        # Normalize frequencies and identify highest CAI codon
        for aa, codons in self.codon_usage_freq.items():
            total = sum(freq for _, freq in codons)
            self.codon_usage_freq[aa] = [(codon, freq / total) for codon, freq in codons]
            # Sort codons by frequency to find the highest CAI codon
            sorted_codons = sorted(self.codon_usage_freq[aa], key=lambda x: x[1], reverse=True)
            self.highest_cai_codons[aa] = sorted_codons[0][0]

        # Build amino_acid_to_codons mapping
        self.amino_acid_to_codons = {
            aa: [codon for codon, _ in codons] for aa, codons in self.codon_usage_freq.items()
        }

    def sample_codon(self, amino_acid, next_codons=None, retries=100):
        """
        Samples a codon for the given amino acid based on usage frequency.
        """
        codons_freqs = self.codon_usage_freq[amino_acid]
        codons, freqs = zip(*codons_freqs)
        
        # attempted resampling to avoid hairpins, no improvement
        # for _ in range(retries):
        #     single_codon = random.choices(codons, weights=freqs, k=1)[0]
        #     if next_codons:
        #         test_seq = [single_codon] + next_codons[:2]
        #         if not self.creates_hairpin(test_seq):
        #             return single_codon

        #falll back 
        return random.choices(codons, weights=freqs, k=1)[0]
    
    def creates_hairpin(self, sequence):
        dna_seq = ''.join(sequence)
        hairpin_ok, _ = hairpin_checker(dna_seq)
        return  not hairpin_ok

    def optimize_sequence(self, peptide):
        """
        Optimize the DNA sequence corresponding to the given peptide.
        """
        # Initial translation using highest CAI codon
        initial_codons = [self.highest_cai_codons[aa] for aa in peptide]

        # Sliding window optimization
        optimized_codons = initial_codons[:1]
        window_size = 3  # Adjust based on performance and effectiveness
        context_size = 6  # Number of codons to consider after the window

        for i in range(1, len(initial_codons), 3):

            left_pane = optimized_codons[-3:] if len(optimized_codons) >= 3 else optimized_codons
            middle_codons = initial_codons[i:i + 3]
            right_pane = initial_codons[i + 3:i + 9]
            optimized_window = self.optimize_window(left_pane, middle_codons, right_pane)

            # Update the optimized_codons list
            optimized_codons.extend(optimized_window)
            optimized_codons.append('TAA')

        return optimized_codons
    

    def optimize_window(self, left_pane, middle_codons, right_pane):
        """
        Optimize a window of codons, considering context before and after.
        """
        best_codons = middle_codons.copy() 
        num_candidates = 6  # Increase to improve chances of finding valid sequences
        failed_samples = 0

        for _ in range(num_candidates):
            candidate_codons = middle_codons.copy() + right_pane.copy()
            
            # Sample new codons for the positions to optimize
            for idx in range(len(candidate_codons)):
                aa = self.codon_to_amino_acid[candidate_codons[idx]]
                next_codons = right_pane[:2] if idx < len(right_pane) else []
                candidate_codons[idx] = self.sample_codon(aa, next_codons)

            # combine static left pane with dynamic middle + right pane
            window = left_pane + candidate_codons

            # Perform checks, if valid exit early---including an early exit yielded no overall improvements, kept for processing speed
            if self.all_checks(window):
                return candidate_codons[:3]
            
            failed_samples += 1

        print(f"\n FAILEed {failed_samples}/{num_candidates} samples")

        return best_codons

    def all_checks(self, codons):

        dna_seq = ''.join(codons) #([bp for codon in codons for bp in codon])

        codons_ok, _, _, _ = self.codon_checker.run(codons)
        forbidden_ok, _ = self.forbidden_checker.run(dna_seq)
        promoter_ok, _ = self.promoter_checker.run(dna_seq)
        hairpin_ok, _ = hairpin_checker(dna_seq)

        if not (hairpin_ok and codons_ok and forbidden_ok and promoter_ok):
            return False

        return True

    def run(self, peptide, ignores):
        """
        Main method to generate the optimized DNA sequence.
        """
        # Optimize the sequence
        optimized_codons = self.optimize_sequence(peptide)  # this is sliding window without overlap 

        # HairpinPostChecker call and hairpin utils
        optimized_codons = self.hpost_checker.replace_hairpins(optimized_codons)

        # Validate final sequence
        cds = ''.join(optimized_codons)
        if len(cds) % 3 != 0:
            raise ValueError(f"Final CDS length is not divisible by 3: {len(cds)} nucleotides")

        # Final checks 
        hairpin_ok, hp_found = hairpin_checker(cds)                                              #Paused for DBUG
        if not hairpin_ok:
            print(f"\n hairpin at ", {hp_found})
            raise ValueError("Failed final check. Unable to generate sequence without hairpins.")

        promoter_ok, _ = self.promoter_checker.run(cds)
        if not promoter_ok:
            raise ValueError("Failed final check. Unable to generate sequence without internal promoters.")

        forbidden_ok, _ = self.forbidden_checker.run(cds)
        if not forbidden_ok:
            raise ValueError("Failed final check. Sequence contains forbidden sequences.")

        # Return the optimized codons
        # Select an RBS

        selectedRBS = self.rbsChooser.run(cds, ignores)
        # Return the Transcript object
        return Transcript(selectedRBS, peptide, optimized_codons)

# Example usage
if __name__ == "__main__":
    peptide = "MYPFIRTARMTV"
    #peptide = "MYPFIRTARMT" # this produced MYPFIRTARFI
    #peptide = "MYPFIRTARM" #this produced MYPFIRTARF
    #peptide = "MYPFIRTAR"

    designer = TranscriptDesigner()
    designer.initiate()
    ignores = set()
    transcript = designer.run(peptide, ignores)

    #Print out the transcript information

    #print(f"input: {peptide}")

    #test_codon_translate = ''.join(transcript.codons)

    # test_output = designer.translate_dna_to_protein(test_codon_translate)

    # print(f"\noutput:{test_output}")

    #print(transcript.codons)