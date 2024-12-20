import itertools
import random
random.seed(42)
class HairpinPostChecker:
    """
    Another checker to post check for the formation of hairpins during codon optimization, 
    after running the provided checks during sampling, it alsoidentifies, removes, and 
    replaces hairpin-causing codons while preserving protein identity.
    """

    def __init__(self):
        """
        Initializes the HairpinPostChecker with empty parameters. Use `initiate` to set them.
        """
        self.hairpin_threshold = None
        self.c2aa = {}
        self.aa2c = {}

    def initiate(self, hairpin_threshold=6, c2aa=None, aa2c=None):
        """
        Initializes the HairpinPostChecker with thresholds and mappings.

        Parameters:
            hairpin_threshold (int): Minimum length of a complementary sequence to form a hairpin. Default is 6.
            codon_to_amino_acid (dict): Mapping from codons to amino acids.
            amino_acid_to_codons (dict): Mapping from amino acids to their codons.
        """
        self.hairpin_threshold = hairpin_threshold
        self.c2aa = c2aa or {}
        self.aa2c = aa2c or {}

    def is_hairpin(self, sequence):
        """
        Checks if a given sequence can form a hairpin structure.

        Parameters:
            sequence (str): The DNA sequence to evaluate.

        Returns:
            bool: True if the sequence forms a hairpin, False otherwise.
        """
        seq_length = len(sequence)

        for size in range(self.hairpin_threshold, seq_length // 2 + 1):
            left = sequence[:size]
            right = sequence[-size:]
            if self.is_complementary(left, right):
                return True

        return False

    def is_complementary(self, seq1, seq2):
        """
        Checks if two DNA sequences are complementary.

        Parameters:
            seq1 (str): The first DNA sequence.
            seq2 (str): The second DNA sequence.

        Returns:
            bool: True if seq1 and seq2 are complementary, False otherwise.
        """
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return all(complement[base1] == base2 for base1, base2 in zip(seq1, reversed(seq2)))

    def identify_hairpins(self, codons):
        """
        Identifies positions of hairpins in a sequence of codons.

        Parameters:
            codons (list[str]): A list of codons (e.g., ["ATG", "GAA", "TTC"]).

        Returns:
            list[tuple[int, int]]: A list of tuples representing start and end indices of hairpin-forming regions.
        """
        dna_sequence = "".join(codons)
        seq_length = len(dna_sequence)

        hairpins = []
        for size in range(self.hairpin_threshold, seq_length // 2 + 1):
            for start in range(seq_length - 2 * size + 1):
                left = dna_sequence[start : start + size]
                right = dna_sequence[start + size : start + 2 * size]
                if self.is_complementary(left, right):
                    hairpins.append((start, start + 2 * size))

        return hairpins

    def replace_hairpins(self, codons):
        """
        Removes hairpins by replacing problematic codons while preserving protein identity.

        Parameters:
            codons (list[str]): A list of codons (e.g., ["ATG", "GAA", "TTC"]).

        Returns:
            list[str]: A modified list of codons with hairpins removed.
        """
        stop_codons = {"TAA", "TAG", "TGA"}  # Define stop codons
        max_attempts = 100  # Limit to prevent infinite loops

        for attempt in range(max_attempts):
            hairpins = self.identify_hairpins(codons)
            if not hairpins:
                print(f"All hairpins resolved after {attempt} attempts.")  # DEBUG
                break  # No hairpins detected

            print(f"Hairpins detected: {hairpins}")  # DEBUG

            for start, end in hairpins:
                # Identify the codons that form the hairpin
                problem_region = codons[start // 3 : end // 3]
                print(f"Problematic region: {problem_region}")  # DEBUG

                # Replace problematic codons
                for idx, codon in enumerate(problem_region):
                    aa = self.c2aa.get(codon)
                    if not aa:
                        continue  # Skip if codon is unmapped

                    # Get valid replacements excluding the original codon and stop codons
                    valid_codons = [
                        alt_codon
                        for alt_codon in self.aa2c.get(aa, [])
                        if alt_codon != codon and alt_codon not in stop_codons
                    ]

                    if valid_codons:
                        replacement = random.choice(valid_codons)
                        codons[start // 3 + idx] = replacement
                        print(
                            f"Replacing {codon} with {replacement} "
                            f"at position {start // 3 + idx} (AA: {aa})."
                        )  # DEBUG

                # Re-run the hairpin check on the modified region
                modified_region = codons[start // 3 : end // 3]
                if self.identify_hairpins(modified_region):
                    print("Warning: Replacement did not resolve the hairpin.")  # DEBUG

        else:
            print(f"Failed to resolve all hairpins after {max_attempts} attempts.")  # DEBUG

        return codons


if __name__ == "__main__":
    # Example usage
    codon_to_aa = {"ATG": "M", "TAC": "Y", "GAT": "D", "GAA": "E", "TTC": "F"}
    aa_to_codons = {
        "M": ["ATG"],
        "Y": ["TAC", "TAT"],
        "D": ["GAT", "GAC"],
        "E": ["GAA", "GAG"],
        "F": ["TTC", "TTT"],
    }

    hpost_checker = HairpinPostChecker()
    hpost_checker.initiate(hairpin_threshold=6, c2aa=codon_to_aa, aa2c=aa_to_codons)

    codons = ["ATG", "TAC", "GAT", "GAA", "TTC", "GAT", "TAC"]
    print("Original Codons:", codons)

    replaced_codons = hpost_checker.replace_hairpins(codons)
    print("Replaced Codons:", replaced_codons)
