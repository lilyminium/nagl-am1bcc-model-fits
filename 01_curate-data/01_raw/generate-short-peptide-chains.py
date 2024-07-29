"""
This script generates peptide chains of length 1 to 4 amino acids.
"""

import itertools

def main():
    natural_codes = [
        # # Positively charged
        "R",
        "H",
        "K",
        # Negatively charged
        "D",
        "E",
        # Polar uncharged
        "S",
        "T",
        "N",
        "Q",
        # # Special cases
        "G",
        "P",
        "C",
        # # Hydrophobic
        "A",
        "V",
        "I",
        "L",
        "M",
        # # Hydrophobic
        "F",
        "Y",
        "W",
    ]

    

    for length in range(1, 5):
        all_codes = []
        for code in itertools.product(natural_codes, repeat=length):
            joined = "".join(code)
            all_codes.append(joined)
        
        output_file = f"input/peptide-chains-{length}.fasta"
        with open(output_file, "w") as f:
            f.write("\n".join(all_codes))
        
        print(f"Saved {len(all_codes)} peptide chains to {output_file}")

    
if __name__ == "__main__":
    main()
