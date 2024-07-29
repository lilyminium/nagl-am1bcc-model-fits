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
    XS = ["R", "E", "S", "C", "A", "F"]

    combinations = []
    for code in itertools.product(natural_codes, repeat=3):
        joined = "".join(code)
        combinations.append(joined)
    
    PATTERNS = {
        5: "{x}{a}{b}{c}{x}",
        6: "{x}{a}{b}{c}{b}{a}",
        7: "{x}{a}{b}{c}{b}{a}{x}",
        8: "{a}{b}{c}{b}{c}{b}{a}{x}",
        9: "{x}{a}{b}{c}{b}{a}{b}{c}{x}",
        10: "{x}{a}{b}{c}{a}{b}{c}{a}{b}{c}",
        11: "{x}{a}{b}{c}{b}{a}{b}{c}{b}{a}{x}",
        12: "{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}",
        15: "{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}",
        18: "{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}",
        21: "{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}",
        30: "{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}{a}{b}{c}"
    }
    
    
    for length, pattern in PATTERNS.items():
        all_codes = []
        for a, b, c in combinations:
            for x in XS:
                peptide = pattern.format(x=x, a=a, b=b, c=c)
                all_codes.append(peptide)
        all_codes = sorted(set(all_codes))
        output_file = f"input/peptide-chains-{length}.fasta"
        with open(output_file, "w") as f:
            f.write("\n".join(all_codes))
        print(f"Saved {len(all_codes)} peptide chains to {output_file}")

    
if __name__ == "__main__":
    main()
