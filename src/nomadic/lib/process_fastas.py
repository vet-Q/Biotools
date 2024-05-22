

def load_fasta_as_dict(fasta_path):
    """
    Load a `.fasta` file as a dictionary

    """

    dt = {}
    with open(fasta_path, "r") as fasta:

        # Initialise
        header = None
        seq = ""
        line = fasta.readline().rstrip()
        
        # Iterate
        while line:
            if line.startswith(">"):
                if header is not None:
                    dt[header] = seq
                header = line
                seq = ""
            else:
                seq += line
            line = fasta.readline().rstrip()
    return dt 


def write_fasta_from_dict(input_dt, output_fasta):
    """
    Write a `.fasta` file to `output_fasta` from an input dictionary
    `input_dt`

    """
    with open(output_fasta, "w") as fasta:
        for header, seq in input_dt.items():
            fasta.write(f">{header}\n")
            fasta.write(f"{seq}\n")

