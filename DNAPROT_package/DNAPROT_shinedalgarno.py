

def find_shine_dalgarno(sequence, shine_dalgarno="AGGAGG"):
    """Find the Shine-Dalgarno sequence in the given DNA sequence."""
    index = sequence.find(shine_dalgarno)
    if index != -1:
        return index
    else:
        return None



def cut_sequence(sequence, shine_dalgarno="AGGAGG"):
    """Cut the DNA sequence based on the Shine-Dalgarno sequence."""
    sections = []
    start_index = find_shine_dalgarno(sequence, shine_dalgarno)
    if start_index is not None:
        start_index += len(shine_dalgarno)
        while True:
            index = find_shine_dalgarno(sequence[start_index:], shine_dalgarno)
            if index is not None:
                sections.append(sequence[start_index:start_index + index])
                start_index += index + len(shine_dalgarno)
            else:
                sections.append(sequence[start_index:])
                break
    return sections



def translate_to_uppercase(sequence):
    """Translate the DNA sequence to uppercase."""
    return sequence.upper()




def filter_dna_sequence(sequence):
    """Filter out characters that are not 'A', 'T', 'C', or 'G'."""
    return ''.join(filter(lambda x: x in 'ATCG', sequence.upper()))




def read_dna_sequence(filename):
    """Read DNA sequence from the document where the first column is empty."""
    sequence = ""
    read_started = False
    with open(filename, 'r') as file:
        for line in file:
            if not read_started:
                if line.startswith(" "):
                    read_started = True  
            else:
                if line.strip() == "//":  
                    read_started = False  
                else:
                    start_index = line.find(' ')
                    while start_index != -1:  
                        stop_index = line.find('//', start_index)
                        if stop_index == -1:
                            stop_index = len(line) 
                        sequence += line[start_index:stop_index]
                        start_index = line.find(' ', stop_index)
    return sequence


def ReadShineDalgarnoFromTxt(filename: str):
    """Read DNA sequence from the document and process it."""
    if not isinstance(filename, str):
        raise TypeError("Filename must be a string")
    
    
    dna_sequence = read_dna_sequence(filename)
    
    
    dna_sequence = filter_dna_sequence(dna_sequence)
    sections = cut_sequence(dna_sequence)
    
    # Write resulting sections to a file
    with open("output.txt", "w") as file:
        for i, section in enumerate(sections):
            file.write(section + "\n")
            if i < len(sections) - 1:
                file.write("//\n") 
        file.write("//\n")
    
    print("Sections written to output.txt")


if __name__ == "__main__":
    ReadShineDalgarnoFromTxt(filename)
