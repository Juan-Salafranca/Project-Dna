import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'DNAPROT_package')))

from DNAPROT_analysis import *
from DNAPROT_shinedalgarno import *

# test find_shine_dalgarno_function

def test_find_shine_dalgarno():
    sequences = [
        ("ATGAGGAGGTAG", 3),
        ("ATGCCCTAG", None),
        ("AGGAGGATG", 0),
        ("ATGAGGAGG", 3),
        ("GGAGGAGG", 2),
    ]
    
    for sequence, expected in sequences:
        result = find_shine_dalgarno(sequence)
        assert result == expected, f"Test failed: did not understand {sequence}, expected {expected}, got {result}"

# test cut_sequence function

def test_cut_sequence():
    sequences = [
        ("ATGAGGAGGTAG", ["TAG"]),
        ("ATGCCCTAG", []),
        ("AGGAGGATG", ["ATG"]),
        ("ATGAGGAGG", [""]),
        ("GGAGGAGGA", ["A"]),
        ("AGGAGGAGGAGGAGGAGG", ["", "", ""])]
    
    for sequence, expected_sections in sequences:
        result = cut_sequence(sequence)
        assert result == expected_sections, f"Test failed: did not understand {sequence}, expected {expected_sections}, got {result}"

# test translate_to_uppercase function

def test_translate_to_uppercase():
    test_cases = [
        ("atgc", "ATGC"),
        ("ATGC-", "ATGC-"),
        ("aTgC3", "ATGC3"),
        ("", "")]
    
    for sequence, expected in test_cases:
        result = translate_to_uppercase(sequence)
        assert result == expected, f"Test failed: did not understand {sequence}, expected {expected}, got {result}"

def test_filter_dna_sequence():
    # Test cases for filter_dna_sequence function
    test_cases = [
        ("ATGC", "ATGC"),
        ("ATGCK", "ATGC"),
        ("atgc", "ATGC"),
        ("atgcn", "ATGC"),
        ("", "")]
    
    # Iterate through test cases
    for sequence, expected in test_cases:
        result = filter_dna_sequence(sequence)
        assert result == expected, f"Test failed: did not understand {sequence}, expected {expected}, got {result}"

def test_read_dna_sequence():
    # Test cases for read_dna_sequence function
    test_files = [
        ("test_file_1.txt", "AGGAGGATGAGGTCTGAT")]

    # Iterate through test cases
    for file_name, expected_sequence in test_files:
        # Get the path to the test file
        file_path = os.path.join("tests", file_name)

        # Call the function to read DNA sequence from the test file
        result = read_dna_sequence(file_path)

        # Compare the result with the expected sequence
        assert result == expected_sequence, f"Test failed: did not understand {file_name}, expected {expected_sequence}, got {result}"

# test for read_dna_sequence

# test for ReadShineDalgarnoFromTxt

# test for seperate_sections

# test for read_genetic_code

# test for transcribe_dna_to_rna

def test_transcribe_dna_to_rna():
    assert transcribe_dna_to_rna("ATCG") == "AUCG"
    assert transcribe_dna_to_rna("ATTGC") == "AUUGC"
    assert transcribe_dna_to_rna("GCTA") == "GCUA"

# test for find_start_codons_rna function

def test_find_start_codons_rna():
    assert find_start_codons_rna("AUGUUUAUG") == [0, 6]
    assert find_start_codons_rna("UUUAUGUUUGA") == [3]
    assert find_start_codons_rna("UAUGUUUA") == [1]

# test for translate_rna_to_protein

def test_translate_rna_to_protein():
    assert translate_rna_to_protein("UUUGAUGUCUCUAGUCGAUCGAUUGCUUU") == "MSLVDRLL"
    assert translate_rna_to_protein("UUUGAUGUCUCUAGCAGUCGAUCAUGUCCUGAUCCGU") == "MSLAVDHVLIR"

# Test for complementary_sequences function

def test_complementary_sequences():
    assert complementary_sequences("ATCG") == "TAGC"
    assert complementary_sequences("GCTA") == "CGAT"
    assert complementary_sequences("TTAAGC") == "AATTCG"

# Test for flip_rna_sequence function

def test_flip_rna_sequence():
    assert flip_rna_sequence("AUGUUUUAG") == "GAUUUUGUA"
    assert flip_rna_sequence("UUUGUUUGA") == "AGUUUGUUU"
    assert flip_rna_sequence("AUAGUUUA") == "AUUUGAUA"
    assert flip_rna_sequence("") == ""







