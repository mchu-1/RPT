# DNA.py

"""Essential functions for DNA cloning."""


def rev_comp(sequence: str) -> str:
    """
    Generate reverse complement of a DNA sequence.
    """
    bp_dict = {"A":"T","T":"A","C":"G","G":"C"} # base-pairing rules for DNA
    bp_table = {ord(i):ord(j) for i,j in bp_dict.items()} # converting rules to ASCII bytes
    return sequence.translate(bp_table)[::-1] # generate reverse complement using rules