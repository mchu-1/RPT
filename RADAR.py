# RADAR.py

"""Generate optimal RNA sensor to detect a target transcript using endogenous ADAR."""

from itertools import combinations
from DNA import rev_comp


def is_editing_site(anticodon: str) -> bool:
    """
    Determine if an anticodon is a possible ADAR editing site.
    """
    edit_sites = ["CCA", "GCA", "TCA"] # base-pairing sites for UAG sensor stop codons with roughly equivalent ADAR editing efficiency (https://doi.org/10.1038/s41587-019-0178-z)
    if anticodon in edit_sites:
        return True
    else:
        return False

def find_editing_sites(frame: str) -> list[int]:
    """
    Find starting positions of all possible ADAR editing sites in a reading frame.
    """
    edit_positions = []
    for i in range(0, len(frame)-2):
        if is_editing_site(frame[i:i+3]):
            edit_positions.append(i) # starting position of edit site
        else:
            continue

    return edit_positions


def distance_between_editing_sites(position_1, position_2: int) -> int:
    """
    Calculate distance between ADAR editing sites.
    """
    return position_2 - position_1

def is_in_frame(position_1, position_2: int) -> bool:
    """
    Determine whether a pair of ADAR editing sites are in frame with each other.
    """
    frameshift = distance_between_editing_sites(position_1, position_2) % 3 # frameshift based on distance between edit sites
    if frameshift == 0:
        return True
    else:
        return False


def choose_editing_sites(sequence: str) -> tuple[int]:
    """
    Choose ADAR editing sites in a target sequence based on edit distances for all pairs of possible editing sites.
    """
    editing_sites = [site for site in find_editing_sites(sequence) if site < len(sequence) - 24]
    # list of all possible edit sites
    # exclude editing sites within 24 nt of the transcript end to facilitate ADAR binding to dsRNA substrates

    pairwise_sites = [sites for sites in combinations(editing_sites, 2) if is_in_frame(*sites)]
    # list all pairs of editing sites which are in frame with each other
    optimal_distance = len(sequence) // 2 # optimal distance between editing sites set at half the length of the target sequence
    ranked_pairs = sorted(pairwise_sites, key = lambda sites: abs(distance_between_editing_sites(*sites) - optimal_distance))
    # rank pairs of editing sites based on how closely the distance between them approximates the optimal distance

    best_edit_sites = ranked_pairs[0] # best pair of editing sites

    return best_edit_sites


def generate_target_frame(sequence: str) -> str:
    """
    Choose frame of sequence to be targeted by RNA sensor based on the positions of the optimal pair of editing sites.
    Mark the editing sites and return the sequence in the chosen frame.
    """
    editing_sites = choose_editing_sites(sequence) # optimal editing sites

    marked_sequence = sequence
    for site in editing_sites:
        marked_sequence = marked_sequence[:site] + "NNN" + marked_sequence[site+3:]
        # mark editing sites with degenerate bases

    left_index = editing_sites[0] % 3  # index of the first base in frame with the editing sites
    right_offset = len(sequence[left_index:]) % 3  # number of bases at the end of the sequence which are out of frame with the editing sites
    right_index = len(sequence) - right_offset # index of last base in frame with the editing sites

    sequence_in_frame = marked_sequence[left_index:right_index]  # sequence in frame with the optimal editing sites

    return sequence_in_frame


def remove_stops(frame: str) -> str:
    """
    Remove stops from the antisense sequence of a reading frame.
    """
    stops = ["TAG", "TAA", "TGA"] # sense sequence of stop codons

    new_frame = frame # copy to new frame
    for i in range(0, len(new_frame), 3):
        if new_frame[i:i+3] in stops:
            new_frame = new_frame[:i] + "C" + new_frame[i+1:]
            # substitute third base for 'C'
        else:
            continue

    return new_frame # new frame with stops removed


def remove_bsmbi_sites(sequence: str) -> str:
    """
    Remove BsmBI sites from a sequence for cloning.
    """
    bsmbi = ["CGTCTC", rev_comp("CGTCTC")]  # sense and antisense BsmBI recognition sites

    new_sequence = sequence
    for i in range(0, len(new_sequence)-5):
        if new_sequence[i:i+6] in bsmbi:
            new_sequence = new_sequence[:i+1] + "C" + new_sequence[i+2:] # substitute second base for 'C'
        else:
            continue

    return new_sequence


def generate_sensor(sequence: str) -> str:
    """
    Generate RNA sensor for a target sequence.
    """
    target_frame = generate_target_frame(sequence) # generate target frame in sequence
    sensor_sequence = rev_comp(target_frame) # generate sensor sequence
    sensor_sequence = remove_stops(sensor_sequence) # remove stops
    sensor_sequence = remove_bsmbi_sites(sensor_sequence)  # remove BsmBI recognition sites
    sensor_sequence = sensor_sequence.replace("NNN", "TAG") # install stop substrates for ADAR editing

    return sensor_sequence

    












    










