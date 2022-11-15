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


def is_in_frame(position_1, position_2: int) -> bool:
    """
    Determine whether a pair of ADAR editing sites are in frame with each other.
    """
    frameshift = (position_2 - position_1) % 3 # frameshift based on distance between edit sites
    if frameshift == 0:
        return True
    else:
        return False


def get_optimal_editing_sites(sequence: str) -> tuple[int, int]:
    """
    Get positions of the optimal editing sites for a sequence.
    """
    position_1 = len(sequence) // 4 # first optimal editing site at approximately the one-quarter position
    position_2 = 3 * position_1 # second optimal editing site at approximately the three-quarter position

    optimal_editing_sites = position_1, position_2
    # distance between optimal editing sites is approximately half the length of the sequence
    # centre position between optimal editing sites is approximately the centre of the sequence

    return optimal_editing_sites


def choose_editing_sites(sequence: str) -> tuple[int, int]:
    """
    Choose ADAR editing sites in a target sequence based on edit distances for all pairs of possible editing sites.
    """
    editing_sites = [site for site in find_editing_sites(sequence) if site < len(sequence) - 24]
    # list of all possible edit sites
    # exclude editing sites within 24 nt of the transcript end to facilitate ADAR binding to dsRNA substrates

    pairwise_sites = [sites for sites in combinations(editing_sites, 2) if is_in_frame(*sites)]
    # list all pairs of editing sites which are in frame with each other

    optimal_sites = get_optimal_editing_sites(sequence)

    ranked_pairs = sorted(pairwise_sites, key = lambda sites: (sites[0] - optimal_sites[0])**2 + (sites[1] - optimal_sites[1])**2)
    # rank pairs of editing sites based on how closely they approximate the optimal editing sites
    # use squared Euclidean distance to determine proximity to the optimum

    best_edit_sites = ranked_pairs[0] # best pair of editing sites

    return best_edit_sites


def generate_target_frame(sequence: str) -> str:
    """
    Choose frame of sequence to be targeted by RNA sensor based on the positions of the best pair of editing sites.
    Mark the editing sites and return the sequence in the chosen frame.
    """
    editing_sites = choose_editing_sites(sequence) # optimal editing sites

    marked_sequence = sequence
    for site in editing_sites:
        marked_sequence = marked_sequence[:site] + "NNN" + marked_sequence[site+3:]
        # mark editing sites with degenerate bases

    left_flank = editing_sites[0] # length of sequence flanking left editing site
    right_flank = len(sequence) - editing_sites[1] # length of sequence flanking right editing site

    left_flank -= left_flank % 3
    right_flank -= right_flank % 3
    # take highest flank lengths divisible by 3

    missing_length = 350 - (editing_sites[1] - editing_sites[0]) # missing length from ideal sensor (350 bp)
    missing_length -= missing_length % 3 # take highest missing length divisible by 3

    if missing_length % 2 == 0:
        left_extension, right_extension = missing_length//2, missing_length//2
    else:
        left_extension, right_extension = (missing_length-3)//2, (missing_length+3)//2
    # length of extensions required to fill missing lengths on left and right of editing sites

    left_index = editing_sites[0]
    right_index = editing_sites[1]
    # begin indexing the ends of the target frame from the editing sites

    if left_extension <= left_flank and right_extension <= right_flank:
        left_index -= left_extension
        right_index += right_extension
    elif left_extension <= left_flank and right_extension > right_flank:
        deficit = right_extension - right_flank
        surplus = left_flank - left_extension
        right_index += right_flank
        if surplus >= deficit:
            left_index -= left_extension + deficit
        else:
            left_index -= left_extension + surplus
    elif left_extension > left_flank and right_extension <= right_flank:
        deficit = left_extension - left_flank
        left_index -= left_flank
        right_index += right_extension + deficit
    else:
        left_index -= left_flank
        right_index += right_flank
    # extend the left and right ends of the target frame as symmetrically as possible
    # use only the available number of bases on the flanks

    sequence_in_frame = marked_sequence[left_index: right_index]
    # completed target frame with marked editing sites

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

    












    










