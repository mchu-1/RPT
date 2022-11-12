# RPT.py

"""Generate sensible transcripts by reverse promoter trapping (RPT) within a target genomic context for the isolation of PASTE recombinants."""

import argparse
from DNA import rev_comp
import RADAR
import pandas as pd


def parse_args():
    """
    Parse input FASTA and output CSV filenames.
    Parse attB sequence for PASTE.
    """
    parser = argparse.ArgumentParser(description = "Generate sensible transcripts by reverse promoter trapping (RPT) within a target genomic context for the isolation of PASTE recombinants.")

    parser.add_argument("-i", "--input", type = str, metavar = "", required = True, help = "Input FASTA file name.")
    parser.add_argument("-o", "--output", type = str, metavar = "", required = True, help = "Output CSV file name.")
    parser.add_argument("-a", "--attb", type = str, metavar = "", required = True, help = "AttB sequence.")
    parser.add_argument("-n", "--number", type = int, metavar = "", required = False, default = 5, help = "Number of transcripts per target (default: 5).")

    args = parser.parse_args()

    return args


def find_terminator(sequence: str, start: int = 0) -> int:
    """
    Find the first RNAPIII termination signal in a sequence.
    """
    terminator = sequence.find("T"*7, start)  # 7 T termination signal

    return terminator


def find_antisense_protospacer(sequence: str, terminator: int) -> int:
    """
    Find the first Cas9 protospacer in the antisense strand of a sequence upstream of an RNAPIII terminator.
    """
    min_length = 200  # minimum length of transcript region for effective RNA sensor targeting
    max_length = 350  # maximum length of transcript region for effective RNA sensor targeting
    pam_proximal = 6  # length of PAM proximal region (excluded from the transcript)
    min_distance_to_pam = min_length + pam_proximal  # minimum distance from terminator to the end of the PAM
    max_distance_to_pam = max_length + pam_proximal  # maximum distance from terminator to the end of the PAM

    start_index = terminator - max_distance_to_pam  # start search here
    end_index = terminator - min_distance_to_pam  # end search here

    if end_index < 0:
        return -1
    elif start_index < 0:
        start_index = 0

    antisense_pam = sequence.rfind("CC", start_index, end_index)  # index of PAM

    if antisense_pam == -1:
        return -1

    antisense_protospacer = antisense_pam + 22  # index of first base in antisense protospacer

    return antisense_protospacer


def find_sense_protospacer(sequence: str, antisense_protospacer: int) -> int:
    """
    Find the first Cas9 protospacer in the sense strand of a sequence 90-150 bp upstream of an antisense protospacer.
    """
    min_cut_distance = 90  # minimum distance between Cas9 cleavage sites for efficient insertion of attB site by twinPE
    max_cut_distance = 150  # maximum distance between Cas9 cleavage sites for efficient insertion of attB site by twinPE
    pam_distal = 17  # length of PAM distal region
    distance_to_GG = 4  # distance from cut site to G dinucleotides in PAM

    min_distance = pam_distal + min_cut_distance - distance_to_GG  # minimum distance to G dinucleotides in distal protospacer
    max_distance = pam_distal + max_cut_distance - distance_to_GG  # maximum distance to G dinucleotides in distal protospacer

    start_index = antisense_protospacer - max_distance  # start search here
    end_index = antisense_protospacer - min_distance  # end search here

    if end_index < 0:
        return -1
    elif start_index < 0:
        start_index = 0

    sense_pam = sequence.rfind("GG", start_index, end_index)  # index of PAM

    if sense_pam == -1:
        return -1

    sense_protospacer = sense_pam - 21  # index of first base in sense protospacer

    if sense_protospacer < 0:
        return -1

    return sense_protospacer


def find_next_antisense_protospacer(sequence: str, antisense_protospacer: int) -> int:
    """
    Find the next antisense protospacer after a given one.
    """
    antisense_pam = antisense_protospacer - 22
    next_antisense_pam = sequence.rfind("CC", 0, antisense_pam)  # index next PAM

    if next_antisense_pam == -1:
        return -1

    next_antisense_protospacer = next_antisense_pam + 22 # index of first base in next protospacer

    return next_antisense_protospacer


def find_end_of_terminator(sequence: str, terminator: int = 0) -> int:
    """
    Find the end of an RNAPIII termination sequence.
    """
    end = 0
    for i in range(terminator, len(sequence)):
        if sequence[i] == "T":
            continue
        else:
            end = terminator + i
            break

    return end


def can_prematurely_terminate(sequence: str) -> bool:
    """
    Determine whether a transcript can prematurely terminate based on the length of internal T homopolymers.
    """
    premature_terminator = "T" * 5  # minimal efficient transcript release signal in yeast is a T homopolymer of length 5 (http://dx.doi.org/10.1016/j.molcel.2015.04.002)

    if premature_terminator in sequence:
        return True
    else:
        return False


def score_stops(transcript: str) -> float:
    """
    Score a transcript based on the number of translational stops in any frame of the antisense sequence.
    """
    stops = [rev_comp(s) for s in ["TAG", "TAA", "TGA"]]  # antisense sequence of stops

    stop_count = 0
    for i in range(0, len(transcript)-2):
        if transcript[i:i+3] in stops:  # count number of stops in any frame
            stop_count += 1
        else:
            continue

    stop_score = stop_count / len(transcript)
    # number of stops as a proportion of total codon space (number of codons across all frames of transcript)
    # a lower score predicts fewer stops that have to be removed from an RNA sensor targeting the transcript

    return stop_score


def find_target_transcripts(sequence: str, k: int) -> list:
    """
    Find target transcripts within a sequence based on stop scores.
    Return the coordinates of up to k transcripts with the lowest stop scores and their protospacers.
    """
    top_coords = [None]*k
    top_scores = [1]*k
    terminator = 0

    while True:

        start_index = find_end_of_terminator(sequence, terminator)  # start search downstream of terminator
        terminator = find_terminator(sequence, start_index)  # find the first terminator

        if terminator == -1:
            break

        antisense_protospacer = find_antisense_protospacer(sequence,
                                                          terminator)  # find proximal protospacer

        while antisense_protospacer >= 0:
            sense_protospacer = find_sense_protospacer(sequence,
                                                        antisense_protospacer)  # find distal protospacer

            if sense_protospacer == -1:
                break

            transcript_start = antisense_protospacer - 16  # start index of new transcript from Cas9 cleavage site of proximal protospacer
            transcript_end = terminator

            new_transcript = sequence[transcript_start:transcript_end]  # new transcript

            if can_prematurely_terminate(new_transcript):
                break

            new_score = score_stops(new_transcript)  # score new transcript by number of stops

            for i in range(k):
                if top_scores[i] < new_score:  # compare stop score to top scores (lower is better)
                    continue
                else:
                    top_scores[i] = new_score # update top scores
                    top_coords[i] = (sense_protospacer, antisense_protospacer,
                                     terminator)  # coordinates of protospacers and transcript termination with top scores
                    break

            antisense_protospacer = find_next_antisense_protospacer(sequence,
                                                                   antisense_protospacer)  # find next proximal protospacer

    top_coords = [coord for coord in top_coords if coord]  # take coordinates of top k or fewer transcripts and their protospacers

    return top_coords


def generate_twinpe_guide(protospacer: str, template: str, pbs_length: int = 13, overlap_length: int = 30) -> str:
    """
    Generate guide RNA to insert attB template at a target protospacer using twinPE.
    """
    motif = "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA"  # tevopreq1 motif
    scaffold = "GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTCGAAAGAGTGGCACCGAGTCGGTGCT"  # F+E cr772 tracrRNA

    if not protospacer[0] == "G":  # add 5' G to spacer if necessary
        spacer = "G" + protospacer
    else:
        spacer = protospacer

    pbs = spacer[:-3][-pbs_length:]
    rtt_length = (len(template)+overlap_length)//2  # length of RTT for twinPE insertion based on overlap between guides
    rtt = template[:rtt_length]

    guide = spacer + scaffold + rev_comp(pbs + rtt) + motif

    return guide

def generate_sensors_and_guides(input_filename, output_filename, attb_sequence: str, number_of_transcripts: int) -> None:
    """
    Generate RNA sensors for selection of PASTE recombinants based on top RNAPIII transcripts.
    """
    fasta_names = []
    fasta_sequences = []

    print(f"Opening file {input_filename} ...")
    with open(input_filename) as input_fasta:
        for line in input_fasta:
            if line[0] == ">":
                fasta_names.append(line[1:].rstrip())
                fasta_sequences.append("")
            else:
                fasta_sequences[-1] += line.rstrip()

    targets = []
    regions = []
    transcripts = []
    sensors = []
    sense_protospacers = []
    antisense_protospacers = []
    sense_guides = []
    antisense_guides = []

    for name, sequence in zip(fasta_names, fasta_sequences):

        print(f"Finding up to {number_of_transcripts} transcripts for {name} ...")
        target_coords = find_target_transcripts(sequence, number_of_transcripts)
        total_transcripts = len(target_coords)

        print(f"{total_transcripts} transcripts found for {name}.")

        print(f"Generating sensors and guides for {name} ...")
        for coord in target_coords:
            targets.append(name)

            sense_index, antisense_index, termination = coord
            # index of sense protospacer, antisense protospacer and termination of transcript

            region = sequence[sense_index: termination+7]
            regions.append(region)

            target_transcript = sequence[antisense_index-16: termination]
            transcripts.append(target_transcript)

            target_sensor = RADAR.generate_sensor(target_transcript) # generate RNA sensor against transcript
            sensors.append(target_sensor)

            sense_protospacer = sequence[sense_index: sense_index+20]
            sense_protospacers.append(sense_protospacer)

            antisense_protospacer = rev_comp(sequence[antisense_index-19: antisense_index+1])
            antisense_protospacers.append(antisense_protospacer)

            sense_guide = generate_twinpe_guide(sense_protospacer, attb_sequence)
            sense_guides.append(sense_guide) # generate sense guide

            antisense_guide = generate_twinpe_guide(antisense_protospacer, rev_comp(attb_sequence))
            antisense_guides.append(antisense_guide) # generate antisense guide

        print(f"Sensors and guides generated for {name}.")

    print("Sensors and guides generated for all transcripts.")

    df = pd.DataFrame({"target": targets,
                       "region": regions,
                       "transcript": transcripts,
                       "sensor": sensors,
                       "sense_protospacer": sense_protospacers,
                       "antisense_protospacer": antisense_protospacers,
                       "sense_guide": sense_guides,
                       "antisense_guide": antisense_guides})

    print(f"Writing to {output_filename}.csv ...")
    df.to_csv(output_filename + ".csv", sep = ",", encoding = "utf-8", index = False)

    print("Completed successfully!")

if __name__ == "__main__":

    args = parse_args()
    input_filename = args.input
    output_filename = args.output
    attb_sequence = args.attb.upper()
    number = args.number

    print("-------------------------------------------------------------------")
    print("Reverse Promoter Trapping (RPT) for Isolation of PASTE Recombinants")
    print("-------------------------------------------------------------------")

    if input_filename.split(".")[1] in ["fasta", "fa"]: # check input filename suffix is FASTA
        generate_sensors_and_guides(input_filename, output_filename, attb_sequence, number)
    else:
        print("ERROR: Input file should be FASTA format.")



















































