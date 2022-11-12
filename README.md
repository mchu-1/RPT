# RPT
Python script to create RNA sensors<sup>1-3</sup> for the isolation of PASTE<sup>4</sup> recombinants by reverse promoter trapping (RPT).

## Description
Given an untranscribed target site in the genome, twinPE<sup>5</sup> guide RNAs are generated to insert an attB sequence.

An orphan RNAPIII promoter is then delivered with the cargo into the genome at the attB site and drives the expression of a short RNA transcript from the untranscribed sequence once integrated, or 'trapped', within the genome. Transcripts are terminated by T homopolymers found within the native genomic sequence.

Using RNA sensors, transcripts generated by this process can be targeted to isolate successful recombinant clones.

## Input
- `-i --input` FASTA file containing untranscribed genomic target sequence(s)
- `-o --output` Name for CSV output file
- `-a --attb` AttB sequence for PASTE
- `-n --number` Maximum number of results desired for each genomic target sequence in input (default: 5)

**NOTE:** The sequence of each genomic target in the input file should be sufficiently long for suitable Cas9 protospacers and internal T homopolymers to be found (ideally 500 bp or more).

## Output
CSV file containing the following:
- **target**: names of genomic target
- **region**: sequence of region selected for twinPE within the target
- **transcript**: transcript generated by RPT following PASTE
- **sensor**: RNA sensor targeting the transcript
- **sense_protospacer, antisense_protospacer**: protospacers used to insert the attB sequence via twinPE
- **sense_guide, antisense_guide**: guide RNAs for the attB insertion

Protospacers and guides for twinPE are labelled relative to the strand identity of the target sequence.

For each genomic target, up to *n* top regions are returned with transcripts, sensors and guides for each region.

## Usage
`RPT.py -i <input_FASTA> -o <output_CSV> -a <attB_sequence> -n <maximum_number_of_results_per_target>`

## References
1. Jiang, K., Koob, J., Chen, X.D., Krajeski, R.N., Zhang, Y., Volf, V., Zhou, W., Sgrizzi, S.R., Villiger, L., Gootenberg, J.S., Chen, F., & Abudayyeh, O.O. 2022. Programmable eukaryotic protein synthesis with RNA sensors by harnessing ADAR. Nature Biotechnology.
2. Qian, Y., Li, J., Zhao, S., Matthews, E.A., Adoff, M., Zhong, W., An, X., Yeo, M., Park, C., Yang, X., Wang, B.S., Southwell, D.G. & Huang, J. 2022. Programmable RNA sensing for cell monitoring and manipulation. Nature.
3. Kaseniit, K.E., Katz, N., Kolber, N.S., Call, C.C., Wengier, D.L., Cody, W.B., Sattely, E.S. & Gao, X.J., 2022. Modular and programmable RNA sensing using ADAR editing in living cells. Nature Biotechnology.
4. Ioannidi, E.I., Yarnall, M.T., Schmitt-Ulms, C., Krajeski, R.N., Lim, J., Villiger, L., Zhou, W., Jiang, K., Roberts, N., Zhang, L., Vakulskas, C.A., Walker, J.A., Kadina, A.P., Zepeda, A.E., Holden, K., Gootenberg, J.S., & Abudayyeh, O.O. 2021. Drag-and-drop genome insertion without DNA cleavage with CRISPR-directed integrases. bioRxiv.
5. Anzalone, A.V., Gao, X.D., Podracky, C.J., Nelson, A.T., Koblan, L.W., Raguram, A., Levy, J.M., Mercer, J.A. & Liu, D.R., 2022. Programmable deletion, replacement, integration and inversion of large DNA sequences with twin prime editing. Nature Biotechnology.


