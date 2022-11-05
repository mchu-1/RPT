# RPT
Python script to generate RNA sensors<sup>1-3</sup> for the isolation of PASTE<sup>4</sup> recombinants by reverse promoter trapping (RPT).

## Description
Given an untranscribed target region in the genome, twinPE<sup>5</sup> guide RNAs are generated to insert a desired attB sequence.

An RNAPIII promoter is then delivered with the cargo into the genome at the attB site and creates a short RNA transcript from the untranscribed sequence.

This can then be targeted using a custom RNA sensor to isolate successful recombinant clones.

## Input
- FASTA file containing untranscribed genomic target sequence(s)
- AttB sequence for PASTE

## Output
CSV file containing the following:
- Names of genomic targets
- Sequences of genomic targets
- Transcripts generated by RPT following PASTE
- RNA sensor targeting the transcripts
- Protospacers used to insert the attB sequence via twinPE
- Guide RNAs for the attB insertion

## Usage
'RPT.py -i <input_FASTA> -o <output_CSV> -a <attB_sequence>'

## References
1. Jiang, K., Koob, J., Chen, X.D., Krajeski, R.N., Zhang, Y., Volf, V., Zhou, W., Sgrizzi, S.R., Villiger, L., Gootenberg, J.S., Chen, F., & Abudayyeh, O.O. 2022. Programmable eukaryotic protein synthesis with RNA sensors by harnessing ADAR. Nature Biotechnology.
2. Qian, Y., Li, J., Zhao, S., Matthews, E.A., Adoff, M., Zhong, W., An, X., Yeo, M., Park, C., Yang, X., Wang, B.S., Southwell, D.G. & Huang, J. 2022. Programmable RNA sensing for cell monitoring and manipulation. Nature.
3. Kaseniit, K.E., Katz, N., Kolber, N.S., Call, C.C., Wengier, D.L., Cody, W.B., Sattely, E.S. & Gao, X.J., 2022. Modular and programmable RNA sensing using ADAR editing in living cells. Nature Biotechnology.
4. Ioannidi, E.I., Yarnall, M.T., Schmitt-Ulms, C., Krajeski, R.N., Lim, J., Villiger, L., Zhou, W., Jiang, K., Roberts, N., Zhang, L., Vakulskas, C.A., Walker, J.A., Kadina, A.P., Zepeda, A.E., Holden, K., Gootenberg, J.S., & Abudayyeh, O.O. 2021. Drag-and-drop genome insertion without DNA cleavage with CRISPR-directed integrases. bioRxiv.
5. Anzalone, A.V., Gao, X.D., Podracky, C.J., Nelson, A.T., Koblan, L.W., Raguram, A., Levy, J.M., Mercer, J.A. & Liu, D.R., 2022. Programmable deletion, replacement, integration and inversion of large DNA sequences with twin prime editing. Nature Biotechnology.


