## Introduction

This is a script to check where the sequences, that are mapped using blastn tool, have mismatches and the mismatches give synonymous or non synonymous amino acid change in the protein sequence. 

## Method

The script takes in the blastn output. While blasting, please use the option -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles'. This will give the mapped query and subject sequences in the blastn output and it will be compatible to the script too.

These data is captured from the blastn output -  queryid, subjectid, mismatch number,  query mapped sequence and subject mapped sequences .

If there is no mismatches (0 mismatch), nothing is not with the alignment data. 

If there is a mismatch, the position of the mismatch is obtained and the few nucleotide sequences (set by variable offset in the script) before and after the mismatch is taken from both query and subject aligned sequences. The nucleotide sequences are then translated to protein sequences in 3 open reading frames (left to right) and then they are compared.

If there is a character "-" in the codon while translating, the amino acid X is used.

If there is a stop codon in the protein sequence, they are not considered. So protein sequences (without X and stop codon) from all open reading frames are compared if there is change in amino acid.

If all open reading frame show there is change in amino acid, the mismatch actually gives Non synonymous mutation. If one of the open reading frame show no change in amino acid, the mismatch gives synonymous mutation.

## Requirement

python version >3.6

## Usage

python scripts/syn_nonsyn_from_blast.py data/effectors1-10_outputb.tab

## Output Sample
```
Qid, subjectid, mismatch_position, qframe1,sframe1,qframe2,sframe2,qframe3,sframe3,result
DP_13,PP1.0_contig14247,3,QQXIKK,QKXIKK,NXQ*KK,KXQ*KK,TXNKK,KXNKK,NS
DP_13,PP1.0_contig14247,55,RV*QQ,RVLQQ,EFNKK,EFYKK,SLTK,SFTK,NS
DP_13,PP1.0_contig14247,292,KTRDD,KTQDD,KLETT,KLKTT,NSRR,NSRR,S
```

Note:

S - synonymous amino acid change

NS - non synonymous amino acid change
