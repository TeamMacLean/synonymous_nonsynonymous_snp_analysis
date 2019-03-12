#!/usr/bin/env python3

# input should be a blast output file, where blast command should have the ouput format set to -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles'
import os, sys

offsets=7

blastfile=open(sys.argv[1])


translate_table={
					"GCA":"A",
					"GCC":"A",
					"GCG":"A",
					"GCT":"A",
					"TGC":"C",
					"TGT":"C",
					"GAC":"D",
					"GAT":"D",
					"GAA":"E",
					"GAG":"E",
					"TTC":"F",
					"TTT":"F",
					"GGA":"G",
					"GGC":"G",
					"GGG":"G",
					"GGT":"G",
					"CAC":"H",
					"CAT":"H",
					"ATA":"I",
					"ATC":"I",
					"ATT":"I",
					"AAA":"K",
					"AAG":"K",
					"CTA":"L",
					"CTC":"L",
					"CTT":"L",
					"CTG":"L",
					"TTG":"L",
					"TTA":"L",
					"ATG":"M",
					"AAC":"N",
					"AAT":"N",
					"CCA":"P",
					"CCC":"P",
					"CCG":"P",
					"CCT":"P",
					"CAA":"Q",
					"CAG":"Q",
					"AGA":"R",
					"AGG":"R",
					"CGA":"R",
					"CGC":"R",
					"CGG":"R",
					"CGT":"R",
					"AGC":"S",
					"AGT":"S",
					"TCA":"S",
					"TCC":"S",
					"TCG":"S",
					"TCT":"S",
					"TAG":"*",
					"TGA":"*",
					"TAA":"*",
					"ACA":"T",
					"ACC":"T",
					"ACG":"T",
					"ACT":"T",
					"GTA":"V",
					"GTC":"V",
					"GTG":"V",
					"GTT":"V",
					"TGG":"W",
					"TAC":"Y",
					"TAT":"Y",
					}
ambiguous_bases={
				"R":['A','G'],
				"Y":['C','T'],
				"S":['G','C'],
				"W":['A','T'],
				"K":['G','T'],
				"M":['A','C'],
				"B":['C','G','T'],
				"D":['A','G','T'],
				"H":['A','C','T'],
				"V":['A','C','G'],
				"N":['A','C','G','T'],
}


def set_syn_or_nonsyn(readingframe_data):

	" checks aminoacid sequence in each reading frame and sets syn or nonsyn"
	result = None
	for query_aa, subject_aa in readingframe_data:
		if ("*" in query_aa) or ("*" in subject_aa):
			continue
		elif query_aa == subject_aa:
			result = "Syn"
		elif query_aa != subject_aa:
			result = "NonSyn"
		else:
			pass

	return result


def get_aminoacid_sequence(queryseq, subjectseq):

	" get 3 reading frame sequence left to right"

	readingframe_aminoacid = []

	for readingframe in range(3):

		queryaminoacidseq=""; subjectaminoacidseq=""
		for position in range(0, len(querysubseq), 3):

			qsubseq=querysubseq[position:position + 3]
			ssubseq=subjectsubseq[position:position + 3]

			if "-" in qsubseq or "-" in ssubseq:
				qaa = 'X'; saa = 'X'
			else:
				# check if there is ambiguous base
				if qsubseq not in translate_table.keys():
					for base in qsubseq:
						if base in ambiguous_bases.keys():
							qsubseq=qsubseq.replace(base, ambiguous_bases[base][0])
				if ssubseq not in translate_table.keys():
					for base in ssubseq:
						if base in ambiguous_bases.keys():
							ssubseq=ssubseq.replace(base, ambiguous_bases[base][0])


				if len(qsubseq) == 3:

					qaa = translate_table[qsubseq]
					saa = translate_table[ssubseq]
				elif len(qsubseq < 3):
					# less than 3 bases not enough for coding amino acid
					pass

			queryaminoacidseq += qaa
			subjectaminoacidseq += saa

			readingframe_aminoacid.append((queryaminoacidseq, subjectaminoacidseq))

	return query_aminoacid, subject_aminoacid


def get_n_bases_around_mismatch(queryseq, subjectseq, offsets, mismatch_positions):
	" get N nucleotide bases on left and right side of a mismatch"
	counter=0
	bases_around_mismatch = []

	for basepos in mismatch_positions:

		if basepos < offsets:
			#print(queryid, subjectid, basepos,, )
			querysubseq =  queryseq[:basepos + offsets*2]; subjectsubseq = subjectseq[:basepos+offsets*2]
		elif basepos + offsets > len(queryseq):
			#print(queryid, subjectid, basepos, )
			querysubseq = queryseq[basepos - offsets*2:]; subjectsubseq = subjectseq[basepos-offsets*2:]
		else:
			#print(queryid, subjectid, basepos, )
			querysubseq = queryseq[basepos - offsets:basepos + offsets]; subjectsubseq = subjectseq[basepos-offsets:basepos+offsets]

		bases_around_mismatch.append((query_bases_around_mismatch, subject_bases_around_mismatch))

		#print_output(queryid, subjectid, basepos, querysubseq, subjectsubseq)
	return bases_around_mismatch


def get_mismatch_positions(queryseq, subjectseq, total_mismatches):

	" get zero (0) based positions of mismatches "

	counter=0; mismatch_positions = []
	for basepos in range(len(queryseq)):

		if (queryseq[basepos] == "-" or subjectseq[basepos] == "-") or queryseq[basepos] == subjectseq[basepos]:
			continue
		else:
			mismatch_positions.append(basepos)

		if counter == total_mismatches:
			break

	return mismatch_positions





if __name__ == "__main__":
	for line in blastfile:
		line=line.rstrip()
		linearray=line.split()
		queryid=linearray[0]
		subjectid=linearray[1]
		mismatch=int(linearray[4])
		queryseq=linearray[20]
		subjectseq=linearray[21]

		#print(queryid, subjectid, mismatch, queryseq, subjectseq)

		if mismatch == 0:
			continue
		else:
			mismatch_positions = get_mismatch_positions(queryseq, subjectseq, mismatch)
			bases_around_mismatch = get_n_bases_around_mismatch(queryseq, subjectseq, offsets, mismatch_positions)
			syn_nonsyn = "None"; resultline = queryseq + "\t" + subjectseq + "\t"
			for querybases, subjectbases in bases_around_mismatch:
				query_aa, subject_aa = get_aminoacid_sequence(querybases, subjectbases)
				syn_nonsyn = set_syn_or_nonsyn(query_aa, subject_aa)
				resultline += query_aa + "\t" + subject_aa + "\t"

			print(resultline + "\t" + syn_nonsyn)
