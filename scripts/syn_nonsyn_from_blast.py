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

def print_output(queryid, subjectid, basepos, querysubseq, subjectsubseq):
	# translate the subsequence to amino acide sequences
	synonymous_result="NS"
	line=queryid + "," + subjectid + "," + str(basepos)
	for readingframe in range(3):
		#print(querysubseq, subjectsubseq)
		queryaminoacidseq=""
		subjectaminoacidseq=""

		for position in range(0, len(querysubseq), 3):

			qsubseq=querysubseq[position:position + 3]
			ssubseq=subjectsubseq[position:position + 3]

			#print('translating ', qsubseq, ssubseq)
			if "-" in qsubseq or "-" in ssubseq:
				qaa = 'X'; saa = 'X'
			else:
				if len(qsubseq) == 3:
					qaa = translate_table[qsubseq]
					saa = translate_table[ssubseq]

			queryaminoacidseq += qaa
			subjectaminoacidseq += saa

			if ("*" in queryaminoacidseq or "*" in subjectaminoacidseq):
				continue
			elif queryaminoacidseq == subjectaminoacidseq:
				synonymous_result="S"
			else:
				synonymous_result="NS"

		line += "," + queryaminoacidseq + "," + subjectaminoacidseq
		querysubseq=querysubseq[1:]
		subjectsubseq=subjectsubseq[1:]


	print(line + "," + synonymous_result)
	#print(queryid, subjectid, basepos, queryaminoacidseq, subjectaminoacidseq)

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
		counter=0
		for basepos in range(len(queryseq)):

			if (queryseq[basepos] == "-" or subjectseq[basepos] == "-") or queryseq[basepos] == subjectseq[basepos]:
				continue
			else:
				counter+=1
				#print('basepos: ', basepos)
				#print("mismatch in position" , basepos, "mismatch counter ", counter, "querybase", queryseq[basepos], "subjectbase", subjectseq[basepos] )
				if basepos < offsets:
					#print(queryid, subjectid, basepos,, )
					querysubseq =  queryseq[:basepos + offsets*2]; subjectsubseq = subjectseq[:basepos+offsets*2]

				elif basepos + offsets > len(queryseq):
					#print(queryid, subjectid, basepos, )
					querysubseq = queryseq[basepos - offsets*2:]; subjectsubseq = subjectseq[basepos-offsets*2:]
				else:
					#print(queryid, subjectid, basepos, )
					querysubseq = queryseq[basepos - offsets:basepos + offsets]; subjectsubseq = subjectseq[basepos-offsets:basepos+offsets]
				print_output(queryid, subjectid, basepos, querysubseq, subjectsubseq)

			if counter == mismatch:
				break
