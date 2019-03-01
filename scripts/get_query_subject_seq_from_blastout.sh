awk -F "\t" 'BEGIN{counter=1}{print ">seq"counter " queryid="$1; print $21; counter+=1}' effectors1-10_outputb.tab > query_seq.fasta
awk -F "\t" 'BEGIN{counter=1}{print ">seq"counter " subjectid="$1; print $22; counter+=1}' effectors1-10_outputb.tab > subject_seq.fasta
