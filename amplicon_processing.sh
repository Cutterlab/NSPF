#!/bin/bash

#read in directories
bioawk_dir=/PathTobioawk/bioawk
read_dir=/PathToRawReads/00_fastq
out_dir=/PathForProcessedReads/processed_reads
blast_dir=/PathtoBestBlast/BestBlast

#create BLAST database fasta file
#>DEL_amplicon
#AGGATAAGAGAGACCCAGACAGTGAGACACTCTAGCTTCCCTGCGTCTCTCTTGTTCTCT
#ACAGGTCTCTAGGTCTCTGTAACTAGGAGTTTAATGATATTTTGTTGTAAGAATACAATT
#TTTCTTTGTGACTTACCGTCTGGTAGGGTGGCAGATCAGTGTTCAGAAGGAAGTGATGGA
#AATGCAAATGTGCTAAAAATTTTGATTGCTCGCCGACTTGGCAAATGGGAAAATTCAGAT
#ATTCCGAAAACCAGCCTCACTTTCCTACTGCCAGA

#read in BLAST database
blast_db=/PathtoBlastDataBase/

#run each file at a time
for f in /PathtoRawReads/*.fastq
do
	base=`basename $f _001.fastq`
	
        #filter reads by mean quality score (Phred Score > 30) and deletion length (>230bp) or all other lengths (<= 230bp & > 50bp)
	${bioawk_dir}/bioawk -c fastx '{ if(meanqual($qual) > 30 && length($seq) > 230) { print "@"$name; print $seq; print "+"; print $qual; }}' $f > ${out_dir}/${base}_FILTEREDforDEL.fastq
	${bioawk_dir}/bioawk -c fastx '{ if(meanqual($qual) > 30 && length($seq) <= 230 && length($seq) > 50) { print "@"$name; print $seq; print "+"; print $qual; }}' $f > ${out_dir}/${base}_FILTEREDelse.fastq

	#convert FASTQ file to FASTA file format
	${bioawk_dir}/bioawk -c fastx '{print ">"$name; print $seq}' ${out_dir}/${base}_FILTEREDforDEl.fastq > ${out_dir}/${base}_FILTEREDforDEL.fasta
	${bioawk_dir}/bioawk -c fastx '{print ">"$name; print $seq}' ${out_dir}/${base}_FILTEREDelse.fastq > ${out_dir}/${base}_FILTEREDelse.fasta

	#blast to deletion db
	blastn -subject ${blast_db}/DEL_amplicon.fasta -query ${out_dir}/${base}_FILTEREDforDEL.fasta -outfmt 6 -out ${blast_dir}/${base}_DEL.txt
	blastn -subject ${blast_db}/unique_amplicon.fasta -query ${out_dir}/${base}_FILTEREDelse.fasta -outfmt 6 -out ${blast_dir}/${base}_unique.txt
done
