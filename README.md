

To compile:
g++ -o outputfile -O3 Source.cpp


To use:

AssemblyStatistics:
./AssemblyStatistics contigs YourAssembly.fasta
or
./AssemblyStatistics contigs [YourAssembly.fasta] GS [Estimated_Genome_Size]

SelectLongestReads:
longest 0: select the first reads that sum to total_length bases. 
longest 1: select the longest reads that sum to total_length bases. 
./SelectLongestReads sum [total_length] longest [0 or 1] o [outfile] f [fa/fq_file] f [fa/fq_file] " 
