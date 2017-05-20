import sys, getopt, math, os.path, re

def main(argv):
    
	concatenated_Seq = ""
	actualStrand = ""
	reverse_kmer = ""
	name = "monte"
	ID = 0
	GC_Count = 0

	#if (sys.argv[1] == '-h'):
		#print("To use this code, enter '-f (filename) -s (sequence name) -k (Kmer name)'")

	
	'''Checks to see if we have all the required -f, -s, -k, arguments.'''
	if ((sys.argv[1] == '-f' or sys.argv[1] == '-F') & (sys.argv[3] == '-s' or sys.argv[3] == '-S') & (sys.argv[5] == '-k' or sys.argv[5] == '-K')):

		'''If this is all true, then we make the filename and sequence and kmer equal to the given inputs.'''
		filename = sys.argv[2]
		sequence = str(">" + sys.argv[4])
		kmer = str(sys.argv[6])
		'''Get the length of the kmer, used for math later on.'''
		kmer_length = len(kmer)

		fileOutput = open(filename)
		boolCheck = False


		'''Here, we create the reverse kmer sequence by getting the complements, then using the [::-1] to reverse the list.'''
		for i in kmer:
			if i == "T":
				reverse_kmer = reverse_kmer + "A"
			if i == "A":
				reverse_kmer = reverse_kmer + "T"
			if i == "C":
				reverse_kmer = reverse_kmer + "G"
			if i == "G":
				reverse_kmer = reverse_kmer + "C"

		'''This reverses the list from above, which gives us the reverse complement to search for.'''
		reverse_kmer = reverse_kmer[::-1]

		#print(kmer, reverse_kmer)
		
		'''For every line, we'll need to check and see if it's a sequence name (starts with ">") or not.
		If it's not, we can concatenate the non-sequence names into a giant list to search through.
		We do this because the list may be seperated by tabs or new lines, but still are part of the same strand.'''
		
		for line in fileOutput:

			if(boolCheck == True):
				if(not line.startswith(">")):
					concatenated_Seq = concatenated_Seq + line.strip()

			'''If the line starts with ">", we designate it as the strand name, and see if the sequence we are looking for is
			contained in this strand name.'''

			if line.startswith(">"):
				strand_name = line
				#print(strand_name)

				if (sequence.strip(">") in strand_name):
					actualStrand = strand_name
					boolCheck = True
				else :
					boolCheck = False
				
		'''To account for some cases, we add .upper, to convert all lowercase letters to uppercase.'''
		
		concatenated_Seq.upper()
		#print len(concatenated_Seq)
				
		for i in range(0,len(concatenated_Seq)):
			if(concatenated_Seq[i] == "C" or concatenated_Seq[i] == "G"):
				GC_Count += 1
		
		if len(concatenated_Seq) == 0:
			print("No sequence to evaluate. Exiting.")
			exit(2)		
		else:
			percentGC = round(100.0 * GC_Count / len(concatenated_Seq),0)	

		'''Writing to the file'''

		outFile = open("output.txt", 'w')
		#outFile.write("%s \t %s \t %s \t %i \t" % (filename, actualStrand.strip(">").replace("\n",""), kmer, len(concatenated_Seq)))
		outFile.write("File:\t" + filename +"\nSeq:\t" + actualStrand.strip(">").replace("\n","") + "\nk-mer:\t" + kmer + "\nReverse k-mer:\t" + reverse_kmer + "\nSeq Length:\t" + str(len(concatenated_Seq)) + "\n")		
		outFile.write("GC Content:\t" + str(percentGC) + "%\n\n")			
		
		ID = 1
		'''For every letter in the length of the new concatenated sequence, we have to loop through and see if our kmer is in it.	'''
		for i in range(0, len(concatenated_Seq) - len(kmer)):
			
			if concatenated_Seq[i:i+len(kmer)] == kmer:
				outFile.write(actualStrand.strip(">").replace("\n","") + "\t" + name + "\tmatch\t" + str(i + 1) + "\t" + str(i + len(kmer)) + "\t100\t" + "+" + "\tID: "+ str(ID) + "\n")
				#print("FOUND" , str(concatenated_Seq[i:i+len(kmer)]), str(kmer))
				ID+=1

			elif concatenated_Seq[i:i+len(kmer)] == reverse_kmer:
				outFile.write(actualStrand.strip(">").replace("\n","") + "\t" + name + "\tmatch\t" + str(i + 1) + "\t" + str(i + len(kmer)) + "\t100\t" + "-" + "\tID: "+ str(ID) + "\n")
				ID+=1
		outFile.close	
		#print(concatenated_Seq[0])

if __name__ == "__main__":
    main(sys.argv[1:])

