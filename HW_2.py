import sys, math, operator
from collections import Counter
###########################################################################################
# Monte Anderson 
# CSCI 4314
# Homework 2 
# February 14, 2017
#
# ASSUMPTIONS MADE:
#	Assuming we need to account for lowercase and uppercase, we convert everything to
#	uppercase using .upper() function.
#
#
# RUN TIME:
#	Our first loop, "for line in inputFile" will take N time, since it loops through all
#	items in the file. From this, we call our function createDictionary, which loops through
#	and creates a new kmer list from those that we've seen. This will take N time as well,
#	but also k time, as we need to create new K-mers. Therefore, this process will take 
#	O(n^2 * k) time. After this, we need to create the dictionaries for both global and unique
#	occurences. Since these are not nested loops, they will only take N time each. However, 
#	since we also have to sort these later, to get the top results, it will take another N time.
#	So, for each dictionary, we have O(n^2) time.
#	Therefore, the predicted run time would be O(n^2 * k) + O(n^2) + O(n^2)
#
# MEMORY USED:
#	For total memory used, we create the first sequence of strings, which takes up (at worst) 8 bytes.
#	For each k-mer, there can be a total of 4^k possible k-mers. Since we create 5 dictionaries and strings,
#	being occurrenceDict, uniqueDict, occurrenceDictSorted, uniqueDictSorted, and the stringSequence, this will
#	take up about 5 * 8 * 4^k space in a worst case scenerio.
###########################################################################################

def main(argv):

	if (sys.argv[1] == '-f' and sys.argv[3] == '-l' and int(sys.argv[4]) <= 8 and int(sys.argv[4]) >= 3):

		#Here, we initialize all the variables from our given command line arguments.
		lengthKmer = int(sys.argv[4])
		filename = sys.argv[2]
		inputFile = open(filename)

		stringSequence = ""
		occurrenceDict = {}
		uniqueDict = {}
			
		#For every line, we either add to the string sequence, or create/add to a dictionary with
		# all the other kmers.

		for line in inputFile:
			if(">" in line):			
				createDictionary(stringSequence, occurrenceDict, uniqueDict, lengthKmer)	
				stringSequence = ""
			else:
				stringSequence += line.upper().strip()
					
		createDictionary(stringSequence, occurrenceDict, uniqueDict, lengthKmer)

		occurrenceDictSorted = sorted(occurrenceDict.items(), key=operator.itemgetter(1), reverse = True)
		uniqueDictSorted = sorted(uniqueDict.items(), key=operator.itemgetter(1), reverse = True)

		
		print("\n\nFile: " + filename + "\nk-mer Length:  " + str(lengthKmer) + "\n")
		printOutOcc(occurrenceDictSorted)
		printOutSeq(uniqueDictSorted)	
		
		
	else:
		print("\nIncorrect format. Correct input format:\n\n<HW_2>.py -f <Filename>.fasta -l <Length of kmer(3-8)>\n\n")
    


#This function will look at the entire sequence given, and create a dictionary for the unique Seq Count and Occurances.
def createDictionary(sequence, occurrenceDict, uniqueDict, kmer):
	#Creates the empty list for us to use later. This will store all the kmers we found.
	listKmers = []

	#This will create a full dictionary of kmers that we will iterate through below.
	for i in range (len(sequence) - kmer):
		listKmers.append(sequence[i:kmer])
		kmer += 1

	#Looking through the entire kmer list, we see if the current kmer is in the dictionary, and if it is, we increment by 1.
	#If not, we make the new count for it and set it at 1.
	for item in listKmers:
		if(item not in occurrenceDict):
			occurrenceDict[item] = 1
		else:
			occurrenceDict[item] += 1

	#Creating a set will eliminate all the other duplicates, so we can count unique entries.
	setKmers = set(listKmers)
	for item in setKmers:
		if(item not in uniqueDict):
			uniqueDict[item] = 1
		else:
			uniqueDict[item] += 1


def printOutOcc(occurrenceDictSorted):
	print("\nk-mer\tOccurrence")
	for i in range(len(occurrenceDictSorted)):
		if(i < 5):
			print(occurrenceDictSorted[i][0] + "\t" + str(occurrenceDictSorted[i][1]))

		elif(occurrenceDictSorted[i][1] == occurrenceDictSorted[i-1][1]):
			print(occurrenceDictSorted[i][0] + "\t" + str(occurrenceDictSorted[i][1]))
		else:
			break

def printOutSeq(uniqueDictSorted):

	print ("\nk-mer\tSeq Count")
	for i in range(len(uniqueDictSorted)):
		if(i < 5):
			print(uniqueDictSorted[i][0] + "\t" + str(uniqueDictSorted[i][1]))

		elif(uniqueDictSorted[i][1] == uniqueDictSorted[i-1][1]):
			print(uniqueDictSorted[i][0] + "\t" + str(uniqueDictSorted[i][1]))
		else:
			break
	
if __name__ == "__main__":
    main(sys.argv[0:])