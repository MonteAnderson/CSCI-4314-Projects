###########################################################################################
#  File            : HW3.py
#
#  Purpose         : Find the best path for aligning sequences of different lengths and same lengths,
#					 as well as creating a score matrix to determine the overall score of the path.
#
#  Developer       : Monte Anderson, Christine Samson, CSCI/MCDB 4314, March 2017
###########################################################################################
#
#   Sample command line arguments to run program:
#   format: python Samson_HW1.py -f <file_name>.fasta
#   example: python Samson_HW1.py -f sample5_mer.fasta
#   -f specifies a <file_name>.fasta input file to be read in and analyzed
#
############################################################################################
#
#The time complexity is O(n^4) + O(n^2), since there are four for loops nested - thus resulting in O(n^4) and the while loop would add another O(n^2). Thus the total would be O(n^4) + O(n^2).
#	The space complexity would be for any (n) matrices (n matrices = n comparisons in the file) of size i * j (where i is the length of the first sequence - can vary from sequence to sequence, 
#	and j is the length of the second sequence - can also vary from sequence to sequence. Then, it would also be the matrix with the length of the list of sequences * length of the list of sequences.
#	Therefore, the total space complexity would be n matrices of size i * j + matrix of length of the list of sequences * length of the list of sequences.
#
############################################################################################

import sys, getopt, math, os.path, re, numpy
	
	#This will create the local matrices for each sequence comparison, which we will use to determine the best path.	
	#The parameters distanceMatrix is our initial empty matrix, sequenceList is our list of all sequences in the file, length_1 is the 
	# length of the first sequence we're looking at, length_2 is the length of the second sequence we're comparing, and c and v are iterators for  sequenceList,
	# which allow us to traverse all sequences in the file.
	
def createMatrix(distanceMatrix, sequenceList, length_1, length_2, c, v):
	#Here, we create our new matrix for each sequence we passed in. This will allow us to determine
	# what is the best path to take next.
	
	matrix = numpy.zeros((length_1+1, length_2+1), dtype=float)
	numpy.set_printoptions(formatter={'float': '{: 0.5f}'.format})
	matrix[1:, 0] = range(1, length_1+1)
	matrix[0, 1:] = range(1, length_2+1)
	
	# Now we need to do this for every sequence.
	Sequence_1 = sequenceList[c]
	Sequence_2 = sequenceList[v]
	
	for i in range(1, len(sequenceList[c])+1):	
		for	j in range(1, len(sequenceList[v])+1):
			if sequenceList[c][i-1] == sequenceList[v][j-1]:
				score = 0
			else:
				score = 1
			
			# This is how we determine the minimum value from our three neighbors.
			matrix[i][j] = min(matrix[i-1][j-1] + score, matrix[i-1][j] + 1, matrix[i][j-1] + 1)

	new_Alignment = findPath(matrix, Sequence_1, Sequence_2)
	
	# The score is the bottom right most cell in our matrix, so we select this to be our score.
	scoreM = matrix[-1][-1]
	createGlobalMatrix(distanceMatrix, scoreM, new_Alignment, c, v)

	#This function will find the path (starting from the bottom right), and create a new sequence that has the best alignment.
	#matrix is the given matrix (from createMatrix), and sequence_1 and sequence_2 are the first and second sequence, respectivley. 
def findPath(matrix, Sequence_1, Sequence_2):
	
	alignment_1 = ""
	alignment_2 = ""
	
	i = len(Sequence_1)
	j = len(Sequence_2)
	
	counter_1 = 0
	counter_2 = 0
	editD = ""	
	
	while i > 0 and j > 0:
		
		# Checking the diagonal first, because we prefer this over the other two.
		if(matrix[i-1][j-1] <= matrix[i][j-1] and matrix[i-1][j-1] <= matrix[i-1][j]):
			alignment_1 += Sequence_1[i - 1]
			alignment_2 += Sequence_2[j - 1] 
			i -= 1
			j -= 1

		# Checking the above cell.
		elif(matrix[i][j-1] <= matrix[i-1][j]):
			alignment_1 += "-" 
			alignment_2 += Sequence_2[j - 1] 
			counter_1 += 1
			j -= 1
			
			
		# Checking the left cell.
		else:
			alignment_1 += Sequence_1[i - 1] 
			alignment_2 += "-" 
			counter_2 += 1
			i -= 1
		editD += "="
						
			
	# Print our alignments.
	print "\n\n", alignment_1[::-1], "\n" ,alignment_2[::-1]
	print editD, len(alignment_2)
	
	# I use this return value to create the part 2 matrix (because we need the length of the new alignment)
	return max(counter_1 + len(Sequence_1), counter_2 + len(Sequence_2))
	
	
	#Creates the final global scoring matrix, and rounds them to 5 decimal places.
	#distanceMatrix is our empty matrix, but we are now filling it with a score (globalScore), the new alignment, and using x and y to iterate through it.
def createGlobalMatrix(distanceMatrix, globalScore, new_Alignment, x, y):
	distanceMatrix[x][y] = globalScore / new_Alignment
	distanceMatrix[x][y] = round(distanceMatrix[x][y], 5)
	distanceMatrix[y][x] = round(distanceMatrix[x][y], 5)
	
	#Calls main function, which runs the entire program.
def main(argv):
	
	#We check to see if the -f is included, and if it is, we make the filename equal to the second argument, after the -f.
	if (sys.argv[1] == '-f'):
		stringSeq = ""
		sequenceList = []
		sequenceCounter = 0
		counter = 1
		length_1 = 0
		length_1 = 0
		
		# Here we are parsing the file to get all the lines/sequences. We will be using these individual sequences later to create our matrix.
		filename = sys.argv[2]
		inputFile = open(filename)
		inputFile.next()
		
		for line in inputFile:
			if line.startswith(">"):
				sequenceList.append(stringSeq.strip())
				stringSeq = ""
				sequenceCounter += 1
			else:
				stringSeq += line
				
		sequenceList.append(stringSeq.strip())
		sequenceCounter += 1

		#Here we create our new distance matrix, passing in how many sequences we have total.
		distanceMatrix = numpy.zeros((sequenceCounter, sequenceCounter), dtype=float)
					
		# This will allow us to run through the pairs of sequences, for all sequences.
		for i in range(sequenceCounter):
			for j in range(counter, sequenceCounter):
				# This will see if we're running the same sequence, and if we're not, then we run our create matrix.
				if(not i == j):
					length_1 = len(sequenceList[i]) 
					length_2 = len(sequenceList[j])
					createMatrix(distanceMatrix, sequenceList, length_1, length_2, i, j)
					
			counter += 1
			
		print "\n\n", distanceMatrix	
		
if __name__ == "__main__":
    main(sys.argv[1:])
