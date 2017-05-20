import sys, math, operator
from random import randint
from collections import Counter
import re

def main(argv):
	print("\nGenerating Random Data...")

	generatedString = ""
	#motif = ["AACGCTTGCACCTTATTCGA","ACCGCGCCTGATTTGCGAA"]
	#motif = ["AACGCTTGCACCTT","ACCGCGCCTGATTT"]
	motif = ["AACGCTTG","ACCGCGCC"]
	#motif = ["--------","________"]

	indexList_1 = []
	indexList_2 = []
	index = 0

	finalString = ""
	print("\nLooking for motifs: " + str(motif) + "\n")
	addedLocations = []
	currentLine = ""
	count = 0


	while len(generatedString) <= 5000:

		currentRandom = randint(0,12)
		if currentRandom == 0 or currentRandom == 1 or currentRandom == 2:
			generatedString += "A"
			index+=1

		if currentRandom == 3 or currentRandom == 4 or currentRandom == 5:
			generatedString += "C"
			index+=1

		if currentRandom == 6 or currentRandom == 7 or currentRandom == 8:
			generatedString += "G"
			index+=1

		if currentRandom == 9 or currentRandom == 10 or currentRandom == 11:
			generatedString += "T"
			index+=1

		if currentRandom == 12:
			currentRandom2 = randint(0,200)
			if currentRandom2 == 5:
				motifToAdd = motif[randint(0,1)]
				generatedString += motifToAdd
				count += 1

			#addedLocations.append(index)
			#index+=len(motifToAdd)	
	
	#print(addedLocations)

	for character in generatedString:
		currentRandom = randint(0,900)

		if currentRandom == 10:
			motifToAdd = motif[randint(0,1)]
			finalString+=motifToAdd
			count+=1
			#addedLocations.append(currentRandom)
		finalString+=character

	#For every insert, I increased a counter. This will let us know at least how many motifs we should find. Depending
	#on how the motifs are split, we may have less or more.
	print("We should have around " + str(count) + " motifs (at least).\n")

	print("\n")
	print(finalString)
	print("\n")

	#Here, we seperate the entire string into smaller parts, in order to evaluate them. Some motifs will be cut off,
	#but this is on purpose, and is an edge case.	
	index = 0
	strandNum = 0
	for character in finalString:
		if index == 180:
			strandNum+=1
			print("Strand " + str(strandNum) + ": " + str(currentLine))
			for m in re.finditer(motif[0], currentLine):
				indexList_1.append([str(m.start()+1)])

			for m in re.finditer(motif[1], currentLine):
				indexList_2.append([str(m.start()+1)])
			print("Motif: " + str(motif[0]),indexList_1)
			print("Motif: " + str(motif[1]), indexList_2)
			print("\n\n")

			indexList_1 = []
			indexList_2 = []
			currentLine=""
			index=0

		else:
			currentLine+=character
			index+=1

	print("Strand " + str(strandNum+1) + ": " + str(currentLine))	

	for m in re.finditer(motif[0], currentLine):
		indexList_1.append([str(m.start()+1)])

	for m in re.finditer(motif[1], currentLine):
		indexList_2.append([str(m.start()+1)])
	print("Motif: " + str(motif[0]),indexList_1)
	print("Motif: " + str(motif[1]), indexList_2)
if __name__ == "__main__":
    main(sys.argv[0:])