import sys
import csv

def makeRef(str):
    result = ""
    blocks = str.split('\n')
    for i in blocks:
        if len(i) > 0 and i[0] != '>':
            result += i
    return result

#read in the reads
#takes in a string and returns a list of just the read sequences (split by >)
def makeReads(str):
    result = []
    blocks = str.split('\n')
    for i in blocks:
        if len(i) > 0 and i[0] != '>':
            result.append(i)
    return result
    
#compute the Hamming distance between two strings
def hammingDistance(a, b):
    if len(a) < len (b):
        length = len(a)
    else:
        length = len(b)
    result = 0
    for i in range(length):
        if a[i] != b[i]:
            result += 1
    return result

#create a dictionary where each key is a k-mer that maps to some place(s) in ref
def indexGenome(ref, k):
    perfectmap = {}
    for i in range(len(ref) - k):
        mer = ref[i:i + k]
        if mer not in perfectmap:
            perfectmap[mer] = [i]
        else:
            perfectmap[mer].append(i)
    return perfectmap
    
def possibleSubstitution(read, ref, matchIndex):
    if len(read) < len (ref):
        length = len(read)
    else:
        length = len(ref)
    subs = []
    for i in range(length):
        if read[i] != ref[i]:
            s = ">S" + str(matchIndex + i) + " " + ref[i] + " " + read[i]
            subs.append(s)
    return subs

def possibleIndel(read, ref, matchIndex):
    results = []
    #read up to where it stops matching
    i = 0
    while i < len(read) and read[i] == ref[i]:
        i += 1
    i -= 1
    readWithDeletion = read[:i] + ref[i] + read[i:]
    readWithInsertion = read[:i] + read[i + 1:]
    originalHD = hammingDistance(read, ref)
    delHD = hammingDistance(readWithDeletion, ref)
    inHD = hammingDistance(readWithInsertion, ref)
    #best delHD and inHD threshold: < 2
    if delHD < 2 and delHD < originalHD and delHD < inHD:
        #possible deletion
        indel = ">D" + str(matchIndex + i) + " " + ref[i]
        results.append(indel)
    elif inHD < 2 and inHD < originalHD and inHD < delHD:
        #possible insertion
        indel = ">I" + str(matchIndex + i) + " " + read[i]
        results.append(indel)
    return results

#maps read r to perfect match index. returns a list of mutations found in that read
def mapRead(index, r, ref):
    results = []
    if len(r) > 48:
        r = r[:48]
    if len(r) < 48:
        return results
    #divide into thirds
    first = r[:16]
    second = r[16:32]
    third = r[32:]
    #find perfect matches for each third
    #add corresponding positions for the start of the read to the list of matches
    matches = []
    if first in index:
        for i in index[first]:
            if i not in matches and i + 48 < len(ref):
                matches.append(i)
    if second in index:
        for i in index[second]:
            if i - 16 not in matches and i - 16 >= 0 and i + 32 < len(ref):
                matches.append(i - 16)
    if third in index:
        for i in index[third]:
            if i - 32 not in matches and i - 32 >= 0 and i + 16 < len(ref):
                matches.append(i - 32)
    #for each match find the hamming distance
    for m in matches:
        refmatch = ref[m : m + 48]
        #best: <= 3
        if hammingDistance(r, refmatch) <= 3:
            #most of read the same, possible substitution
            results += possibleSubstitution(r, refmatch, m)
        else:
            #possible indel
            results += possibleIndel(r, refmatch, m)
    return results

def isSub(mutation: str):
    if mutation[1:2] == "S":
        return True
    else:
        return False

def isIndel(mutation: str):
    if mutation[1:2] == "I" or mutation [1:2] == "D":
        return True
    else:
        return False

def filterErrors(freqs, subThreshold, indelThreshold):
    mutations = []
    for m in freqs:
        if (isSub(m) and freqs[m] >= subThreshold) or (isIndel(m) and freqs[m] >= indelThreshold):
            mutations.append(m)
    return mutations

#main

from sys import argv

#read in arguments
#reference genome
input1 = open(argv[1], "r")
strinput1 = str(input1.read())
#paired reads
input2 = open(argv[2], "r")
strinput2 = str(input2.read())

#read in and index the reference genome using 16mers
ref = makeRef(strinput1)
perfectMatches = indexGenome(ref, 16)

#read in and map the reads, storing possible mutations (treat paired reads as individual for now)
reads = makeReads(strinput2)
mutationFreqs = {}

#count how many times each possible mutation occurs
for r in reads:
    possibleMutations = mapRead(perfectMatches, r, ref)
    for m in possibleMutations:
        if m in mutationFreqs:
            mutationFreqs[m] += 1
        else:
            mutationFreqs[m] = 1

#filter out sequencing errors
subThreshold = 4
indelThreshold = 2
trueMutations = filterErrors(mutationFreqs, subThreshold, indelThreshold)


#write output to csv (run on Hoffman)
with open("predictions.csv", "w", newline = "") as f:
    writer = csv.writer(f)
    for m in trueMutations:
        writer.writerow([m])

