import sys, os, argparse, re, csv

def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

def fasta(file_path):
    if os.path.isfile(file_path) and file_path.lower().endswith('.fasta'):
        return file_path
    else:
        raise argparse.ArgumentTypeError(f"readable_fasta: {file_path} is not a valid fasta file")

#takes in input directory of genomes and creates a dictionary where each key is the genome number and each value is the genome sequence
def extractGenomeSequences(directory):
    genome_sequences = {}  # Use a dictionary for easy access

    files = os.listdir(directory)
    for file_name in files:
        file_path = os.path.join(directory, file_name)
        with open(file_path, 'r') as file:
            genome_number = None
            sequence = ''
            for line in file:
                if line.startswith('>'):
                    match = re.search(r'Genome_Number_(\d+)', line)
                    if match:
                        genome_number = int(match.group(1))
                else:
                    sequence += line.strip()
            if genome_number is not None:
                genome_sequences[genome_number] = sequence
    return genome_sequences
    
#takes in fasta file with reads and creates a dictionary where each key is the read number and each value is the read
def makeReads(file_path):
    read_sequences = {}
    current_read_number = None
    current_sequence = ''

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>>'):
                if current_read_number is not None:
                    read_sequences[current_read_number] = current_sequence
                current_read_number = re.search(r'read_(\d+)', line).group(1)
                current_sequence = ''
            else:
                current_sequence += line.strip()

    if current_read_number is not None:
        read_sequences[current_read_number] = current_sequence

    return read_sequences

#create a dictionary where each key is a k-mer and each value is a list of the genomes it maps to
def indexGenomes(refs, k):
    perfectmap = {}
    for r in refs:
        for i in range(len(refs[r]) - k):
            mer = refs[r][i:i + k]
            if mer not in perfectmap:
                perfectmap[mer] = [r]
            else:
                perfectmap[mer].append(r)
    return perfectmap
    
#maps read r to perfect match index. returns a dict of the genomes that the read matches with how many times it maps to each one
def mapRead(index, r):
    matches = {}
    if len(r) > 45:
        r = r[:45]
    #divide into thirds
    first = r[:15]
    second = r[15:30]
    third = r[30:]
    #find perfect matches for each third
    #add corresponding genomes
    if first in index:
        for genome in index[first]:
            if genome not in matches:
                matches[genome] = 1
            else:
                matches[genome] += 1
    if second in index:
        for genome in index[second]:
            if genome not in matches:
                matches[genome] = 1
            else:
                matches[genome] += 1
    if third in index:
        for genome in index[third]:
            if genome not in matches:
                matches[genome] = 1
            else:
                matches[genome] += 1
    return matches

#return a narrowed down index with only the genomes appearing more than threshold times in the full index
#threshold = 0.4 for sample
def narrow(index, numGenomes, refs, reads, threshold):
    counts = [0] * numGenomes
    for r in reads:
        matches = mapRead(index, reads[r])
        for m in matches:
            counts[m] += 1
    abovethreshold = []
    #find which genomes pass the reads threshold
    for i in range(len(counts)):
        if counts[i]/sum(counts) > threshold:
            abovethreshold.append(i)
    #create new index with just the above-threshold genomes
    newdict = {}
    for genome in abovethreshold:
        for i in index:
            if genome in index[i]:
                if i not in newdict:
                    newdict[i] = [genome]
                else:
                    newdict[i].append(genome)
    return newdict
    
#given a read and an index, identify which genome in the index the read best maps to
def pickGenome(index, read):
    matches = mapRead(index, read)
    if len(matches) == 1:
        for m in matches:
            return m
    elif len(matches) > 1:
        bestmatch = -1
        maxcount = 0
        for genome in matches:
            if matches[genome] > maxcount:
                maxcount = matches[genome]
                bestmatch = genome
        return bestmatch #will always return last best match in case of tie
    #if can't find return 0
    return -1

#main
def main():

    parser = argparse.ArgumentParser(description='Map reads to reference genomes.')
    parser.add_argument('--referencedir', type=dir_path, required=True, help='Path to the reference genome directory')
    parser.add_argument('--reads', type=fasta, required=True, help='Path to the reads file')

    args = parser.parse_args()
    
    refs = extractGenomeSequences(args.referencedir)
    numGenomes = len(refs)
    reads = makeReads(args.reads)
    
    perfectMatches = indexGenomes(refs, 15)
    smallerMatches = narrow(perfectMatches, numGenomes, refs, reads, 0.007)

    results = {}
    for r in reads:
        matchingGenome = pickGenome(smallerMatches, reads[r])
        results[r] = matchingGenome
    
    #write output to csv
    with open("predictions.csv", "w", newline = "") as f:
        writer = csv.writer(f)
        for r in results:
            writer.writerow([f">read_{r} Genome_Number_{results[r]}"])
    

if __name__ == "__main__":
    main()
