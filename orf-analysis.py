import csv
from numpy import median
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns 
from scipy import stats
from collections import Counter
import copy
import scipy
import random


orf_counter = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0}
norf_counter = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0, 10:0, 11:0, 12:0, 13:0, 14:0, 15:0, 16:0}
chromosomes = {'chrI': '01', 'chrII': '02', 'chrIII': '03', 'chrIV': '04', 'chrV': '05', 'chrVI': '06', 'chrVII': '07', 'chrVIII': '08', 'chrIX': '09', 'chrX': '10', 'chrXI': '11', 'chrXII': '12', 'chrXIII': '13', 'chrXIV': '14', 'chrXV': '15', 'chrXVI': '16'}
w_chromosomes = {'1': '01', '2': '02', '3': '03', '4': '04', '5': '05', '6': '06', '7': '07', '8': '08', '9': '09', '10': '10', '11': '11', '12': '12', '13': '13', '14': '14', '15': '15', '16': '16'}


def chromosome_distribution(list, dict):
    if list[6] == 'chrI': dict[1] += 1
    elif list[6] == 'chrII': dict[2] += 1
    elif list[6] == 'chrIII': dict[3] += 1
    elif list[6] == 'chrIV': dict[4] += 1
    elif list[6] == 'chrV': dict[5] += 1
    elif list[6] == 'chrVI': dict[6] += 1
    elif list[6] == 'chrVII': dict[7] += 1
    elif list[6] == 'chrVIII': dict[8] += 1
    elif list[6] == 'chrIX': dict[9] += 1
    elif list[6] == 'chrX': dict[10] += 1
    elif list[6] == 'chrXI': dict[11] += 1
    elif list[6] == 'chrXII': dict[12] += 1
    elif list[6] == 'chrXIII': dict[13] += 1
    elif list[6] == 'chrXIV': dict[14] += 1
    elif list[6] == 'chrXV': dict[15] += 1
    elif list[6] == 'chrXVI': dict[16] += 1

def reverse_complement(seq:str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

def summary_statistics():
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        lengths = []
        norf_count = 0
        corf_count = 0
        longnorf_count = 0
        extralongnorf_count = 0
        for row in reader:
            chromosome_distribution(row, orf_counter)
            if 'noncanonical' in row:
                chromosome_distribution(row, norf_counter)
                norf_count += 1
                coor1 = row[7]
                coor2 = row[8]
                length = abs(int(coor1)-int(coor2))
                lengths.append(length)
                if length > 300: longnorf_count += 1
                if length > 1000: extralongnorf_count += 1
            else:
                corf_count += 1
        chromosomes_percent_norfs = {key: round((norf_counter[key] / orf_counter[key] * 100) if orf_counter[key] > 0 else 0, 0) for key in orf_counter}
        print("Total cORFs:", corf_count, "Total nORFs:", norf_count)
        print("Average nORF Length:", sum(lengths)/len(lengths))
        print("Median nORF Length:", median(lengths))
        print("Longest nORF:", max(lengths) if lengths else 0)
        print("Long nORF Count:", longnorf_count)
        print("Extra Long nORF Count:", extralongnorf_count)
        print("Chromosome Distribution (Non-canonical) %:", chromosomes_percent_norfs)



def reverse_complement(seq:str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

def open_chromosome (chrom:int):
    with open (chrom) as f:
        lines = f.readlines()
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def query_rich(orf_id):
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        orf = ''
        for row in reader:
            if row[1] == orf_id or row[0] == orf_id:
                orf = row
        if orf == '': return ''
        chrom = chromosomes[orf[6]] 
        coor1 = int(orf[7])
        coor2 = int(orf[8])
        sequence = open_chromosome(f'SGD/{chrom}.fsa')
        excerpt = sequence[coor1-1:coor2]
        if orf[9] == '-': excerpt = reverse_complement(excerpt)
        #print(f"ORF ID: {orf[0]}\nChromosome: {chrom}\nCoordinates: {coor1}-{coor2}\nStrand: {orf[9]}\nType: {'noncanonical' if 'noncanonical' in orf else 'canonical'}")
        #print(f"Excerpt: {excerpt}")
        return excerpt
    
def query_wacholder(orf_id):
    with open('Wacholder_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        orf = ''
        for row in reader:
            if row[0] == orf_id:
                orf = row
        if orf == '': return "" 
        chrom = w_chromosomes[orf[1]]
        coor1 = int(orf[2])
        coor2 = int(orf[3])
        sequence = open_chromosome(f'SGD/{chrom}.fsa')
        excerpt = sequence[coor1-1:coor2]
        if orf[4] == '-': excerpt = reverse_complement(excerpt)
        return excerpt

def query_sequence(orf_id:str):
    if orf_id.startswith('chr'):
        return query_rich(orf_id)
    elif orf_id.startswith('orf'):
        return query_wacholder(orf_id)

def rich_to_wacholder(richname):
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if row[1] == richname:
                return row[0]
        return 'none'

def wacholder_survery():
    with open('Wacholder_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        lengths = []
        for row in reader:
            if row[6] == 'Transient':   
                coor1 = int(row[2])
                coor2 = int(row[3])
                length = abs(coor1-coor2)
                lengths.append(length)
    print("Wacholder ORF Count:", len(lengths))
    print("Average Wacholder ORF Length:", sum(lengths)/len(lengths))
    print("Median Wacholder ORF Length:", median(lengths))
    print("Longest Wacholder ORF Length:", max(lengths))
    return lengths

#wacholder_survery()

def longest_noncanonical_orf():
    longest_seq = ''
    longest_orf_id = ''
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if 'noncanonical' in row:
                chrom = chromosomes[row[6]]
                coor1 = int(row[7])
                coor2 = int(row[8])
                sequence = open_chromosome(f'SGD/{chrom}.fsa')
                excerpt = sequence[coor1-1:coor2]
                if row[9] == '-':
                    excerpt = reverse_complement(excerpt)
                if len(excerpt) > len(longest_seq):
                    longest_seq = excerpt
                    longest_orf_id = row[1]
    return longest_seq

def gc_content(sequence: str) -> float:
    count = 0
    for i in sequence:
        if i == 'G' or i == 'C':
            count += 1
        percentage = (count / len(sequence))
    return percentage

def codonize(sequence: str) -> list:
    return [sequence[i:i+3] for i in range(0, len(sequence)-len(sequence)%3, 3)]

def decodonize(codon_list: list) -> str:
    return ''.join(codon_list)

def start_codon (sequence: str) -> int:
    for i in range (len(sequence)-2):
        codon = sequence[i:i+3]
        if codon == 'ATG':
            return i

def first_stop_codon(sequence: str) -> int:
    from_atg = sequence[start_codon(sequence):]
    for i in range (len(from_atg)-2):
        codon = from_atg[i:i+3]
        if (i)%3 == 0 and codon in ('TAG', 'TAA', 'TGA'):
            return i
    return  

#print(first_stop_codon('CCCCATGCCCCCCCCCTAGGGGTAG')) #<- should return 12

def start_to_stop (sequence:str) -> str:
    if start_codon(sequence) == None or first_stop_codon(sequence) == None:
        return ""
    else: 
        from_atg = sequence[(start_codon(sequence)):]
        return from_atg[:(first_stop_codon(from_atg)+3)]

def first_start_last_stop (sequence, list):
    if start_to_stop (sequence) == "":
        return [sequence, list]
    else:
        to_atg = sequence[:start_codon(sequence)]
        first_orf = start_to_stop(sequence)
        excerpt_length = len(to_atg) + len(first_orf)
        list.append (first_orf)
        newseq = sequence[excerpt_length:]
        return first_start_last_stop (newseq, list)

def intron_survey():
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        total_lengths = []
        exon_counts = []
        for row in reader:
            if row[0] == 'NA':
                seq = query_sequence(row[1])
                orfs = first_start_last_stop(seq, [])[1]
                total_length = len(''.join(orfs))
                total_lengths.append(total_length)
                exon_count = len(orfs)
                exon_counts.append(exon_count)
        length_average = sum(total_lengths)/len(total_lengths)
        exon_average = sum(exon_counts)/len(exon_counts)
    print (f"Average Total Length: {length_average} bp. Median Total Length: {median(total_lengths)} bp.")
    print (f"Average Total Exon Count: {exon_average}. Median Total Exon Count: {median(exon_counts)}.")
    
    
#print(start_to_stop('CCCCATGCCCCCCCCCCTAGGGGTAG'))  #should return an empty string
#print(start_to_stop('CCCCATGCCCCCCCCCTAGGGGTAG'))  #should return ATGCCCCCCCCCTAG


def orf_check(id):
    fullstring = query_sequence(id)
    orfstring = start_to_stop(fullstring)
    return len(fullstring) - len(orfstring)

def orf_check_all ():
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        differences = []
        zerocount = 0
        for row in reader:
            id = row[1]
            fullstring = query_sequence(id)
            orfstring = start_to_stop(fullstring)
            diff = (len(fullstring) - len(orfstring))
            if diff == 0:
                zerocount +=1
            else: differences.append (diff)
    return differences


nucleotides = [ "A", "T", "C", "G"]

#output a random string of nucleotides of length n
def random_dna (n):
    str = ''
    for i in range (0, n-1):
        str += nucleotides[random.randint(0,3)]
    return str

#generate a single random ORF
def random_orf():
    str = ''
    for i in range (0, 2000):
        str += random.choices(nucleotides, weights = [0.31, 0.31, 0.19, 0.19], k = 1)[0]
        #check if the sequence contains a valid ORF
        if len(start_to_stop(str)) >= 1:
            return start_to_stop(str)
        elif i == 1000:
            return ""
    return

def average_orf(n):
    lens = []
    for j in range (0, n):
        seq = random_orf()
        lens.append(len(seq))
    print ("Median Length:", median (lens))
    print ("Mean length:", np.mean(lens))
    return lens

def intron_orf():
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        list = []
        for row in reader:
            if row[0] == 'NA' and 'noncanonical' in row:
                list.append(row[1])
    return list

def summary_output(numbers):
    if not numbers:
        print("No data provided.")
        return
    mean = np.mean(numbers)
    std = np.std(numbers)
    med = np.median(numbers)
    rng = max(numbers) - min(numbers)
    count = len(numbers)
    return(f"Mean: {mean}, Standard Deviation: {std}, Median: {med}, Range: {rng}, n = {count}")


def orf_length_plot():
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        corflengths = []
        norflengths = []
        rorflengths = average_orf(10000)
        for row in reader:
            if 'noncanonical' in row:
                coor1 = row[7]
                coor2 = row[8]
                length = abs(int(coor1)-int(coor2))
                norflengths.append(length)
            elif row[3] == 'canonical' and row[1] != 'NA':
                coor1 = row[7]
                coor2 = row[8]
                length = abs(int(coor1)-int(coor2))
                corflengths.append(length)
    print("nORF:", summary_output(norflengths))
    print("cORF:", summary_output(corflengths))
    print("rORF:", summary_output(rorflengths))
    plt.boxplot([rorflengths, norflengths, corflengths], labels=["gORFs", "nORFs", 'cORFs'])
    plt.yscale("log")
    plt.title('Box Plot of ORF Lengths')
    plt.xlabel('ORF Set')
    plt.ylabel('Length (bp)')
    plt.legend()
    plt.grid(True)
    plt.show()

def strong_mutants():
    with open('Deletion_Data.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        for i in range(3):
            next(reader)
        orfids = []
        for row in reader:
            fitness = float(row[2])
            qval = float(row[4])
            if fitness < 1 and qval <0.05:
                orfids.append(row[1])
        orfids = set(orfids)
    print(orfids)
    with open('Wacholder_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)   
        seqs = []
        lens = []
        for id in orfids:
            seq = query_sequence(id)
            seqs.append(seq)
            lens.append(len(seq))
    return (lens)

#print(strong_mutants())


def orf_length_stats():
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        corflengths = []
        norflengths = []
        rorflengths = average_orf(10000)
        for row in reader:
            if 'noncanonical' in row:
                coor1 = row[7]
                coor2 = row[8]
                length = abs(int(coor1)-int(coor2))
                norflengths.append(length)
            elif row[3] == 'canonical' and row[1] != 'NA':
                coor1 = row[7]
                coor2 = row[8]
                length = abs(int(coor1)-int(coor2))
                corflengths.append(length)
    norf_corf, pval1 = stats.ttest_ind(norflengths, corflengths)
    norf_rorf, pval2 = stats.ttest_ind(norflengths, rorflengths)
    corf_rorf, pval3 = stats.ttest_ind(corflengths, rorflengths)
    rorf_90 = np.percentile(rorflengths, 90)
    corf_10 = np.percentile(corflengths, 10)
    lower_10 = 0
    upper_90 = 0
    for i in norflengths:
        if i > corf_10: lower_10 +=1
        elif i < rorf_90: upper_90 +=1
    print("rorf90:", rorf_90, "corf90:", corf_10)
    print(upper_90, "nORFs shorter than 90th-percentile gORF")
    print(lower_10, "nORFs longer than the 10th-percentile cORF")
    print("norf shapiro:", stats.shapiro(np.log(norflengths)))
    print("gorf shapiro:", stats.shapiro((rorflengths)))
    print("corf shapiro:", stats.shapiro(np.log(corflengths)))
    print ("n/c:", norf_corf, pval1, "n/r:", norf_rorf, pval2, "c/r:", pval3)

orf_length_stats()

def codon_analyze(sequence):
    nucleotides = ["A", "T", "C", "G"]
    codons = [a + b + c for a in nucleotides for b in nucleotides for c in nucleotides]
    #print(codon_index)
    codon_counts = {codon: 0 for codon in codons}
    for i in codonize(sequence):
        codon_counts[i] +=1
    return codon_counts

#print(codon_analyze(query_sequence('orf169601')))

def random_orfs_codons(n):
    seqs = []
    for i in range(n):
        seqs.append(random_orf())
    rorf_dicts = []
    for rorf in seqs:
        rorf_dicts.append(codon_analyze(rorf))
    rorf_final_count = dict(sum((Counter(d) for d in rorf_dicts), Counter()))
    rorf_acid_sum = sum(rorf_final_count[f] for f in rorf_final_count)
    rorf_adjusted_count = codon_analyze('')
    for key in rorf_final_count:
        rorf_adjusted_count[key] = ((rorf_final_count[key] / rorf_acid_sum) * 100)
    return rorf_adjusted_count

def codon_survey():
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        #counting codons
        norfdicts = []
        corfdicts = []
        reader = csv.reader(csvfile)
        reader.__next__
        for row in reader:
            rowseq = query_sequence(row[1])
            if row[3] == 'noncanonical':
                norfdicts.append(codon_analyze(rowseq))
            elif row[3] == 'canonical' and row[0] != 'NA':
                corfdicts.append(codon_analyze(rowseq))
    norf_final_count = dict(sum((Counter(d) for d in norfdicts), Counter()))
    corf_final_count = dict(sum((Counter(d) for d in corfdicts), Counter()))
    #accounting for difference in total codon count
    norf_codon_sum = sum(norf_final_count[f] for f in norf_final_count)
    corf_codon_sum = sum(corf_final_count[f] for f in corf_final_count)
    norf_adjusted_count = codon_analyze('')
    corf_adjusted_count = codon_analyze('')
    for key in corf_final_count:
        norf_adjusted_count[key] = ((norf_final_count[key] / norf_codon_sum) * 100)
        corf_adjusted_count[key] = ((corf_final_count[key] / corf_codon_sum) * 100)
    #integrating random orfs
    rorf_adjusted_count = random_orfs_codons(10000)
    #finding differences
    diff_arit = codon_analyze('')
    diff_geo = codon_analyze('')
    for key in norf_final_count:
        diff_arit[key] = (corf_adjusted_count[key] - norf_adjusted_count[key])
        diff_geo[key] = (corf_adjusted_count[key] / norf_adjusted_count[key])
    #plottinga
    cats = []
    vals1 = []
    vals2 = []
    vals3 = []
    for key in norf_final_count:
        cats.append(key)
        vals1.append(norf_adjusted_count[key])
        vals2.append(corf_adjusted_count[key])
        vals3.append(rorf_adjusted_count[key])
    x = np.arange(len(cats))   
    width = 0.25
    plt.bar(x - width, vals2, width, label='Canonical ORFs')
    plt.bar(x, vals1, width, label='Noncanonical ORFs')
    plt.bar(x + width, vals3, width, label='Generated ORFs')    
    plt.xlabel('Codon')
    plt.ylabel('Codon Frequency (percent of total codons)')
    plt.title('Codon Usage of nORFs and cORFs')
    plt.legend()
    plt.xticks(rotation=90)
    plt.xticks(x, cats)
    plt.show()
    #print (f"Codon Counts in nORFs: {norf_final_count}.\nCodon Counts in cORFs: {corf_final_count} \nAdjusted Arithmetic Difference: {diff_arit}\nAdjusted Geometric Difference: {diff_geo}")

#codon_survey()

def conserved_orfs():
    with open('Wacholder_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        conserved_norfs = []
        for row in reader:
            if row[5] == 'X' and row[6] == 'Conserved':
                conserved_norfs.append(row[0])
    return(conserved_norfs)


the_thirteen = ['orf8709', 'orf11872', 'orf9589', 'orf7803', 'orf7969', 'orf66062', 'orf82574', 'orf94142', 'orf22494', 'orf16465', 'orf140368', 'orf180746', 'orf213157']

def thirteen_survey():
    lens = []
    phenos = []
    for i in the_thirteen:
        lens.append(len(query_sequence(i)))
    with open('Deletion_Data.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        for i in range(3): next(reader)
        for row in reader:
            if row[1] in the_thirteen:
                phenos.append(row[1])
                print(row)
    print (lens)
    print (phenos)

#thirteen_survey()
    
def long_orfsR(n):
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        list = []
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            if row[3] == 'noncanonical' and int(row[8]) - int(row[7]) > n:
                list.append(row[1])
    return list

#print(len(long_orfsR(1000)))

def long_orfsW(n):
    with open('Wacholder_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        list = []
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            if row[5] == 'X' and abs(int(row[3]) - int(row[2])) > n:
                list.append(row[0])
    return list

#print(len(long_orfsW(1000)))
    

def fix_my_shit():
    for i in long_orfsR(300):
        print(rich_to_wacholder(i), len(query_sequence(i)))
        
#fix_my_shit()

def fix_my_other_shit():
    with open('Wacholder_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        longlist = []
        next(reader)
        for row in reader:
            orfid = row[0]
            seq = query_sequence(orfid)
            if row[5] == 'X' and len(seq) > 1000:
                longlist.append(orfid)
                print(orfid, len(seq))
    print(longlist)

#fix_my_other_shit()

onekbp_norfs = ['orf467', 'orf8263', 'orf8452', 'orf7785', 'orf7924', 'orf32240', 'orf38826', 'orf40437', 'orf46489', 'orf68385', 'orf66887', 'orf65992', 'orf90558', 'orf86678', 'orf87953', 'orf84903', 'orf120976', 'orf120025', 'orf132862', 'orf129911', 'orf135775', 'orf131220', 'orf148259', 'orf149781', 'orf149897', 'orf147822', 'orf167341', 'orf189809', 'orf186549', 'orf210476', 'orf205083', 'orf205198', 'orf211102', 'orf211320', 'orf211989', 'orf226303', 'orf227274', 'orf231517', 'orf227570', 'orf223907', 'orf232165', 'orf228428', 'orf249248', 'orf248117', 'orf252456', 'orf252501', 'orf3967', 'orf22937', 'orf15806', 'orf21331', 'orf20866', 'orf17242', 'orf29803', 'orf29024', 'orf63683', 'orf51918', 'orf62652', 'orf59042', 'orf58640', 'orf74596', 'orf105969', 'orf97014', 'orf114364', 'orf114360', 'orf126910', 'orf126327', 'orf144992', 'orf141741', 'orf138366', 'orf159687', 'orf153359', 'orf153356', 'orf175359', 'orf200985', 'orf203622', 'orf215710', 'orf245401', 'orf244196', 'orf234838', 'orf234171', 'orf262587', 'orf264111']

codon_table = {
        'TTT':'Phe','TTC':'Phe','TTA':'Leu','TTG':'Leu',
        'TCT':'Ser','TCC':'Ser','TCA':'Ser','TCG':'Ser',
        'TAT':'Tyr','TAC':'Tyr','TAA':'Stop','TAG':'Stop',
        'TGT':'Cys','TGC':'Cys','TGA':'Stop','TGG':'Trp',
        'CTT':'Leu','CTC':'Leu','CTA':'Leu','CTG':'Leu',
        'CCT':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
        'CAT':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
        'CGT':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
        'ATT':'Ile','ATC':'Ile','ATA':'Ile','ATG':'Met',
        'ACT':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr',
        'AAT':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys',
        'AGT':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg',
        'GTT':'Val','GTC':'Val','GTA':'Val','GTG':'Val',
        'GCT':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
        'GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
        'GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly'
    }

aa_list = [
    'Ala', 'Arg', 'Asn', 'Asp', 'Cys',
    'Gln', 'Glu', 'Gly', 'His', 'Ile',
    'Leu', 'Lys', 'Met', 'Phe', 'Pro',
    'Ser', 'Thr', 'Trp', 'Tyr', 'Val',
    'Stop']

codon_aa_map = {aa: {codon: 0 for codon, a in codon_table.items() if a == aa}for aa in aa_list}

def acidize_sequence(seq):
    cd_seq = codonize(seq)
    aa_seq = []
    for codon in cd_seq:
        aa_seq.append(codon_table[codon])
    return aa_seq

def acid_analyze(seq):
    aa_seq = acidize_sequence(seq)
    aa_dict = {aa: 0 for aa in aa_list}
    for aa in aa_seq: aa_dict[aa] +=1
    return aa_dict

def random_orfs_acids(n):
    seqs = []
    for i in range(n):
        seqs.append(random_orf())
    rorf_dicts = []
    for rorf in seqs:
        rorf_dicts.append(acid_analyze(rorf))
    rorf_final_count = dict(sum((Counter(d) for d in rorf_dicts), Counter()))
    rorf_acid_sum = sum(rorf_final_count[f] for f in rorf_final_count)
    rorf_adjusted_count = acid_analyze('')
    for key in rorf_final_count:
        rorf_adjusted_count[key] = ((rorf_final_count[key] / rorf_acid_sum) * 100)
    return rorf_adjusted_count

def acid_survey():
    corf_dicts = []
    norf_dicts = []
    #counting acids
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            rowseq = query_sequence(row[1])
            if row[3] == 'noncanonical':
                norf_dicts.append(acid_analyze(rowseq))
            elif row[3] == 'canonical' and row[0] != 'NA':
                corf_dicts.append(acid_analyze(rowseq))
    #summing, adjusting
    norf_final_count = dict(sum((Counter(d) for d in norf_dicts), Counter()))
    corf_final_count = dict(sum((Counter(d) for d in corf_dicts), Counter()))
    norf_acid_sum = sum(norf_final_count[f] for f in norf_final_count)
    corf_acid_sum = sum(corf_final_count[f] for f in corf_final_count)
    norf_adjusted_count = acid_analyze('')
    corf_adjusted_count = acid_analyze('')
    for key in corf_final_count:
        norf_adjusted_count[key] = ((norf_final_count[key] / norf_acid_sum) * 100)
        corf_adjusted_count[key] = ((corf_final_count[key] / corf_acid_sum) * 100)
    #integrating random orfs
    rorf_adjusted_count = random_orfs_acids(10000)
    #plotting
    cats = []
    vals1 = []
    vals2 = []
    vals3 = []
    for key in norf_final_count:
        cats.append(key)
        vals1.append(norf_adjusted_count[key])
        vals2.append(corf_adjusted_count[key])
        vals3.append(rorf_adjusted_count[key])
    x = np.arange(len(cats))   
    width = 0.25
    plt.bar(x - width, vals2, width, label='Canonical ORFs', color = 'plum')  
    plt.bar(x, vals1, width, label='Noncanonical ORFs', color = 'teal')
    plt.bar(x + width, vals3, width, label='Generated ORFs', color = 'red')
    plt.xlabel('Amino Acid')
    plt.ylabel('Frequency (percent of total acids)')
    plt.title('Amino Acid Usage of nORFs and cORFs')
    plt.legend()
    plt.xticks(rotation=90)
    plt.xticks(x, cats)
    plt.show()

#acid_survey()

def meta_acid_analyze(seq):
    cd_seq = codonize (seq)
    map = copy.deepcopy(codon_aa_map)
    for codon in cd_seq:
        aa = codon_table[codon]
        map[aa][codon] +=1
    return map

#print(acid_analyze('ATGTACTTACTCATGCTATCTATCATGCAGTAA'))
#print(meta_acid_analyze('ATGTACTTACTCATGCTATCTATCATGCAGTAA'))

def total_variation_distance (x, y):
    return sum(abs(x[i]-y[i]) for i in range(len(x))) / 2   

# print(total_variation_distance([1, 2, 3, 4, 5], [5, 4, 3, 2, 1])) #should output 6 as TVD   

#print(stats.permutation_test(([0.1, 0.2, 0.5, 0.2], [0.2, 0.2, 0.3, 0.3]), total_variation_distance, permutation_type="samples", n_resamples=10000, vectorized=False, alternative="greater"))


#The actual codon bias function in all its glory
def meta_acid_survey():
    norf_maps = []
    corf_maps = []
    #counting acids
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            print(row[1])
            rowseq = query_sequence(row[1])
            rowmap = meta_acid_analyze(rowseq)
            if row[3] == 'noncanonical':
                norf_maps.append(rowmap)
            elif row[3] == 'canonical' and row[0] != 'NA':
                corf_maps.append(rowmap)
    #summing maps together
    amino_acids = norf_maps[0].keys()
    norf_final_count = {}
    corf_final_count = {}
    for acid in amino_acids:
        norf_sum = sum((Counter(sequence[acid]) for sequence in norf_maps), Counter())
        corf_sum = sum((Counter(sequence[acid]) for sequence in corf_maps), Counter())
        norf_final_count[acid] = dict(norf_sum)
        corf_final_count[acid] = dict(corf_sum)
    #generation of rORF data
    rorf_maps = []
    for i in range(100000):
        seq = random_orf()
        rorf_maps.append(meta_acid_analyze(seq))

    rorf_final_count = {}
    for acid in norf_final_count.keys():
        rorf_sum = sum((Counter(m[acid]) for m in rorf_maps), Counter())
        rorf_final_count[acid] = dict(rorf_sum)

    # data prep for plotting
    corf_diffsums = []
    norf_diffsums = []
    sig_corfs = []
    sig_norfs = []
    corf_x2vals = []
    norf_x2vals = []
    polycodonic_acids = []
    for acid, codons in norf_final_count.items():
        corfcodons = corf_final_count[acid]
        rorfcodons = rorf_final_count[acid]
        cats = []
        vals1raw = []
        vals2raw = []
        vals3raw = []
        for codon, count in codons.items():
            cats.append(codon)
            vals1raw.append(count)
            vals2raw.append(corfcodons.get(codon, 0))
            vals3raw.append(rorfcodons.get(codon, 0))
        rawsum1 = sum(vals1raw) or 1
        rawsum2 = sum(vals2raw) or 1
        rawsum3 = sum(vals3raw) or 1
        vals1adj = [v / rawsum1 for v in vals1raw]
        vals2adj = [v / rawsum2 for v in vals2raw]
        vals3adj = [v / rawsum3 for v in vals3raw]
        #plot for each acid
        '''width = 0.25
        x = np.arange(len(cats))
        plt.bar(x - width, vals2adj, width, label='Canonical ORFs', color='plum')
        plt.bar(x, vals1adj, width, label='Noncanonical ORFs', color='teal')
        plt.bar(x + width, vals3adj, width, label='Generated ORFs', color='red')
        plt.xlabel('Codon')
        plt.ylabel('Frequency (fraction of codons)')
        plt.title(f'Relative Codon Frequency in {acid}')
        plt.legend()
        plt.xticks(rotation=90)
        plt.xticks(x, cats)
        plt.show()'''
        #chi-square test + cramer's v
        corf_observed = np.array([vals2raw, vals3raw])
        corf_chi2, corf_pval, corf_dof, corf_expected = stats.chi2_contingency(corf_observed)
        norf_observed = np.array([vals1raw, vals3raw])
        norf_chi2, norf_pval, norf_dof, norf_expected = stats.chi2_contingency(norf_observed)
        corf_observed_sum = np.sum(corf_observed)
        norf_observed_sum = np.sum(norf_observed)
        corf_dims = min(corf_observed.shape)-1
        norf_dims = min(norf_observed.shape)-1
        corf_cv = np.sqrt((corf_chi2/corf_observed_sum) / corf_dims)
        norf_cv = np.sqrt((norf_chi2/norf_observed_sum) / norf_dims)
        corf_x2vals.append(corf_chi2)
        norf_x2vals.append(norf_chi2)
        #statistical analysis of differences from random
        print(acid, cats)
        if acid in ["Stop", "Trp", "Met"]: print(acid, "exlcuded")
        elif acid not in ["Stop", "Trp", "Met"]:
            polycodonic_acids.append(acid)
            c_r_diff = []
            n_r_diff = []
            for i in range(len(vals1adj)):
                c_r_diff.append(abs(vals2adj[i] - vals3adj[i]))
                n_r_diff.append(abs(vals1adj[i] - vals3adj[i]))
            c_r_summdiff = sum(c_r_diff)/2
            n_r_summdiff = sum(n_r_diff)/2
            print("CORFS TVD:", c_r_summdiff, "NORFS TVD:", n_r_summdiff)
            print("CORFS pval:", corf_pval, "df", corf_dof)
            print("NORFS pval:", norf_pval, "df", norf_dof)
            print ("CORFS Cramer's V:", corf_cv)
            print ("NORFS Cramer's V:", norf_cv)
            corf_diffsums.append(c_r_summdiff)
            norf_diffsums.append(n_r_summdiff)
            if corf_pval <0.001 and c_r_summdiff > 0.05: sig_corfs.append(acid)
            if norf_pval <0.001 and n_r_summdiff > 0.05: sig_norfs.append(acid)
    #plotting the stats
    print("significant cORFs:", len(sig_corfs), sig_corfs)
    print("significant nORFs:", len(sig_norfs), sig_norfs)
    print("t-test on the TVDs:", stats.ttest_rel(norf_diffsums, corf_diffsums))
    print("CORF mean TVD:", np.mean(corf_diffsums))
    print("NORF mean TVD:", np.mean(norf_diffsums))
    boxplot_data = [corf_diffsums, norf_diffsums]
    """plt.boxplot(boxplot_data, labels=['cORFs', 'nORFs'], meanline = True)
    plt.xlabel("ORF set")
    plt.ylabel("Codon variance from generated ORFs")
    plt.title("Codon Bias of Canonical and Noncanonical ORFs")
    plt.show()"""
    #plotting the chi2 values
    x = np.arange(len(polycodonic_acids))
    width = 0.4
    plt.bar(x - width, corf_diffsums, width, label='Canonical ORFs', color='plum')
    plt.bar(x, norf_diffsums, width, label='Noncanonical ORFs', color='teal')
    plt.xticks(rotation=90)
    plt.xlabel("Amino Acid")
    plt.ylabel("Total Variation Distance")
    plt.legend()
    plt.title("Codon Bias of Canonical and Noncanonical ORFs by Amino Acid")
    plt.xticks(x-0.2, polycodonic_acids)
    #plt.show()

#meta_acid_survey()

def check_noncanonical(orf):
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0] == orf:
                return row[3]

def check_noncanonical_set(srclist):
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        n_list = []
        c_list = []
        na_list = []
        for row in reader:
            if row[0] in srclist:
                orf = row[0]
                if row[3] == 'canonical': c_list.append(orf)
                elif row[3] == 'noncaonical': n_list.append(orf)
        for orf in srclist:
            if orf not in c_list or n_list: na_list.append(orf)
        print(f"corfs: {len(c_list)} norfs: {len(n_list)} n/a: {len(na_list)}")
        return n_list


def sort_by_length():
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        short_norfs = []
        medium_norfs = []
        long_norfs = []
        for row in reader:
            if 'noncanonical' in row:
                coor1 = row[7]
                coor2 = row[8]
                length = abs(int(coor1)-int(coor2))
                if length < 160: short_norfs.append(row[1])
                if length > 159: medium_norfs.append(row[1])
                if length > 1000: long_norfs.append(row[1])
    #print("Short nORFs:", len(short_norfs))
    #print("Medium nORFs:", len(medium_norfs))
    #print("Long nORFs:", len(long_norfs))
    return long_norfs
#sort_by_length()

def strong_mutants(x, y):
    with open('Deletion_Data.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        list = []
        for i in range(3):
            next(reader)
        qvals = []
        for row in reader:
            fitness = float(row[2])
            qval = float(row[4])
            if fitness < x and qval < y: list.append(row[1])
    return list
#from matplotlib_venn import venn2, venn3

def significant_sets():
    long_norfs = check_noncanonical_set(long_orfsW(300))
    purified_norfs = the_thirteen
    phenotype_norfs = strong_mutants(0.95, 0.05)
    sets = [long_norfs, purified_norfs, phenotype_norfs]
    sumlist = []
    for lst in sets:
        print(len(lst))
        sumlist.extend(lst)
    print("significant orf count:", len(set(sumlist)))
    duplicates = set(x for x in sumlist if sumlist.count(x) > 1)
    print("Duplicates:", duplicates, len(duplicates))
    for gene in duplicates:
        list = []
        for i in range(len(sets)):
            if gene in sets[i]: list.append(i)
        print(gene, list, len(query_sequence(gene)))
    return [duplicates, sets]


#significant_sets()

def test_my_shit():
    with open('Wacholder_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        conservedlist = []
        unclassifiedlist = []
        transientlist = []
        for row in reader:
            if row[0] in onekbp_norfs:
                if row[6] == 'Conserved': conservedlist.append(row[0])
                elif row[6] == "Unclassified": unclassifiedlist.append(row[0])
                else: transientlist.append(row[0])
        print("conserved count:", len(conservedlist))
        print("unclassified count:", len(unclassifiedlist))
        print("transient count:", len(transientlist))

#test_my_shit()

def test_my_other_shit():
    conservedlist = []
    unclassifiedlist = []
    unconservedlist = []
    with open('Wacholder_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        overlaps = significant_sets()[0]
        for row in reader: 
            if row[0] in overlaps:
                if row[6] == 'Conserved': conservedlist.append(row[0])
                elif row[6] == "Unclassified": unclassifiedlist.append(row[0])
                else: unconservedlist.append(row[0])
        print("conserved count:", len(conservedlist))
        print("unclassified count:", len(unclassifiedlist))
        print("unconserved count:", len(unconservedlist))

#test_my_other_shit()

def test_my_third_shit():
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        overlaps = significant_sets()[0]
        richoverlaps = []
        richnorfoverlaps = []
        for row in reader: 
            if row[0] in overlaps:
                richoverlaps.append(row[1])
                if row[3] != "canonical": richnorfoverlaps.append(row[0])
        print("richnorfoverlaps:", richnorfoverlaps)
        for orf in richnorfoverlaps:
            print(len(query_sequence(orf)))

#test_my_third_shit()

def test_my_fourth_shit(n):
    with open('Wacholder_ORFs.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        orflist = []
        for row in reader:
            orf = row[0]
            if len(query_sequence(orf)) > n:
                    orflist.append(orf)
    return orflist

#print(len(long_orfsR(300)))

import sys
import csv

csv.field_size_limit(sys.maxsize)

def coexpression_analysis():
    norflist = []
    corflist = [] 
    with open('Rich_ORFs.csv', newline='', encoding='utf-8') as richorfs:
        reader = csv.reader(richorfs)
        next(reader)
        for row in reader:
            if row[3] == 'noncanonical': norflist.append(row[1])
            elif row[3] == 'canonical': corflist.append(row[1])
    with open('Transposed_TOM.csv', newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        norf_rowsums = []
        corf_rowsums = []
        for row in reader: 
            orfname = row[0]
            if orfname in norflist: 
                rawvals = row[1:]
                floatvals = []
                for item in rawvals:
                    floatvals.append(float(item))
                norf_rowsums.append(sum(floatvals))
            elif orfname in corflist:
                rawvals = row[1:]
                floatvals = []
                for item in rawvals:
                    floatvals.append(float(item))
                corf_rowsums.append(sum(floatvals))
        print("corfs:", summary_output(corf_rowsums))
        print("norfs:", summary_output(norf_rowsums))
        connected_norf_count = 0
        for i in norf_rowsums: 
            if i > 192: connected_norf_count +=1
        print(connected_norf_count, "nORFs more connected than the average cORF")
    #plt.hist(norf_rowsums, bins=50)
    #plt.show()
    plt.hist(corf_rowsums, bins=100)
    plt.xlim(0, 1000)
    plt.xlabel("Total ORF Connectivity (arbitrary units)")
    plt.ylabel("Count")
    plt.title("Connectivity Distribution of Canonical ORFs")
    plt.show()

#coexpression_analysis()


             