###############################################################################

### Bio SeqIO is to parse a FASTA file, SeqUtils is for nt_search, which scans 
### a file for a desired string of nts, reverse_complement is used one time to
### search specifically for wc base pair potential in an Ori
### RNA is the Vienna RNA package that we used to fold the sequence surrounding
### a predicted Ori sequence

from Bio import SeqIO, SeqUtils
from Bio.Seq import reverse_complement
import RNA

###############################################################################

### function takes the combined FASTA containing all CRESS-DNA genomes, list of 
### Oris defined below, and a list of false positives that were manually determined
### to be false predictions after running this script

### Important note - we 'reorganize' the sequence by splitting it in half and appending
### the back half to the front because a lot of these genomes are amplified via
### RCR prior to sequencing, which makes the Ori span the begining and end of the 
### sequenced genome

### Once a genome has a prediction it is added to the predicted_genomes dictionary
### where it is excluded from search with the next Ori in the predefined list

def Ori_Seek(file, ori_list, false_ps):
    
    final = []
    predicted_genomes = {}
        
    for w in ori_list:
        for seq_record in SeqIO.parse(file, "fasta"):
            if seq_record.description not in false_ps:
                if seq_record.description not in predicted_genomes:
                    length = len(seq_record.seq)
                    reorg = ""
                    reorg += seq_record.seq[(round(length/2)):]
                    reorg += seq_record.seq[:(round(length/2))]
                    
                    OriMotifPositions = SeqUtils.nt_search(str(reorg), w)
                    
                    OriMotifPositions.pop(0)
                    
                    ProcessedOriMotifs = []
                    
                    for i in OriMotifPositions:
                        if reorg[i+7] == reverse_complement(reorg[i+3]):
                            temp = []
                            temp.append(i)
                            temp.append(str(reorg[i:i+9]))
                            temp.append(seq_record.description)
                            ProcessedOriMotifs.append(temp)
                    
                    Fold = []
                    
                    for i in ProcessedOriMotifs:
                        pStemLoop = str(reorg[i[0] - 15:i[0]+25])
                        if len(pStemLoop) == 40:
                            (ss, mfe) = RNA.fold(pStemLoop)
                            if mfe < -15.0:
                                ss_temp = ss
                                sl_temp = ss_temp.split(".")
                                for z in sl_temp:
                                    if z.count("(") >= 6:
                                        for z in sl_temp:
                                            if z.count(")") >= 6:
                                                if ss[14:21].count(".") >= 3:
                                                    if i[2] not in predicted_genomes:
                                                        temp = []
                                                        temp.append(i[2])
                                                        temp.append(i[1])
                                                        temp.append(mfe)
                                                        temp.append(ss)
                                                        Fold.append(temp)
                                                        final.append(Fold)
                                                        predicted_genomes[i[2]] = 1
                                                    
    return(final, predicted_genomes)

###############################################################################

### Same function as above, but does not use the reverse_complement function to 
### make a requirement for WC base pairing in the defined positions

def Ori_Seek_No_WC(file, predicted_genomes, ori_list, false_ps):
    
    final = []
    
    for w in ori_list:
        for seq_record in SeqIO.parse(file, "fasta"):
            if seq_record.description not in false_ps:
                if seq_record.description not in predicted_genomes:
                    length = len(seq_record.seq)
                    reorg = ""
                    reorg += seq_record.seq[(round(length/2)):]
                    reorg += seq_record.seq[:(round(length/2))]
                    
                    OriMotifPositions = SeqUtils.nt_search(str(reorg), w)
                    
                    OriMotifPositions.pop(0)
                    
                    ProcessedOriMotifs = []
                    
                    for i in OriMotifPositions:
                        temp = []
                        temp.append(i)
                        temp.append(str(reorg[i:i+9]))
                        temp.append(seq_record.description)
                        ProcessedOriMotifs.append(temp)
                    
                    Fold = []
                    
                    for i in ProcessedOriMotifs:
                        pStemLoop = str(reorg[i[0] - 15:i[0]+25])
                        if len(pStemLoop) == 40:
                            (ss, mfe) = RNA.fold(pStemLoop)
                            if mfe < -15.0:
                                ss_temp = ss
                                sl_temp = ss_temp.split(".")
                                for z in sl_temp:
                                    if z.count("(") >= 6:
                                        for z in sl_temp:
                                            if z.count(")") >= 6:
                                                if ss[14:21].count(".") >= 3:
                                                    if i[2] not in predicted_genomes:
                                                        temp = []
                                                        temp.append(i[2])
                                                        temp.append(i[1])
                                                        temp.append(mfe)
                                                        temp.append(ss)
                                                        Fold.append(temp)
                                                        final.append(Fold)
                                                        predicted_genomes[i[2]] = 1
                                                    
    return(final)

###############################################################################

### Bringing in the manually generated 'false positive' list and using it to
### exclude those genomes from search

false_pos_file = open("FalsePositives.txt", "r")

false_pos_dict = {}

for i in false_pos_file:
    if i in false_pos_dict:
        false_pos_dict[i.strip()] += 1
    else:
        false_pos_dict[i.strip()] = 1

ori_list = ["TAANATTNC", "NATNATTNC", "NAGNATTNC"]

wc_data, predicted = Ori_Seek("CRESSGenomes.FASTA", ori_list, false_pos_dict)

no_wc_data = Ori_Seek_No_WC("CRESSGenomes.FASTA", predicted, ori_list, false_pos_dict)

###############################################################################

### Function to count the number of WC-containing Oris

def seq_counts_wc(data, num, restriction):
    new_dict = {}
    for i in data:
        if i[0][1][num] != restriction:
            if i[0][1] in new_dict:
                new_dict[i[0][1]] += 1
            else:
                new_dict[i[0][1]] = 1
    
    return(new_dict)

### Function to count the number of Oris per search term

def seq_counts_species(data, num, requirement):
    new_dict = {}
    for i in data:
        if i[0][1][num] == requirement:
            if i[0][1] in new_dict:
                new_dict[i[0][1]] += 1
            else:
                new_dict[i[0][1]] = 1
    
    return(new_dict)

### Function to count the total number of Oris

def totals(dictionary, restriction):
    total = 0
    for i in dictionary:
        if i[3] != restriction:
            total += dictionary[i]
    
    return(total)

###############################################################################

###

wc_ori_dict = seq_counts_wc(wc_data, 3, "N")

no_wc_ori_dict = seq_counts_wc(no_wc_data, 3, "N")

wc_swaps_ori_dict = seq_counts_wc(wc_data, 3, "T")

swaps_count = totals(wc_swaps_ori_dict, "N")

###############################################################################

###

gemini_wc_counts = seq_counts_species(wc_data, 2, "A")

gemini_wc = totals(gemini_wc_counts, "N")

gemini_swaps = totals(gemini_wc_counts, "T")

gemini_no_wc_counts = seq_counts_species(no_wc_data, 2, "A")

gemini_no_wc = totals(gemini_no_wc_counts, "N")

###

circo_wc_counts = seq_counts_species(wc_data, 2, "T")

circo_wc = totals(circo_wc_counts, "N")

circo_swaps = totals(circo_wc_counts, "T")

circo_no_wc_counts = seq_counts_species(no_wc_data, 2, "T")

circo_no_wc = totals(circo_no_wc_counts, "N")

###

nano_wc_counts = seq_counts_species(wc_data, 2, "G")

nano_wc = totals(nano_wc_counts, "N")

nano_swaps = totals(nano_wc_counts, "T")

nano_no_wc_counts = seq_counts_species(no_wc_data, 2, "G")

nano_no_wc = totals(nano_no_wc_counts, "N")

###############################################################################

###

nonredundant_counts = {}

genomes = open("CRESSGenomes.FASTA", "r")     

for seq_record in SeqIO.parse(genomes, "fasta"):
    if seq_record.description not in nonredundant_counts:
        nonredundant_counts[seq_record.description] = 1
        
###############################################################################

###
     
percent_predictions = round(((len(wc_data) + len(no_wc_data))/(len(nonredundant_counts))* 100), 2)

percent_nonprediction = round((100 - percent_predictions), 2)

###############################################################################

###

percent_wc = round(((len(wc_data)/(len(wc_data) + len(no_wc_data)))*100), 2)

percent_no_wc = round(((len(no_wc_data)/(len(wc_data) + len(no_wc_data)))*100), 2)

percent_swaps = round(((swaps_count/(len(wc_data) + len(no_wc_data)))*100), 2)

###############################################################################

###

percent_gemini_predictions = round((gemini_wc + gemini_no_wc)/(len(wc_data) + len(no_wc_data))*100, 2)

percent_gemini_wc = round((gemini_wc)/(gemini_wc + gemini_no_wc)*100, 2)

percent_gemini_swap = round((gemini_swaps)/(gemini_wc + gemini_no_wc)*100, 2)

percent_gemini_no_wc = round((gemini_no_wc)/(gemini_wc + gemini_no_wc)*100, 2)

###############################################################################

###

percent_circo_predictions = round((circo_wc + circo_no_wc)/(len(wc_data) + len(no_wc_data))*100, 2)

percent_circo_wc = round((circo_wc)/(circo_wc + circo_no_wc)*100, 2)

percent_circo_swap = round((circo_swaps)/(circo_wc + circo_no_wc)*100, 2)

percent_circo_no_wc = round((circo_no_wc)/(circo_wc + circo_no_wc)*100, 2)

###############################################################################

###

percent_nano_predictions = round((nano_wc + nano_no_wc)/(len(wc_data) + len(no_wc_data))*100, 2)

percent_nano_wc = round((nano_wc)/(nano_wc + nano_no_wc)*100, 2)

percent_nano_swap = round((nano_swaps)/(nano_wc + nano_no_wc)*100, 2)

percent_nano_no_wc = round((nano_no_wc)/(nano_wc + nano_no_wc)*100, 2)

###############################################################################

###

export_predictions_all = open("Ori_Predictions_ALL.txt", "w")
export_predictions_wc = open("Ori_Predictions_WC.txt", "w")
export_predictions_wc_swaps = open("Ori_Predictions_WC_swaps.txt", "w")
export_predictions_no_wc = open("Ori_Predictions_NO_WC.txt", "w")
export_counts = open("Ori_Predictions_Counts.txt", "w")
export_supplemental = open("Ori_Predictions_Supplemental_Information.txt", "w")

### File containing all predictions

print("Genome Description" +"\t" + "Putative Ori" + "\t" + "Putative Structure" + "\t" + "Putative Folding Energy", file = export_predictions_all)

for i in wc_data:
    print(str(i[0][0]) + "\t" + str(i[0][1]) + "\t" + str(i[0][3]) + "\t" + str(i[0][2]), file = export_predictions_all)
    
for i in no_wc_data:
    print(str(i[0][0]) + "\t" + str(i[0][1]) + "\t" + str(i[0][3]) + "\t" + str(i[0][2]), file = export_predictions_all)

### File containing only WC-containing Ori predictions

print("Genome Description" +"\t" + "Putative Ori" + "\t" + "Putative Structure" + "\t" + "Putative Folding Energy", file = export_predictions_wc)

for i in wc_data:
    print(str(i[0][0]) + "\t" + str(i[0][1]) + "\t" + str(i[0][3]) + "\t" + str(i[0][2]), file = export_predictions_wc)

### File containing only WC-Swaps Ori predictions

print("Genome Description" +"\t" + "Putative Ori" + "\t" + "Putative Structure" + "\t" + "Putative Folding Energy", file = export_predictions_wc_swaps)

for i in wc_data:
    if i[0][1][3] != "T":
        print(str(i[0][0]) + "\t" + str(i[0][1]) + "\t" + str(i[0][3]) + "\t" + str(i[0][2]), file = export_predictions_wc_swaps)

### File containing only non-WC Ori predictions

print("Genome Description" +"\t" + "Putative Ori" + "\t" + "Putative Structure" + "\t" + "Putative Folding Energy", file = export_predictions_no_wc)

for i in no_wc_data:
    print(str(i[0][0]) + "\t" + str(i[0][1]) + "\t" + str(i[0][3]) + "\t" + str(i[0][2]), file = export_predictions_no_wc)

### File containing counts

print("Identity" + "\t" + "Total Counts" + "\t" + "WC Counts" + "\t" + "WC Swap Counts" + "\t" + "non-WC Counts", file = export_counts)

print("All Search Terms" + "\t" + str((len(wc_data) + len(no_wc_data))) + "\t" + str(len(wc_data)) + "\t" + str(swaps_count) + "\t" + str(len(no_wc_data)), file = export_counts)

print("Gemini Search Term" + "\t" + str((gemini_wc + gemini_no_wc)) + "\t" + str(gemini_wc) + "\t" + str(gemini_swaps) + "\t" + str(gemini_no_wc), file = export_counts)

print("Circo Search Term" + "\t" + str((circo_wc + circo_no_wc)) + "\t" + str(circo_wc) + "\t" + str(circo_swaps) + "\t" + str(circo_no_wc), file = export_counts)

print("Nano Search Term" + "\t" + str((nano_wc + nano_no_wc)) + "\t" + str(nano_wc) + "\t" + str(nano_swaps) + "\t" + str(nano_no_wc), file = export_counts)

### File containing supplemental information

print("Identity" + "\t" + "% Predicted" + "\t" + "% WC" + "\t" + "% WC Swap" + "\t" + "% non-WC", file = export_supplemental)

print("All Search Terms" + "\t" + str(percent_predictions) + "\t" + str(percent_wc) + "\t" + str(percent_swaps) + "\t" + str(percent_no_wc), file = export_supplemental)

print("Gemini Search Term" + "\t" + str(percent_gemini_predictions) + "\t" + str(percent_gemini_wc) + "\t" + str(percent_gemini_swap) + "\t" + str(percent_gemini_no_wc), file = export_supplemental)

print("Circo Search Term" + "\t" + str(percent_circo_predictions) + "\t" + str(percent_circo_wc) + "\t" + str(percent_circo_swap) + "\t" + str(percent_circo_no_wc), file = export_supplemental)

print("Nano Search Term" + "\t" + str(percent_nano_predictions) + "\t" + str(percent_nano_wc) + "\t" + str(percent_nano_swap) + "\t" + str(percent_nano_no_wc), file = export_supplemental)

export_predictions_all.close()
export_predictions_wc.close()
export_predictions_wc_swaps.close()
export_predictions_no_wc.close()
export_counts.close()
export_supplemental.close()

###############################################################################

###

putative_oris = open("AllPutativeOris.txt", "w")

print("Search Term" + "\t" + "Putative Ori" + "\t" + "Count", file = putative_oris)

for i in gemini_wc_counts:
    print("Gemini_WC" + "\t" + str(i) + "\t" + str(gemini_wc_counts[i]), file = putative_oris)

for i in gemini_no_wc_counts:
    print("Gemini_no_WC" + "\t" + str(i) + "\t" + str(gemini_no_wc_counts[i]), file = putative_oris)

for i in circo_wc_counts:
    print("Circo_WC" + "\t" + str(i) + "\t" + str(circo_wc_counts[i]), file = putative_oris)

for i in circo_no_wc_counts:
    print("Circo_no_WC" + "\t" + str(i) + "\t" + str(circo_no_wc_counts[i]), file = putative_oris)

for i in nano_wc_counts:
    print("Nano_WC" + "\t" + str(i) + "\t" + str(nano_wc_counts[i]), file = putative_oris)

for i in nano_no_wc_counts:
    print("Nano_no_WC" + "\t" + str(i) + "\t" + str(nano_no_wc_counts[i]), file = putative_oris)

putative_oris.close()

###############################################################################
