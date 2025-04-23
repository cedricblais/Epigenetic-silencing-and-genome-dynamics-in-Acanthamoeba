#python version: 3.11.9 ('biopython')

"""When looking for which viral genes overlap with unaligned regions in the mummer alignment, a few additional files are worth creating.

The primary aim of this script is to concatenate all overlapping regions in mummer4's show-diff file. In particular, this script does so in a way that accounts
for the fact that only the gaps on the edges of translocations, and not the transcolations themselves, are listed in the show-diff file (under the name SEQ). This
script uses the coordinates of consecutive SEQ features to infer the approximate coordinates of translocations. This way, concatenated divergent regions are not 
broken up into several regions because of short alignments corresponding to translocations.

The output is a gff-like file with the coordinates of each concatenated divergent region, which specifies in the "TYPE column whether this regions fully overlaps
with a viral gene. Also of interest is "

In addition to a file listing the concatenated divergence regions, the script produces the following useful files:
- a gff file of show-diff features which includes translocations (labelled as INTER-SEQ)
- a bed file of show-diff features which includes unaligned contigs (labelled as NO_HIT)
- a gff file of genes indicating which are conserved, degraded or unconserved based on whether they didn't overlap, artially overlapped, or fully overlapped with gaps
(note that not all show-diff features are gaps! so this is not the same as saying they are part of concatenated gap)

The script requires the following inputs:
- the .unref file produced by mumer4
- a .gff file where viral genes have the type 'virus_gene', and all other genes have the type 'gene'
- the mummer4 .out.mdelta file
- the mummer4 .out.1delta file
-the mummer 4 show-diff file
"""




import pandas as pd

from io import StringIO

from collections import Counter

my_types = {
    "[SEQ]": "string",
    "[TYPE]": "string",
    "[S1]": "Int64",
    "[E1]": "Int64",
    "[LEN 1]": "Int64",
}


"""SET VARIABLES"""

output_prefix = 'Neff_vs_C3'

output_directory = '/Users/cedricblais/Documents/PhD_work/paper_work/data_files/'

#What is the maximum size that a short local alignment can have to be included in the concatenation
threshold = 30000

#What is the distance allowed between regions before they are concatenated (distance + 1). E.g.: 1 will concatenate regions that are flush ([1,4][5,7] -> [1,7]), 2 will allow for for 1 bp of separation ([6,10][12,15] -> [6,15])
overlap_value = 1

#list of types included in the analysis
types_list = ['BRK', 'DUP', 'INTER-SEQ', 'INV', 'SEQ', 'GAP', 'JMP']

"""Get files"""

#open .unref file, which tells us which contigs have no alignment to the other genome and should be treated as gaps
#Neff: find unref in mummer directory
#C3: find unref in mummer directory
unref = pd.read_csv("out.unref", sep='\t', na_values=["-"], header = None)

#print(unref)
unref.columns = ['SEQ', 'FEAT', 'S1', 'E1', 'LEN']

#open virus gff file; can be found in the directory Custom_python_input
#Neff: Neff_genes_info.gff
#C3: C3_genes_info.gff 
virus_df = pd.read_csv("Neff_genes_info.gff", sep='\t', na_values=["-"], usecols=range(9), header = None)

columns_list = ['SEQ', 'SOURCE', 'TYPE', 'S1', 'E1', 'SCORE', 'STRAND', 'FRAME', 'ATTRIBUTE']

virus_df.columns = columns_list

#We set aside a version that will not be purged of non-viral things
all_genes = virus_df

#print(all_genes)

#open mdelta and1delta files, which we use to get length of scaffolds


#Neff: find out.mdelta in mummer directory
#C3: find out.mdelta in mummer directory
mdelta = pd.read_csv("out.mdelta", sep=' ', na_values=["-"], names=range(7), header = None)





#We retrieve header lines, strip them of '>', and remove redudant lines

mdelta = mdelta[mdelta[0].str.contains(">", na=False)]

mdelta[0] = mdelta[0].str.replace(r'>', '', regex=True)

first_row_m = mdelta.groupby(0).first()

#Neff: find out.1delta in mummer directory
#C3: find out.1delta in mummer directory
delta1 = pd.read_csv("out.1delta", sep=' ', na_values=["-"], names=range(7), header = None)

delta1 = delta1[delta1[0].str.contains(">", na=False)]

delta1[0] = delta1[0].str.replace(r'>', '', regex=True)

first_row_1 = delta1.groupby(0).first()

#We merge the two lists of lengths and create a dictionary for chromosome lengths
to_merge = [first_row_1, first_row_m]

merged_info = pd.concat(to_merge)

#print(merged_info)

merged_info.columns = ['1', '2', '3', '4', '5', '6']

merged_info.set_index(merged_info['1'])
length_dict = merged_info['2'].to_dict()


#Open mummer diff file, skipping the first four lines (non tsv info)
#/Users/cedricblais/Documents/PhD_work/mummer/Neff_v_c3_alignments/neff_c3_defaults_breaklen1000/Colp_Neff_vs_C3_defaults_breaklen1000.nucmer.delta.show-diff
#/Users/cedricblais/Documents/PhD_work/mummer/C3_v_Neff_alignments/Colp_C3_vs_Neff_defaults_breaklen_1000.nucmer.delta.show-diff
diff_input = "/Users/cedricblais/Documents/PhD_work/mummer/C3_v_Neff_alignments/Colp_C3_vs_Neff_defaults_breaklen_1000.nucmer.delta.show-diff"
mummer_diff = pd.read_csv(diff_input, sep='\t', na_values=["-"], skiprows=[0, 1, 2, 3], names=range(7), header=None)

mummer_diff.columns = ['SEQ', 'TYPE', 'S1', 'E1', 'LEN', 'SUP1', 'SUP2']


#For every feature in the file, which check whether start is smaller than end. If not, we reverse them (mummer allows for reversed coordinates and negative lengths in some cases)

for index, row in mummer_diff.iterrows():
    if row['S1'] > row['E1']:
        mummer_diff.at[index, 'S1'], mummer_diff.at[index, 'E1'] = row['E1'], row['S1']

#We get a list of all scaffolds in the file, and create a dictionary with each scaffold as a key
scaffolds = mummer_diff['SEQ'].drop_duplicates().to_list()

dict_scaff = dict.fromkeys(scaffolds, [])


#We separate our diff data frame by scaffold using the directory
for scaffold in dict_scaff:
    dict_scaff[scaffold] = mummer_diff[(mummer_diff['SEQ'] == str(scaffold))]

#Outputs will be stored here
out_list = ''

#We will store all outputted diffs here
diff = ['SEQ', 'TYPE', 'S1', 'E1', 'LEN', 'SUP1', 'SUP2']
CONTIG_diff_store = pd.DataFrame(columns=diff)


#We iterate through every contig represented in the file:
for scaffold in dict_scaff:
    #print(scaffold)
    scaff_list = []

    #We get all lines containing SEQ
    SEQ = dict_scaff[scaffold][(dict_scaff[scaffold]['TYPE'] == "SEQ")]

    #We set a variable for the dataframe containing only diffs in this particular contig
    CONTIG_diff = dict_scaff[scaffold]

    #We check if this contig contains SEQ: if it does, we use it to generate our desired small translocation coordinates
    if SEQ.empty:
        pass
    else:
        
        SEQ = SEQ.reset_index(drop=True)

        new = SEQ['E1'].shift(1)

        start_align = SEQ['E1'].shift(1) + 1
        end_align = SEQ['S1'] - 1

        seq = SEQ['SUP1']
        prev_seq = SEQ['SUP2'].shift(1)

        SEQ['AL_START'] = start_align

        SEQ['AL_END'] = end_align

        SEQ['SEQUENCE'] = SEQ['SUP1']

        SEQ['DISTANCE'] = SEQ['S1']- SEQ['E1'].shift(1)

        #We correct the coordinates of the first entry in each dataframe

        SEQ.loc[0, 'AL_START'] = 1

        SEQ.loc[0, 'DISTANCE'] = SEQ['AL_END'].loc[0]

        #We write the output to a new dataframe

        columns = ['SEQ', 'TYPE', 'S1', 'E1', 'LEN', 'SUP1', 'SUP2']

        output_df = pd.DataFrame(columns=columns)

        #print(output_df)

        output_df['SEQ'] = SEQ['SEQ']
        output_df['TYPE'] = "INTER-SEQ"
        output_df['S1'] = SEQ['AL_START']
        output_df['E1'] = SEQ['AL_END']
        output_df['LEN'] = SEQ['DISTANCE']
        output_df['SUP1'] = SEQ['SUP1']
        output_df['SUP2'] = 'NA'
        

        #We add a final line, for the aligning region after the final SEQ

        if length_dict[SEQ['SEQ'].loc[0]] > SEQ['E1'].loc[0]:

            #print(output_df)
            last_line = [SEQ['SEQ'].loc[0], "INTER-SEQ", SEQ['E1'].iloc[-1] + 1, length_dict[SEQ['SEQ'].loc[0]], length_dict[SEQ['SEQ'].loc[0]] - SEQ['E1'].iloc[-1], SEQ['SUP2'].iloc[-1], 'NA']
            output_df.loc[len(output_df)] = last_line

        #We create a df with the total lengths of alignments to each contig:
   
        alignment_sums_df = output_df.groupby('SUP1', as_index=True).sum('LEN')

        #We setup a cutoff as one third the total length of the query contig

        cutoff = length_dict[SEQ['SEQ'].loc[0]]/3



        #We use this cutoff to determine which contigs are considered to be part of the true contiguous alignment:
       
        alignment_sums_df = alignment_sums_df[(alignment_sums_df['LEN'] > cutoff)]

        exclusions = alignment_sums_df.index.values.tolist()


        #We filter each list based on the exclusions, keeping only alignments that do not hit the main homologous contig:

        filtered_output_df = output_df[~output_df['SUP1'].isin(exclusions)]

        
        #We further filter the remaining alignments to non-contiguous contigs by size

        filtered_output_df = filtered_output_df[(filtered_output_df['LEN'] < threshold)]


        #We add the identified alignments to the mummer diff

        add_to_diff = [filtered_output_df, CONTIG_diff]

        CONTIG_diff = pd.concat(add_to_diff)


    
    #Because of some overlap between SEQs, we make sure the start is the before the end
    for index, row in CONTIG_diff.iterrows():
        if row['S1'] > row['E1']:
            #print(str(row[2]) +"\t"+ str(row[3]))
            CONTIG_diff.at[index, 'S1'], CONTIG_diff.at[index, 'E1'] = row['E1'], row['S1']

    
    add_to_store = [CONTIG_diff_store, CONTIG_diff]

    CONTIG_diff_store = pd.concat(add_to_store)


    #we sort for good measure
    CONTIG_diff = CONTIG_diff.sort_values(by=['SEQ', 'S1'])

    list_of_coords = []

    for index, rows in CONTIG_diff.iterrows():
        coord_list = [rows.S1, rows.E1] 
        list_of_coords.append(coord_list)

    if len(list_of_coords) > 1000:
        print('DANGER MIGHT BREAK RECURSION LIMIT')

    #We collapse all overlapping lists into one list
    def recursive_merge(inter, start_index = 0):
        for i in range(start_index, len(inter) - 1):
            if inter[i][1] >= (inter[i+1][0] - overlap_value):
                new_start = inter[i][0]
                if inter[i+1][1] >= inter[i][1]:
                    new_end = inter[i+1][1]
                else:
                    new_end = inter[i][1]
                inter[i] = [new_start, new_end]
                del inter[i+1]
                return recursive_merge(inter.copy(), start_index=i)
        return inter    

    sorted_on_start = sorted(list_of_coords)
    merged = recursive_merge(sorted_on_start.copy())

    ranges_concat_gaps = []

    for list in merged:
        list.insert(0, scaffold)
        scaff_list.append(list)
        val = (list[1], list[2], list[0])

        overlapping_features = CONTIG_diff[(CONTIG_diff['S1'] >= list[1]) & (CONTIG_diff['S1'] <= list[2])]

        features = overlapping_features['TYPE'].tolist()

        features_number = len(overlapping_features) - 1

    
        length_kb = round(float(list[2] - (list[1] - 1))/float(1000), 3)

        complexity = round(features_number/length_kb, 2)

        list.append(str(length_kb)+"kb")
        list.append(features_number)
        list.append(complexity)

        counter_out = (Counter(features))

        sorted_dict = {key: value for key, value in sorted(counter_out.items())}

        out_features = 'DIFF='

        for entry in sorted_dict:
            out_features += str(entry)+":"+str(sorted_dict[entry])+";"

        out_features = out_features.strip(";")

        id = 'ID='+str(scaffold) +"_diff_"+ str(merged.index(list)+1)

        list.append(id)

        #print(list)

        line = "\t".join(str(item) for item in list)

        #print(line)

        gff_line = str(list[0])+"\t"+"mummer4_custom"+"\t"+"diff_concat"+"\t"+str(list[1])+"\t"+str(list[2])+"\t"+str(complexity)+"\t"+"."+"\t"+"."+"\t"+str(id)+";"+str(out_features)+"\n"

        out_list += str(gff_line)


#We convert the string gff into a new dataframe

out_list_IO = StringIO(out_list)

gff = pd.read_csv(out_list_IO, sep='\t', na_values=["-"], header = None)

columns_list = ['SEQ', 'SOURCE', 'TYPE', 'S1', 'E1', 'SCORE', 'STRAND', 'FRAME', 'ATTRIBUTE']

gff.columns = columns_list


#We now list all contigs not aligning to the reference as concat gaps

columns = ['SEQ', 'SOURCE', 'TYPE', 'S1', 'E1', 'SCORE', 'STRAND', 'FRAME', 'ATTRIBUTE']

unref_gff = pd.DataFrame(columns=columns_list)


unref_gff['SEQ'] = unref['SEQ']
unref_gff['SOURCE'] = 'mummer4_custom'
unref_gff['TYPE'] = 'diff_concat'
unref_gff['S1'] = unref['S1']
unref_gff['E1'] = unref['E1']
unref_gff['SCORE'] = 0
unref_gff['STRAND'] = '.'
unref_gff['FRAME'] = '.'
unref_gff['ATTRIBUTE'] = "ID="+unref['SEQ'].astype(str)+";DIFF=full_unaligned_contig"

to_merge = [gff, unref_gff]

#The following dataframe contains the information for all concatenated gaps:

all_concat_gaps_gff = pd.concat(to_merge)

"""We now determine which of these concatenated gaps contain viruses"""

#We make sure that index matches the number of each line:

all_concat_gaps_gff = all_concat_gaps_gff.reset_index()


columns_list = ['SEQ', 'SOURCE', 'TYPE', 'S1', 'E1', 'SCORE', 'STRAND', 'FRAME', 'ATTRIBUTE']

virus_df.columns = columns_list

virus_df = virus_df[(virus_df['TYPE'] == 'virus_gene')]


#We loop through every gap and check to see whether virus genes fall in it
for index, row in all_concat_gaps_gff.iterrows():
    #print(row['S1'])
    start = row['S1']
    end = row['E1']
    chrom = row['SEQ']
    #We filter for viral candidates fully overlapping those regions
    get_virus = virus_df[(virus_df['S1'] >= start) & (virus_df['S1'] <= end) & (virus_df['E1'] >= start) & (virus_df['E1'] <= end) & (virus_df['SEQ'].isin([chrom]))]
    get_virus = get_virus.reset_index(drop=True)

    #We filter for viral candidates which have any overlap with those regions
    get_virus_inc = virus_df[((((virus_df['S1'] >= start) & (virus_df['S1'] <= end)) | ((virus_df['E1'] >= start) & (virus_df['E1'] <= end))) | ((virus_df['S1'] <= start) & (virus_df['E1'] >= end))) & (virus_df['SEQ'].isin([chrom]))]
    get_virus_inc = get_virus_inc.reset_index(drop=True)

    full_overlap = get_virus['ATTRIBUTE'].tolist()
    partial_overlap = get_virus_inc['ATTRIBUTE'].tolist()

    Virus_number = len(full_overlap)
    Partial_virus_number = len(partial_overlap) - len(full_overlap)

    #we make list of overlapping ids, partial and full
    full_ids = [item.strip("ID=").strip(";") for item in full_overlap]
    partial_ids = []
    
    for entry in partial_overlap:
        if str(entry) not in str(full_ids):
            partial_ids.append(entry.strip("ID=").strip(";"))

    #If there are viral candidates completely fitting with the region, it is labelled as viral:
    if len(full_overlap) > 0:
        all_concat_gaps_gff.loc[index, 'TYPE'] = 'viral_diff_concat'
    
    if len(partial_overlap) > 0:
        all_concat_gaps_gff.loc[index, 'ATTRIBUTE'] = all_concat_gaps_gff.loc[index, 'ATTRIBUTE']+";FULL_VIR:"+str(len(full_overlap))+";PART_VIR:"+str(len(partial_overlap) - len(full_overlap))

        #unref_gff['ATTRIBUTE'] = "ID="+unref['SEQ'].astype(str)+";DIFF=full_unaligned_contig"

#print(all_concat_gaps_gff)

header = ['SEQ', 'SOURCE', 'TYPE', 'S1', 'E1', 'SCORE', 'STRAND', 'FRAME', 'ATTRIBUTE']

#print('THE DIFF')
#print(CONTIG_diff)
#print(CONTIG_diff_store)

CONTIG_diff_gff = pd.DataFrame(columns=header)

CONTIG_diff_gff['SEQ'] = CONTIG_diff_store['SEQ']
CONTIG_diff_gff['SOURCE'] = 'mummer4'
CONTIG_diff_gff['TYPE'] = CONTIG_diff_store['TYPE']
CONTIG_diff_gff['S1'] = CONTIG_diff_store['S1'].astype(int)
CONTIG_diff_gff['E1'] = CONTIG_diff_store['E1'].astype(int)
CONTIG_diff_gff['SCORE'] = '.'
CONTIG_diff_gff['STRAND'] = '.'
CONTIG_diff_gff['FRAME'] = '.'
CONTIG_diff_gff['ATTRIBUTE'] = '.'
all_concat_gaps_gff['S1'] = all_concat_gaps_gff['S1'].astype(int)
all_concat_gaps_gff['E1'] = all_concat_gaps_gff['E1'].astype(int)

#We add the unaligned contigs to the diff file gff


add_unaligned_to_diff = [CONTIG_diff_gff, unref_gff]

CONTIG_diff_gff = pd.concat(add_unaligned_to_diff)



"""We now check to see which genes are conserved, degraded or unconserved based on overlap with gaps alone"""

#We remove DUPs from the pre interSEQ diff file to have a file which includes only gaps

#We store a version with DUP

mummer_diff_w_DUP = mummer_diff

mummer_diff = mummer_diff[(mummer_diff['TYPE'] != 'DUP')]

mummer_diff_test = mummer_diff[(mummer_diff['SEQ'] == 'neff_scaffold_166')]


#We go through each chromosome present in the file and concatenate overlapping gaps

merged_gaps = []

#print(dict_scaff['neff_scaffold_166'])

for scaffold in dict_scaff:

    scaff_list = []

    #We get all lines containing SEQ

    mummer_gap_chrom = mummer_diff[(mummer_diff['SEQ'] == scaffold)]

    list_of_coords = []

    for index, rows in mummer_gap_chrom.iterrows():
        coord_list = [rows.S1, rows.E1] 
        list_of_coords.append(coord_list)

    def recursive_merge(inter, start_index = 0):
        for i in range(start_index, len(inter) - 1):
            if inter[i][1] >= (inter[i+1][0] - overlap_value):
                new_start = inter[i][0]
                if inter[i+1][1] >= inter[i][1]:
                    new_end = inter[i+1][1]
                else:
                    new_end = inter[i][1]
                inter[i] = [new_start, new_end]
                del inter[i+1]
                return recursive_merge(inter.copy(), start_index=i)
        return inter    

    sorted_on_start = sorted(list_of_coords)
    merged = recursive_merge(sorted_on_start.copy())

    for item in merged: item.append(str(scaffold))

    merged_gaps.extend(merged)

merged_gaps_df = pd.DataFrame(merged_gaps, columns=['S1', 'E1', 'SEQ'])

#We add unaligned contigs to the dataframe

header2 = ['S1', 'E1', 'SEQ']

add_unaligned_df = pd.DataFrame(columns=header2)

add_unaligned_df['SEQ'] = unref['SEQ']
add_unaligned_df['S1'] = unref['S1']
add_unaligned_df['E1'] = unref['E1']

for_merge = [add_unaligned_df, merged_gaps_df]

all_gaps = pd.concat(for_merge)

#We go through viral genes and check whether they overlap with gaps:

for index, row in all_genes.iterrows():
    start = row['S1']
    end = row['E1']
    chrom = row['SEQ']

    #We grab any line where any item has any overlap with the query
    any_overlap = all_gaps[((((all_gaps['S1'] >= start) & (all_gaps['S1'] <= end)) | ((all_gaps['E1'] >= start) & (all_gaps['E1'] <= end))) | ((all_gaps['S1'] <= start) & (all_gaps['E1'] >= end))) & (all_gaps['SEQ'].isin([chrom]))]
    any_overlap = any_overlap.reset_index(drop=True)

    #We separate subject coordinates that fully encompass the query. These are the UNCONSERVED gens
    englobed_items = any_overlap[((any_overlap['S1'] <= start) & (any_overlap['E1'] >= end)) & (any_overlap['SEQ'].isin([chrom]))]
    #print(englobed_items)

    #Any coordinates that do not fully encompass are then extracted:
    non_englobed_items = any_overlap[~((any_overlap['S1'] <= start) & (any_overlap['E1'] >= end)) & (any_overlap['SEQ'].isin([chrom]))]

    status = ''

    if len(englobed_items) > 0:
        status='UNCONSERVED'
    elif len(non_englobed_items) > 0:
        status='DEGRADED'
    else:
        status='CONSERVED'
    
    all_genes.loc[index, 'ATTRIBUTE'] = all_genes.loc[index, 'ATTRIBUTE'].split(";")[0]+";"+str(status)



test = all_genes[all_genes['ATTRIBUTE'].str.contains('DEGRADED')]


#We make a bed file of diff and all differences

bed_differences_diff = pd.DataFrame(columns=['SEQ', 'S1', 'E1', 'TYPE'])

bed_differences_diff['S1'] = mummer_diff_w_DUP['S1'] - 1
bed_differences_diff['E1'] = mummer_diff_w_DUP['E1']
bed_differences_diff['SEQ'] = mummer_diff_w_DUP['SEQ']
bed_differences_diff['TYPE'] = mummer_diff_w_DUP['TYPE']

bed_differences_unaligned = pd.DataFrame(columns=['SEQ', 'S1', 'E1', 'TYPE'])

bed_headers = ['SEQ', 'S1', 'E1', 'TYPE']

bed_differences_unaligned['S1'] = add_unaligned_df['S1'] -1
bed_differences_unaligned['E1'] = add_unaligned_df['E1']
bed_differences_unaligned['SEQ'] = add_unaligned_df['SEQ']
bed_differences_unaligned['TYPE'] = 'NO_HIT'

bed_merge = [bed_differences_unaligned, bed_differences_diff]

bed_with_un = pd.concat(bed_merge)

#Concat gap output.gff (THE IMPORTANT ONE)
all_concat_gaps_gff.to_csv(str(output_directory)+str(output_prefix)+"_divergent_regions.gff", sep="\t", index=False, columns=header, header=None)

#Diff file WITHOUT the interseq but with unaligned
bed_with_un.to_csv(str(output_directory)+str(diff_input.split("/")[-1])+"_with_unaligned_scaffolds.bed", sep="\t", index=False, columns=bed_headers, header=None)

#Diff with interseq and unaligned scaffolds output.gff
CONTIG_diff_gff.to_csv(str(output_directory)+str(diff_input.split("/")[-1])+"_with_unaligned_scaffolds_and_translocations.gff", sep="\t", index=False, columns=header, header=None)

#All genes with conservation info output.gff
all_genes.to_csv(str(output_directory)+str(output_prefix.split("/")[-1])+"_genes_conservation.gff", sep="\t", index=False, columns=header, header=None)


