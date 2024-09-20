#!/usr/bin/env python

# This is code takes input as your interested STR regions and SNP loci and check if any read in the a perticular STR region is enclosing the read and the same time the Snp loci (matching the snp allele with the reference or alternative allele) finally it gives the variation occured in STR region which is associating with the Snp allele.

# packages

import gzip
from tqdm import tqdm
from collections import Counter
import pysam
import sys



def Determine_PofO(NHsnps_regions, Aligned_file_loc):
    
    samfile = pysam.AlignmentFile(f'{Aligned_file_loc}', "rc")

    Determined_PofO = []
    
    for repeat_region in tqdm(NHsnps_regions):
        
        # info
        start =  int(repeat_region[1])  
        end =    int(repeat_region[2])
        chrom = repeat_region[0]
        rgn = f'{chrom}-{start}'

        
        snp_pos = int(repeat_region[-5])-1 # comparing with read.get_aligned_pairs(with_seq=True)
        snp_gts = [repeat_region[-4], repeat_region[-3]]
        snp_gts_representation = [repeat_region[-8], repeat_region[-7], repeat_region[-6]]
        denovo_gts = repeat_region[5].split('/')
        repeat_len = int(repeat_region[3])

        # here we are calling the function to store the snp allele parent of origin.
        snp_allele_origin = snp_allele_origin_classification(snp_gts, snp_gts_representation)

        # dictionary to maintian the snp_loc and repeat region enclosing the reads.
        read_info = {}

        
        # Here we will itterate the reads which enclose the repeat region and also NHsnp.
        for read in samfile.fetch(chrom, snp_pos, snp_pos):

            #
            if read.mapping_quality < 1: continue

            #
            read_query = read.query_name
            read_start = read.reference_start
            read_end = read.reference_end

            # checking if the ( read is enclosed) and ( covering the interested snp position in the flank). { R1:[s,e],[reference pos],[read base at reference pos] }
            if not read_query in read_info: 
                read_info[read_query] = [] # [STR_alleles, NHsnp]

                #
                base_info = sum(read.get_aligned_pairs(with_seq=True), ())
                snp_on_read_info = snp_check(read, snp_pos, snp_gts, base_info)           
                read_info[read_query].append(snp_on_read_info)


        # # Here we will itterate the reads which enclose the NHsnp.
        for read in samfile.fetch(chrom, start, end):

            # Here we check if the reads in the STR region, are enclosing the STR region and NHsnp pos.
            if read.mapping_quality < 1: continue
            if read.reference_start <= start and read.reference_end >= end :

                # Here we update the "read_info" dictionary with reads only enclosing the NHsnps
                str_read_query = read.query_name

                #
                if str_read_query in read_info:
                    #
                    Read_lenvar_info = STR_variation_from_read(read, start, end, repeat_len)
                    read_info[str_read_query].append(Read_lenvar_info)

        
        # # Calling the function which gives Snps associated with STR alleles.
        Snp_AssociatedSTR_GT_print = Snp_AssociatedSTR_GT(read_info, snp_gts,denovo_gts, start, end, snp_pos)
        

        #
        if Snp_AssociatedSTR_GT_print != None and len(Snp_AssociatedSTR_GT_print) > 1:
            Determined_PofO.append([*repeat_region[:5],repeat_region[5],repeat_region[6],repeat_region[7],repeat_region[8],f'{snp_gts[0]}|{snp_gts[1]}',f'{snp_allele_origin[0]}:{snp_allele_origin[1]};{snp_allele_origin[2]}:{snp_allele_origin[3]}', *Snp_AssociatedSTR_GT_print[:-1],repeat_region[-1], repeat_region[-2]])

    
    samfile.close()

    # Here we are returning list of determined regions.
    return Determined_PofO


