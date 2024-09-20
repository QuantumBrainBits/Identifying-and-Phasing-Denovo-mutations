# Function n : STR_variation_from_read()
# Considering only Repeat overlapped CIGAR variants.

def STR_variation_from_read(read, start, end, repeat_len):

        
    # Finding out the way to get the CIGAR string pos in read, reference & type of variation.
    # Then figuring out the way to include the way to fecth the info in our set of reads.

    repeat_start  = start
    repeat_end    = end


    cigar_s = read.cigarstring
    cigar_t = sum(read.cigartuples, ())

    num_bases = list(cigar_t)[1::2]
    alignment_type = list(cigar_t)[::2]

    # We update while itterating through the read positions to cover the entire read positions with any len and capture the variations which encounter within the repeat region.
    current_ref_pos = read.reference_start-1

    # saving Insertions and deletions counts which are part of repeat region in the list.
    change_in_len_insert = {'1':0}
    change_in_len_del = []

    for indx in range(len(alignment_type)): # this loop runs for every read, itterating through M's I's D's.

        # Here we get the actual reference pos by substracting the reference pos with hard or soft clip bases.
        if alignment_type[indx] == 1:
            # adding the insertion counts to the dictionary change_in_len
            if current_ref_pos in range(repeat_start, repeat_end +  1 + 1) :
                change_in_len_insert['1'] += num_bases[indx]


        elif alignment_type[indx] == 2:
            # adding the deletion counts to the dict.
            del_end_pos = (current_ref_pos + num_bases[indx])
            l = [(change_in_len_del.append(1)) for del_pos in range(current_ref_pos+1, del_end_pos+1) if del_pos in range(repeat_start, repeat_end+1)] 
            current_ref_pos += num_bases[indx]


        elif alignment_type[indx] not in [4,5]:
            current_ref_pos += num_bases[indx]


    # change in the len of read overlapping the repeat.
    Insertion_substract_Deletion = int(change_in_len_insert['1']) - sum(change_in_len_del)
    len_diff_in_read = repeat_len + Insertion_substract_Deletion


    return len_diff_in_read
