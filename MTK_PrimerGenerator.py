import numpy as np
import readline
from datetime import date

######################################################
# Basic tools for sequences and primer generation
######################################################

complement = {'A': 'T', 'T': 'A', 'C':'G', 'G':'C'}
translate_dict = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V', 'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V', 'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V', 'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V', 'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A', 'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A', 'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D', 'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G', 'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G', 'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', 'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
part_end_dict = {'1forward': 'GCATCGTCTCATCGGTCTCACCCT', '1reverse': 'ATGCCGTCTCAGGTCTCACGTT', '2forward': 'GCATCGTCTCATCGGTCTCAAACG', '2reverse': 'ATGCCGTCTCAGGTCTCACATA', '3forward': 'GCATCGTCTCATCGGTCTCATATG', '3reverse': 'ATGCCGTCTCAGGTCTCAGGATCC', '3aforward': 'GCATCGTCTCATCGGTCTCATATG', '3areverse': 'ATGCCGTCTCAGGTCTCAAGAACC', '3bforward': 'GCATCGTCTCATCGGTCTCATTCT', '3breverse': 'ATGCCGTCTCAGGTCTCAGGATCC', '4forward': 'GCATCGTCTCATCGGTCTCAATCCTAA', '4reverse': 'ATGCCGTCTCAGGTCTCACAGC', '4aforward': 'GCATCGTCTCATCGGTCTCAATCC', '4areverse': 'ATGCCGTCTCAGGTCTCAGCCATTA', '4aIIforward': 'TCGCGTCTCATCCA', '4aIIreverse': 'ATGCCGTCTCAGGTCTCAGCCATTA', '4bforward': 'GCATCGTCTCATCGGTCTCATGGC', '4breverse': 'ATGCCGTCTCAGGTCTCACAGC', '5forward': 'GCATCGTCTCATCGGTCTCAGCTG', '5reverse': 'ATGCCGTCTCAGGTCTCATGTA', '6forward': 'GCATCGTCTCATCGGTCTCATACA', '6reverse': 'ATGCCGTCTCAGGTCTCAACTC', '7forward': 'GCATCGTCTCATCGGTCTCAGAGT', '7reverse': 'ATGCCGTCTCAGGTCTCATCGG', '8forward': 'GCATCGTCTCATCGGTCTCACCGA', '8reverse': 'ATGCCGTCTCAGGTCTCAAGGG', '8aforward': 'GCATCGTCTCATCGGTCTCACCGA', '8areverse': 'ATGCCGTCTCAGGTCTCAATTG', '8bforward': 'GCATCGTCTCATCGGTCTCACAAT', '8breverse': 'ATGCCGTCTCAGGTCTCAAGGG', '3bIforward': 'GCATCGTCTCATCGGTCTCATTCT', '3bIreverse': 'GAACGTCTCATGCG', '3bIIforward': 'GAACGTCTCACGCA', '3bIIreverse': 'ATGCCGTCTCAGGTCTCAGGATCC'}


def reverse_complement(seq):
    n = len(seq)
    rc_seq = ''
    for i in np.arange(n)[::-1]:
        rc_seq += complement[seq[i]]
    return(rc_seq)

def translate(dna_seq):

    if len(dna_seq)%3 != 0:
        dna_seq = dna_seq[:-(len(dna_seq)%3)].upper()
    else:
        dna_seq = dna_seq.upper()

    codons = [dna_seq[i : i + 3] for i in range(0, len(dna_seq), 3)]
    aa_seq = ''
    for i in codons:
        aa_seq += translate_dict[i]

    return(aa_seq)

def calculate_melting_temp(dna_seq):
    # calculate melting temp of a given sequence using a simple formula
    A = dna_seq.count('A')
    T = dna_seq.count('T')
    G = dna_seq.count('G')
    C = dna_seq.count('C')

    if len(dna_seq) < 14:
        Tm = (A + T) * 2 + (G + C) * 4

    else:
        Tm = 64.9 +41*(G + C - 16.4)/(A + T + G + C)

    return(round(Tm, 2))

def calculate_optimal_primer_length(seq, starting_ix, direction):
    # find the length of a primer at which the temp first exceeds 57 degrees centigrade
    n = 10
    sequence = seq[starting_ix : starting_ix + n]
    melt_temp = 0

    while (melt_temp < 57.): # changed to 57.0 on 2/16/19
        if direction == 'forward':
            sequence = seq[starting_ix : starting_ix + n]
        elif direction == 'reverse':
            sequence = seq[starting_ix - n + 1 : starting_ix + 1]
        melt_temp = calculate_melting_temp(sequence)
        n += 1

    return(n)


######################################################
# BSMBI and BSAI restriction site detection
######################################################


def find_silent_mutations(seq, n):
    # takes dna sequence (in frame) and nucleotide position to mutate
    # returns current codon, position of first nucleotide in codon, and options for silent mutation

    current_codon = seq[3 * (n//3) : 3 * (n//3 + 1)]
    possible_codons = []

    # find amino acid associated with position n
    translation = translate(seq)
    aa_of_interest = translation[n//3]

    # collect translation data
    keys = np.array(list(translate_dict.keys()))
    vals = np.array(list(translate_dict.values()))

    # which codons could represent the amino acid at this position?
    ix = np.where(vals == aa_of_interest)[0]
    possible_codons = [keys[i] for i in ix if keys[i] != current_codon]




    return(current_codon, 3 * (n//3), possible_codons)



def find_restriction_sites(dna_seq):

    dna_seq = dna_seq.upper()

    bsmbi_recog_seq_f = 'CGTCTC'
    bsmbi_recog_seq_r = 'GAGACG'
    bsai_recog_seq_f = 'GGTCTC'
    bsai_recog_seq_r = 'GAGACC'

    potential_recognition_sites = np.array([dna_seq[i : i + 6] for i in range(0, len(dna_seq) - 6)])

    bsmbi_for_sites = np.where(potential_recognition_sites == bsmbi_recog_seq_f)[0]
    bsmbi_rev_sites = np.where(potential_recognition_sites == bsmbi_recog_seq_r)[0]
    bsai_for_sites = np.where(potential_recognition_sites == bsai_recog_seq_f)[0]
    bsai_rev_sites = np.where(potential_recognition_sites == bsai_recog_seq_r)[0]

    return(bsmbi_for_sites, bsmbi_rev_sites, bsai_for_sites, bsai_rev_sites)


def number_restriction_sites(dna_seq):
    a, b, c, d = find_restriction_sites(dna_seq)
    return(len(a), len(b), len(c), len(d))

def number_reactions_needed(dna_seq):
    RS = find_restriction_sites(dna_seq)
    total = 0
    for rs in RS:
        total += len(rs)

    return(total + 1)

def expected_product_sizes(dna_seq):

    RS = find_restriction_sites(dna_seq)

    all_rs = [0]
    for i in RS:
        for j in i:
            all_rs.append(j)
    all_rs.append(len(dna_seq))

    incomplete_prod_sizes = np.diff(np.sort(np.array(all_rs)))
    # need to add 24 nts at each edge
    # and 15 nts for each internal site
    # default to 30
    pcr_add = np.zeros_like(incomplete_prod_sizes) + 30

    # except for ends where we need and additional 9 nts
    pcr_add[0] += 9
    pcr_add[-1] += 9

    complete_prod_sizes = incomplete_prod_sizes + pcr_add

    return(complete_prod_sizes)

# write a function to find overhangs (top strand) in sequences containing a restriction site
def find_overhangs(dna_seq):

    overhangs = []
    bsmbi_for_sites, bsmbi_rev_sites, bsai_for_sites, bsai_rev_sites = find_restriction_sites(dna_seq)

    if len(bsmbi_for_sites) > 0:
        for i in bsmbi_for_sites:
            overhangs.append(dna_seq[i + 7 : i + 11])

    if len(bsmbi_rev_sites) > 0:
        for i in bsmbi_rev_sites:
            overhangs.append(dna_seq[i - 5 : i - 1])

    if len(bsai_for_sites) > 0:
        for i in bsai_for_sites:
            overhangs.append(dna_seq[i + 7 : i + 11])

    if len(bsai_rev_sites) > 0:
        for i in bsai_rev_sites:
            overhangs.append(dna_seq[i - 5 : i - 1])

    return(overhangs)

def find_bsmbi_overhangs(dna_seq):

    overhangs = []
    bsmbi_for_sites, bsmbi_rev_sites, bsai_for_sites, bsai_rev_sites = find_restriction_sites(dna_seq)

    if len(bsmbi_for_sites) > 0:
        for i in bsmbi_for_sites:
            overhangs.append(dna_seq[i + 7 : i + 11])

    if len(bsmbi_rev_sites) > 0:
        for i in bsmbi_rev_sites:
            overhangs.append(dna_seq[i - 5 : i - 1])

    return(overhangs)


######################################################
# Putting it all together
######################################################

def find_silent_mutations_in_RS(seq, ix_0):

    a1, b1, c1, d1 = number_restriction_sites(seq)
    num_sites_start = a1 + b1 + c1 + d1

    recog_site_ix = np.arange(ix_0, ix_0 + 6)

    proposed_silent_mutations = []

    for ix in recog_site_ix:

        previous_codon, position_start, sil_muts = find_silent_mutations(seq, ix)

        for j in sil_muts:
            # check to make sure we have removed the restrcition site
            candidate_seq = seq[ : position_start] + j + seq[position_start + 3:]
            a2, b2, c2, d2 = number_restriction_sites(candidate_seq)
            num_sites_end = a2 + b2 + c2 + d2
            removal_condition = num_sites_end < num_sites_start

            # check to make sure all requested changes are within the recognition site
            changed_nts = np.where(np.array([candidate_seq[i] == seq[i] for i in range(len(seq))]) == False)[0]

            within_bounds_condition = len(np.intersect1d(recog_site_ix, changed_nts)) == len(changed_nts)

            if (removal_condition & within_bounds_condition):
                proposed_silent_mutations.append(previous_codon + str(position_start) + j)

    return(np.unique(proposed_silent_mutations))


def point_mutation_generator(seq, n, new_aa):
    # takes dna sequence (in frame) and nucleotide position to mutate
    # returns current codon, position of first nucleotide in codon, and options for silent mutation

    current_codon = seq[3 * (n//3) : 3 * (n//3 + 1)]
    possible_codons = []

    # find amino acid associated with position n
    translation = translate(seq)
    aa_of_interest = translation[n//3]

    # collect translation data
    keys = np.array(list(translate_dict.keys()))
    vals = np.array(list(translate_dict.values()))

    # which codons could represent the amino acid at this position?
    ix = np.where(vals == new_aa)[0]
    possible_codons = [keys[i] for i in ix if keys[i] != current_codon]

    return(current_codon, 3 * (n//3), possible_codons)


def generate_GG_PMut_primers(seq, ix, mutate_to):

    ##################################################################
    # these primers can be used to mutate internal restriction sites.
    # they won't interfere with the edge overhangs, but we
    # still need to check for compatibility in overall reaction
    ##################################################################

    # BASIC STRUCTURE OF ONE OF THESE MUTATION PRIMERS:
    # SPACER, BSMBI_SITE, 6 NUCLEOTIDES OF OUR CHOICE, BINDING SEQUENCE
    # THE SIX NUCLEOTIDES CAN BE USED TO REMOVE RESTRICTION SITES
    # OR TO INTRODUCE MUTATIONS

    # find overlaps that are compatible with the GG rxn
    spacer = 'GAA'
    bsmbi_site = 'CGTCTC'

    # FIND ALL POSSIBLE PRIMERS THAT CAN BE USED TO MUTATE seq STARTING AT ix
    # AND SUBSTITUTING mutate_to IN ITS PLACE
    forward_primers = []
    reverse_primers = []

    target_seq = seq[ : ix] + mutate_to + seq[ix + len(mutate_to) : ]

    for shift in np.arange(6):
        left = ix - (6 - len(mutate_to))  + shift
        right = ix + (6 - len(mutate_to))  + shift
        six_nuc_seq = target_seq[ left : right ]

        n_R = calculate_optimal_primer_length(seq, right, 'forward')
        n_L = calculate_optimal_primer_length(seq, left, 'reverse')

        binding_seq_for = seq[right : right + n_R]
        binding_seq_rev = reverse_complement(seq[left - n_L : left ])

        fp = spacer + bsmbi_site + six_nuc_seq + binding_seq_for
        rp = spacer + bsmbi_site + reverse_complement(six_nuc_seq) + binding_seq_rev

        # check to verify that overhangs do not interfere with part assembly overhangs "TCGG" and "GACC"
        fp_overhang = find_overhangs(fp)
        rp_overhang = find_overhangs(rp)
        overhang_condition = (fp != "TCGG") & (fp != "GACC") & (rp != "TCGG") & (rp != "GACC")

        # and that we have just one bsmbi site in the primer
        just_one_site_condition = (number_restriction_sites(fp) == (1,0,0,0)) & (number_restriction_sites(rp) == (1,0,0,0))

        if overhang_condition & just_one_site_condition:
            forward_primers.append(fp)
            reverse_primers.append(rp)

    return(forward_primers, reverse_primers)


def generate_GG_edge_primers(seq, part_num):

    part_specific_f = part_end_dict[str(part_num) + 'forward']
    part_specific_r = part_end_dict[str(part_num) + 'reverse']

    n_R = calculate_optimal_primer_length(seq, 0, 'forward')
    n_L = calculate_optimal_primer_length(seq, len(seq), 'reverse')

    forward = part_specific_f + seq[ : n_R]
    reverse = part_specific_r + reverse_complement(seq[- n_L - 1 : ])

    return(forward, reverse)


def generate_GG_protocol(seq, part_num, verbose):

    seq = seq.upper()

    if part_num in ['3', '3a', '3b', '4a']:
        if len(seq) % 3 != 0:
            contin = input('Sequence appears to be out of frame. Continue?\n')
            if contin:
                pass
            else:
                return()

        aa_seq = translate(seq)
        if aa_seq[-1] == '*':
            print('Warning! Translated sequence ends with a stop codon (or is out of frame)\n')
            y_n = input('Should I remove it? (y/n)')

            if y_n in ['y', 'Y', 'yes', 'Yes', 'YES']:
                print('Stop codon removed!')
                seq = seq[:-3]
            else:
                print('Not removed...')
                pass

        if aa_seq[0] == 'M':
            print('Warning! Start codon is not needed at beginning of sequence\n')
            y_n = input('Should I remove it? (y/n)')

            if y_n in ['y', 'Y', 'yes', 'Yes', 'YES']:
                print('Start codon removed!')
                seq = seq[3:]
            else:
                print('Not removed...')
                pass

    #########################################################
    # preliminary summary of restriction sites to be removed
    #########################################################
    if verbose == True:
        print('\n=====================================================')
        print('Sequence summary:')
        print('=====================================================\n')

    sites = find_restriction_sites(seq)
    site_types = ['BsmBI (forward)', 'BsmBI (reverse)', 'BsaI (forward)', 'BsaI (reverse)']
    sites_flat = np.sort(np.concatenate(sites))

    for i in range(4):
        if len(sites[i]) > 0 :
            if verbose == True:
                print('Sequence contains ' + str(len(sites[i])) + ' ' + site_types[i] + ' sites beginning at:')
                print(sites[i])

    #########################################################
    # check for very early or late restriction sites
    #########################################################
    # if restriction sites occurs within 50 nts of beginning
    # or within 50 nts of the end of the target sequence,
    # generate a pair of oligos that can be annealed and
    # phosphorylated


    #########################################################
    # initialize dictionaries that contain primerset information
    #########################################################
    forward_primer_seq = {}
    reverse_primer_seq = {}
    overhang_seq = {}

    #########################################################
    # make edge primers for part assembly
    #########################################################

    f1, r1 = generate_GG_edge_primers(seq, part_num)

    primer_sets = []
    primer_sets.append(['FOR'])
    overhang_seq['FOR'] = find_bsmbi_overhangs(f1)[0]


    #########################################################
    # For each restriction site:
    #    find silent mutations
    #    for each possible silent mutation:
    #        find the 6 primer sets that can be used to generate the mutations
    #        and find overhangs
    # Now select best combination of overhangs for assembly
    #########################################################

    if verbose == True:
        print('\n=====================================================')
        print('Silent Mutation Details:')
        print('=====================================================\n')

    for site in sites_flat:

        sub_primer_set = []

        potential_mutations = find_silent_mutations_in_RS(seq, site)

        if len(potential_mutations) > 0:
            if verbose == True:
                print('allowable mutations for recognition sequence beginning at ' + str(site) + ':')
                for pm in potential_mutations:
                    print(pm[:3], pm[3:-3], pm[-3:])
                    #print(seq[pm[0]] + str(pm[0]) + pm[1])
                print('...Generating primers for ' + str(len(potential_mutations)) + ' potential silent mutation(s)\n')


            for pm in potential_mutations:

                mut_site, new_nt = int(pm[3:-3]), pm[-3:]

                potential_primers_f, potential_primers_r = generate_GG_PMut_primers(seq, mut_site, new_nt)
                potential_overhangs = np.array([find_bsmbi_overhangs(i)[0] for i in potential_primers_f])

                for primer_number in range(len(potential_primers_f)):

                    str_id = pm + '_' + str(primer_number + 1)

                    po = potential_overhangs[primer_number]
                    overhang_seq[str_id] = po

                    f_prim = potential_primers_f[primer_number]
                    forward_primer_seq[str_id] = f_prim

                    r_prim = potential_primers_r[primer_number]
                    reverse_primer_seq[str_id] = r_prim

                    sub_primer_set.append(str_id)

            primer_sets.append(sub_primer_set)

    primer_sets.append(['REV'])
    overhang_seq['REV'] = find_bsmbi_overhangs(r1)[0]


    #########################################################
    # Test 10000 random combinations for overhang compatibility
    # and choose best option
    #########################################################

    if verbose == True:
        print('\n=====================================================')
        print('Designing ' + str(number_reactions_needed(seq)) + ' PCR reaction(s)...')
        print('=====================================================\n')

    np.random.seed(0)
    still_looking = True
    count = 0

    while (still_looking == True) and (count <= 1000):
        if verbose:
            print('Testing iteration: ' + str(count + 1))

        rand_prim_set = [i[np.random.randint(len(i))] for i in primer_sets]

        rand_prim_set_OH_for = np.array([overhang_seq[i] for i in rand_prim_set])
        rand_prim_set_OH_rev = np.array([reverse_complement(overhang_seq[i]) for i in rand_prim_set])
        rand_prim_set_OH = np.concatenate([rand_prim_set_OH_for, rand_prim_set_OH_rev])

        N = len(rand_prim_set_OH)

        # check that no overhangs share three consecutive bases that are the same
        # this also covers the case of identitical primers
        fail_condition_1 = []
        for a in range(N):
            rp_a = rand_prim_set_OH[a]
            for b in range(a + 1, N):
                rp_b = rand_prim_set_OH[b]
                fail_condition_1.append((rp_a[:-1] in rp_b) | (rp_a[1:] in rp_b))


        # check that no overhangs differ in only one base pair (e.g. TAAG and TTAG)?
        fail_condition_2 = []
        for a in range(N):
            rp_a = rand_prim_set_OH[a]
            for b in range(a + 1, N):
                rp_b = rand_prim_set_OH[b]
                b = np.sum(np.array([nuc for nuc in rp_a]) == np.array([nuc for nuc in rp_b])) >= 3
                fail_condition_2.append(b)


        # check that no overhang has GC content of 0% or 100%.
        fail_condition_3 = []
        for rp in rand_prim_set_OH:
            gc_content = np.sum(np.array([(nuc in 'CGcg') for nuc in rp]))
            fail_condition_3.append(gc_content == 0)
        if verbose:
            print('Found ' + str((np.sum(fail_condition_1) + np.sum(fail_condition_2) + np.sum(fail_condition_3))) + ' exceptions')
        if (np.sum(fail_condition_1) + np.sum(fail_condition_2) + np.sum(fail_condition_3)) > 0:
            if verbose:
             print(' Failed, Trying Next Primer Set...')
            still_looking = True
            count += 1

        else:
            if verbose:
                print(' Found a set!')
            best_set = rand_prim_set
            still_looking = False


    decision_forward_primers = [f1]
    decision_reverse_primers = []
    if verbose:
        print('Results: ', best_set, '\n')
    for i in best_set[1:-1]:
        decision_forward_primers.append(forward_primer_seq[i])
        decision_reverse_primers.append(reverse_primer_seq[i])
    decision_reverse_primers.append(r1)

    if verbose:
        print('Overhangs: ', rand_prim_set_OH, '\n')


    prod_sizes = expected_product_sizes(seq)

    if verbose == True:
        for i in range(len(prod_sizes)):
            print('PCR Reaction ' + str(i + 1) + ', Expected Size: ' + str(prod_sizes[i]) + ' bp')
            print('Forward Primer:')
            print(decision_forward_primers[i])
            print('Reverse Primer:')
            print(decision_reverse_primers[i])
            print()
            
    # figure out the new sequence for optional generation of a genbank file
    output_seq = seq
    
    for i in best_set[1:-1]:
        pm = i.split('_')[0]
        print(pm[:3], pm[3:-3], pm[-3:])
        N = int(pm[3:-3])
        output_seq = output_seq[:N] + pm[-3:] + output_seq[N + 3:]
        
    output_seq = rand_prim_set_OH[0] + output_seq + rand_prim_set_OH[-1]
    
    return(decision_forward_primers, decision_reverse_primers, output_seq)



def generate_order_form(primers, prefix):
    # meant to take as input the output of the function above
    primers_f, primers_r = primers

    n = len(primers_f)

    for i in range(n):
        print(prefix + '_P' + str('{:02d}'.format(i+1)) + '_F\t' + primers_f[i][0:60] + '\t25nm\tSTD')
        print(prefix + '_P' + str('{:02d}'.format(i+1)) + '_R\t' + primers_r[i][0:60] + '\t25nm\tSTD')

    print('\n')

    for i in range(n):
        print(prefix + '_P' + str('{:02d}'.format(i+1)) + '_F, ' + primers_f[i][0:60])
        print(prefix + '_P' + str('{:02d}'.format(i+1)) + '_R, ' + primers_r[i][0:60])

        
def generate_gb_file(insert_seq, insert_name, plasmid_name):
    
    #######################################################################
    # Generate a genbank file with your domesticated insert of interest in
    # the MTK0_027 backbone
    #######################################################################
    

    save_to_directory = '' # <- REPLACE THIS WITH AN APPROPRIATE LOCATION
    
    
    pre_insert = 'tcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtgGGTCTCa'
    post_insert = 'tGAGACCagaccaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaattacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaa'

    full_sequence = pre_insert + insert_seq + post_insert

    bp_length = str(len(full_sequence))
    insert_length = len(insert_seq)
    today = date.today()
    date_string = today.strftime("%d-%b-%Y").upper()
    part_type = '3b'

    camRTerm_seq = 'accaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaa'
    camR_seq = 'ttacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccat'
    camRProm_seq = 'tttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaa'
    camRTerm_position = full_sequence.find(camRTerm_seq) + 1
    camR_position = full_sequence.find(camR_seq) + 1
    camRProm_position = full_sequence.find(camRProm_seq) + 1

    L0 = f'LOCUS       {plasmid_name}        {bp_length} bp ds-DNA     circular     {date_string}' #19-MAR-2021
    L1 = f'DEFINITION  Mammalian toolkit cloning part containing {insert_name}'
    L2 = 'ACCESSION   <unknown id>                                                       '
    L3 = 'VERSION     <unknown id>                                                       '

    L4 = 'FEATURES             Location/Qualifiers'
    L5 = '     rep_origin      complement(1..764)'
    L6 = '                     /label="ColE1"'
    L7 = '                     /ApEinfo_revcolor=#7f7f7f'
    L8 = '                     /ApEinfo_fwdcolor=#7f7f7f' 
    L9 = f'     misc_feature    777..{str(777 + int(insert_length - 8 - 1))}'
    L10 = f'                     /label="{insert_name}"'
    L11 = f'     misc_feature    complement({camRTerm_position}..{camRTerm_position + len(camRTerm_seq) - 1})'
    L12 = '                     /label="CamR Terminator"'
    L13 = '                     /ApEinfo_revcolor=#84b0dc'
    L14 = '                     /ApEinfo_fwdcolor=#84b0dc'       
    L15 = f'     CDS             complement({camR_position}..{camR_position + len(camR_seq) - 1})'
    L16 = '                     /label="CamR"'
    L17 = '                     /ApEinfo_revcolor=#0000ff'
    L18 = '                     /ApEinfo_fwdcolor=#0000ff'

    L19 = f'     promoter        complement({camRProm_position}..{camRProm_position + len(camRProm_seq) - 1})'
    L20 = '                     /label="CamR Promoter"'
    L21 = '                     /ApEinfo_revcolor=#84b0dc'
    L22 = '                     /ApEinfo_fwdcolor=#84b0dc'
    L23 = '        '
    L24 = 'ORIGIN'

    gb_line_compilation = [L0, L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12, L13, L14, L15, L16, L17, L18, L19, L20, L21, L22, L23, L24]

    seq_to_print = full_sequence.lower()
    groups_of_10 = [seq_to_print[(10*i):(10*(i+1))] for i in range(1 + len(seq_to_print)//10)]

    for i in range(len(groups_of_10)//6 + 1):
        nucs_so_far = str(60*i + 1)
        spacer = ''
        for s in range(5 - len(nucs_so_far)):
            spacer += ' '
        line_string = spacer + nucs_so_far + ' '
        for k in range(6):
            try:
                line_string += groups_of_10[6*i + k] + ' '
            except:
                pass
        gb_line_compilation.append(line_string)

    gb_line_compilation.append('//')
    
    file_path = save_to_directory + plasmid_name + '.gb'
    f = open(file_path, 'w')
    f.writelines(line + '\n' for line in gb_line_compilation)
    f.close()
    
    
# get a sorted list of all currently supported parts
def get_parts():
    parts = []
    for key in part_end_dict.keys():
        if "forward" in key:
            parts.append(key.replace("forward", ""))
    parts.sort()
    return parts

if __name__ == "__main__":
  seq = input("Enter desired nucleotide sequence (in frame if CDS)\n")
  seq_name = input('What is the name (annotation) of the nucleotide sequence\n')
  part_type = input("Enter desired part type ({})\n".format(", ".join(get_parts())))
  prefix = input("Enter a prefix for primer order form\n")
  print('Sequence is ' + str(len(seq)) + ' nucleotides')
  primers_f, primers_r, output_sequence = generate_GG_protocol(seq, part_type, True)
  primers = [primers_f, primers_r]
  generate_order_form(primers, prefix)
  
  plasmid_name = 'MTK' + part_type + '_' + seq_name
  generate_gb_file(output_sequence, seq_name, plasmid_name)
