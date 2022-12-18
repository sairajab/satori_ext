from random import choice, randint
from readpwms import get_motif_proteins, get_motif_seq
import random
from mm_seq import GCRichSeq
import numpy as np

def embed_motifs(seq, motifs, pos = True):
    positions = [0,21,42,65,87,108, 129, 150,171]#np.arange(0,175).tolist()
    #position = choice(positions)
    # for j in range(no_of_motifs):
    #     positions.append()
    #
    # position =  randint(1, len(seq)-20)
    random.shuffle(positions)
    #random.shuffle(dist)
    i = 0
    k = 0
    last_tf = ""
    motifs_count = {}
    for e in motifs:
        if i == 9:
            i = 0
        for cons in motifs[e]:
            position = positions[i]
            bp = len(motifs[e][cons])
            # while not set(np.arange(position, position + bp + 1)).issubset(set(positions)):
            #     print("finding")
            #     position = choice(positions)

            if pos:
                print("pos is" + str(position))
                seq = seq[:position] + motifs[e][cons] + seq[position + bp:]
            else:
                seq = seq[:position] + motifs[e][cons] + seq[position + bp:]
                #print(position)

        i = i + 1
        if len(seq) > 200:
            raise Exception("Sorry, length above 200");
    return seq

def embed_motifs_dist(seq, motifs, pos = True):
    '''
    This function embeds given motifs-pair in the sequence with in <dist> : [10,8,9,11] and returns sequence
    '''
    positions = [0,50,100,150]
    dist = [10,8,9,11]
    # for j in range(no_of_motifs):
    #     positions.append()
    #
    # position =  randint(1, len(seq)-20)
    random.shuffle(positions)
    random.shuffle(dist)
    i = 0
    k = 0
    last_tf = ""
    motif_to_pos = {} #Which motifs are embedded at which positions
    for e in motifs:
        if i == 4:
            i = 0
        position = positions[i]
        for cons in motifs[e]:
            bp = len(motifs[e][cons])
            motif_to_pos[e][cons] = position 
            '''If the sequence is positive sequence, embed with distance'''
            if pos: 
                seq = seq[:position] + motifs[e][cons] + seq[position + bp:]
                position = position + bp + dist[i]

            else:
                '''If sequence is not positive embed motifs randomly using <positions> variable'''
                seq = seq[:positions[i]] + motifs[e][cons] + seq[positions[i] + bp:]
                print(positions[i])
        i = i + 1
        if len(seq) > 200:
            raise Exception("Sorry, length above 200");
    return seq, motif_to_pos

def motifs2embed(pos = True):
    '''
    This method returns the ids of motifs to be embedded in positive and negative sequences. 
    Number of pairs : 80
    Number of unique TFs : 93
    '''
    pairs = []
    pair_choices = [1,2] # add 1 pair or 2 pairs
    pair_prob = [0.5,0.5] # Probability of picking 1 pair or 2
    pair_copy = [0.5, 0.5] # Probability of copying pair 
    if pos: #Positive sequence
        #randomly pick number of pairs
        no_pairs = pair_choices[np.random.choice(2, p=pair_prob)]
        for i in range(no_pairs):
            #Randomly pick which pair id from 80
            value = randint(0, 79)
            #randomly pick if you wanna copy a pair
            copy = pair_choices[np.random.choice(2, p=pair_copy)]
            if copy == 2:
                    pairs.append(value)
            pairs.append(value)
    else: #Negative sequences
        #Randomly pick number of motifs to embed in negative sequences
        no_tfs = randint(1,4)
        for i in range(no_tfs):
            #randomly pick motif to embed from 93 unique motifs
            value = randint(0, 92)
            pairs.append(value)
    return pairs

def generateData(no_of_positive, no_of_negative, seq_len, pwms, outfile_name = "seqs93_tf_2_a"):
    '''
    '''
    label_output = open(outfile_name + ".txt", "w")
    fasta_output = open(outfile_name + ".fa", "w")
    annotations_output = open(outfile_name + "_info.txt", "w")
    head = 1
    for i in range(no_of_positive):
        # Generates GC Rich sequence of a given length 
        seq = GCRichSeq(seq_len)
        # Get ids of motif-pairs to embed in positive sequence
        ids = motifs2embed()
        #get motif sequence everytime using letter probability matrix
        tf_motif = get_motif_seq(pwms, idx=ids, pair=True)
        #Embed motifs in positive sequence 
        seq_embed, motif_to_pos = embed_motifs_dist(seq, tf_motif)
        print(len(seq_embed))

        #Get info about which motifs are embedded at which positions and write in file
        info_str = ""
        for num in motif_to_pos:
            for motif in motif_to_pos[num]:
                info_str = motif + "[" + str(motif_to_pos[num][motif]) + "]:" + info_str
        # Write .fa and .txt files
        line = ">Pos:1-" + str(head) + "(+)\n"
        fasta_output.write(line)
        fasta_output.write(seq_embed + "\n")
        label_output.writelines("Pos\t1\t" + str(head) + '\t' + '1' + '\n')
        annotations_output.writelines("Pos\t1\t" + str(head) + '\t' + '1:' + info_str[:-1] + '\n')

        head = head + 1

    for i in range(no_of_negative):
        #Generate GC Rich Sequence
        seq = GCRichSeq(seq_len)
        # Get ids of motifs for negative sequence
        idxs = motifs2embed(pos=False)
        #get motif sequence everytime using letter probability matrix
        tf_motif = get_motif_seq(pwms, idxs, pair=False)
        #Embed motifs in negative sequence
        seq_embed = embed_motifs_dist(seq, tf_motif, pos=False) #embed_motifs
        #Get info about which motifs are embedded at which positions and write in file
        info_str = ""
        for num in motif_to_pos:
            for motif in motif_to_pos[num]:
                info_str = motif + "[" + str(motif_to_pos[num][motif]) + "]:" + info_str
        # Write .fa and .txt files
        line = ">Neg:2-" + str(head) + "(+)\n"
        fasta_output.write(line)
        fasta_output.write(seq_embed + "\n")
        label_output.writelines("Neg\t2\t" + str(head) + '\t' + '0' + '\n')
        annotations_output.writelines("Neg\t2\t" + str(head) + '\t' + '0:' + info_str[:-1] + '\n')

        head = head + 1

    label_output.close()
    annotations_output.close()
    fasta_output.close()
    return True


if __name__ == "__main__":

    path_to_meme = 'C:/Users/Saira/Documents/Genomics/deepRAM-master/deepRAM-master/motifs_database/Jaspar.meme'

    num_of_seq = 30000
    seq_len = 200
    protiens, pwms = get_motif_proteins(path_to_meme)
    success = generateData(30000,30000, 200,pwms) 



