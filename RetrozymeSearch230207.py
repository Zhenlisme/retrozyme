import os, re, subprocess, argparse, time,sys, shutil, gc
from Bio import SeqIO
from collections import defaultdict
from multiprocessing.pool import ThreadPool

## Try to avoid pairing the constraint loops of minimum hh1 (CUGANGA and GAAANNUH). Revised at 30/08/2022
## Fix the problem of the name system of hhm motifs. For example, in the previous version, the motifs NC_030677.2-516402-516508-p should be in the negative strand.
## To fix the problem that the LTR type of retrozyme contains much satellite sequence.  Revised on 13/09/2022.
## To solve the problem that some retrozymes don't contain hammerhead at all.  1) Find active retrozyme firstly, and 2) find their homologous counterparts. 15/09/2022
## To detect the full copies of retrozyme in genomes. If it is a retrozyme, then its full copies should be detectable in different regions (distances greater than 1000 bp) of the genome. 29/09/2022
                ########## 1. To make consensus sequences of all repeated units or just use the first unit of tandem repeat as query
                ########## 2. To use the consensus sequences to blastn to find full copies.
## To split the blastn table into multiple pieces to accelerate the processing (as some genomes are too big to be handled well by the previous versions).   13/10/2022
class Hmm_Identify:
    def __init__(self, genome_file, hammerhead_description):
        self.genome_file = genome_file
        self.descr = hammerhead_description
    def hammerhead_search(self, optfile):
        with open(optfile, 'w') as rnabob_result:
            rnabob_program = subprocess.Popen(["rnabob", "-c", "-q", "-F", "-s", self.descr, self.genome_file],
                                              stderr=subprocess.DEVNULL, stdout=rnabob_result)
            rnabob_program.wait()
    def revise_rnabob_opt(self, rnabob_opt, topology_describe, ignore_list):
        """
        rnabob sometimes does not return as long as possible of helix. The loop length should be at least 4 nt long.
        """
        complement = {'A': "T", "T": "A", "G": "C", "C": "G",
                      "K": "M", "M": "K", "Y": "R", "R": "Y", "S": "S", "W": "W",
                      "B": "V", "V": "B", "H": "D", "D": "H", "N": "N", "X": "X"}
        rnabob_opt = rnabob_opt.upper().replace('U', 'T')

        rnabob_list = rnabob_opt.strip('|').split('|')
        topology_list = topology_describe.split(' ')
        topology_dict = {topology_list[i]: rnabob_list[i] for i in range(len(rnabob_list))}
        left_helix = set(re.findall('(h\d+)', topology_describe))
        for lh in left_helix:
            lh_index = topology_list.index(lh)  ## the index for left helix
            rh_index = topology_list.index(lh + "'")  ## the index for right helix

            left_drop, right_drop = [], []
            if lh_index + 2 == rh_index:  ### to find the single loop that between paired helices
                loop_index = lh_index + 1  ## the index for loop
                loop_seq = rnabob_list[loop_index]  ## the sequence for the loop
                if len(loop_seq) <= 4:
                    continue
                loop_name = topology_list[loop_index]
                drop_length = 0
                mid_point = int((len(loop_seq) - 4) / 2)
                for i in range(mid_point):
                    left_nuc = loop_seq[i]
                    right_nuc = loop_seq[-i - 1]
                    if complement[left_nuc] == right_nuc or {left_nuc, right_nuc} == {'G', 'T'}:
                        drop_length += 1
                        left_drop.append(left_nuc)
                        right_drop.append(right_nuc)
                    else:
                        break
                right_drop.reverse()
                topology_dict[loop_name] = topology_dict[loop_name][drop_length:len(loop_seq) - drop_length]
                topology_dict[lh] = topology_dict[lh] + ''.join(left_drop)
                topology_dict[lh + "'"] = ''.join(right_drop) + topology_dict[lh + "'"]
            else:
                ls_name = topology_list[lh_index + 1]
                rs_name = topology_list[rh_index - 1]

                ls_seq = rnabob_list[lh_index + 1]
                rs_seq = rnabob_list[rh_index - 1]

                short_length = len(ls_seq) if len(ls_seq) <= len(rs_seq) else len(rs_seq)
                left_seq_range = short_length - 4
                if short_length <= 4:
                    continue
                if ls_name in ignore_list:  #### specifially for hh1.minimal
                    ## To obtain the location of key motifs that should not be paired.
                    constraint1_loc = [(i.start(), i.end()) for i in re.finditer('CTGA.GA', ls_seq)]
                    constraint2_loc = [(i.start(), i.end()) for i in re.finditer('GAAA..{0,6}?.T[ATC]', rs_seq)]
                    left_pair_length = constraint1_loc[0][0]
                    right_pair_length = len(rs_seq) - constraint2_loc[0][1]
                    left_seq_range = left_pair_length if left_pair_length < right_pair_length else right_pair_length
                    ### only the 5' end of first constraint and 3 ' end of second constraint could be paired. So only round 1 pair was considered.
                    ############### Round 1 to pair 5' end
                drop_length = 0
                for i in range(left_seq_range):
                    left_nuc = ls_seq[i]
                    right_nuc = rs_seq[-i - 1]
                    if complement[left_nuc] == right_nuc or {left_nuc, right_nuc} == {'G', 'T'}:
                        drop_length += 1
                        left_drop.append(left_nuc)
                        right_drop.append(right_nuc)
                    else:
                        break
                ############### Round 2 to pair 3' end
                left_drop2, right_drop2 = [], []
                lh2_name = topology_list[lh_index + 2]  ## try to pair the single strand into another end of helix
                short_length2 = short_length - len(left_drop)
                drop_length2 = 0
                if short_length2 > 4 and ls_name not in ignore_list:  ### skip round 2 pairing if loop is constrained
                    for i in range(short_length2 - 4):
                        left_nuc = ls_seq[-i - 1]
                        right_nuc = rs_seq[i]
                        if complement[left_nuc] == right_nuc or {left_nuc, right_nuc} == {'G', 'T'}:
                            drop_length2 += 1
                            left_drop2.append(left_nuc)
                            right_drop2.append(right_nuc)
                        else:
                            break

                right_drop.reverse()
                left_drop2.reverse()
                topology_dict[ls_name] = ls_seq[drop_length:len(ls_seq) - drop_length2]
                topology_dict[rs_name] = rs_seq[drop_length2:len(rs_seq) - drop_length]
                topology_dict[lh] = topology_dict[lh] + ''.join(left_drop)
                topology_dict[lh + "'"] = ''.join(right_drop) + topology_dict[lh + "'"]
                topology_dict[lh2_name] = ''.join(left_drop2) + topology_dict[lh2_name]
                topology_dict[lh2_name + "'"] = topology_dict[lh2_name + "'"] + ''.join(right_drop2)


        full_topology = [topology_dict[i] for i in topology_list]
        if '' in full_topology:
            return rnabob_opt
        else:
            full_topology = ''.join(['|', '|'.join(full_topology), '|'])
            return full_topology
    def convert(self, rnabob_opt, topology_describe, ignore_list):
        rnabob_opt = rnabob_opt.strip('|').split('|')
        topology_list = topology_describe.strip().split(' ')

        full_loops = []
        for idx in range(len(topology_list) - 1):
            tpname = topology_list[idx]
            if tpname.startswith('s') and tpname not in ignore_list:
                # if topology_list[idx - 1] == topology_list[idx + 1].strip("'"):
                seq = rnabob_opt[idx]
                if len(seq) >= 5:
                    full_loops.append(''.join(['>', topology_list[idx], '\n', seq, '\n']))

        loop_name_dict = {}
        if full_loops:
            echo_seq_program = subprocess.Popen(["echo", '-e', "".join(full_loops)], stdout=subprocess.PIPE)
            rnafold_program = subprocess.Popen(["RNAfold", '--noPS'], stdin=echo_seq_program.stdout,
                                               stdout=subprocess.PIPE)
            rnafold_program.wait()
            for i in rnafold_program.stdout.readlines():
                line = repr(i).replace("b'", "").replace("\\n'", "")
                if line.startswith('>'):
                    loop_name = line.strip('>\n')
                elif re.findall('\d', line):
                    loop_name_dict[loop_name] = re.split('\s+', line)[0]


        opt_topology = []
        for idx in range(len(rnabob_opt)):
            topology_name = topology_list[idx]
            length_element = len(rnabob_opt[idx])
            if re.fullmatch('h\d+', topology_name):
                opt_topology.append((idx, '(' * length_element))
            elif topology_name.endswith("'"):
                opt_topology.append((idx, ')' * length_element))
            else:
                if topology_name in loop_name_dict:
                    opt_topology.append((idx, loop_name_dict[topology_name]))
                else:
                    opt_topology.append((idx, '.' * length_element))

        opt_topology = [i[1] for i in sorted(opt_topology, key=lambda x: x[0])]
        return ''.join(opt_topology)
    def rnabob_opt_to_fa(self, input_file, topology_describe, ignore_list, optfile):
        oplines = []
        seq_dict = {}
        with open(input_file, 'r') as F:
            for line in F:
                if re.match('\d', line.strip()):
                    location = tuple(re.split('\s+', line.strip())[:3])
                    seq_dict[location] = ''
                else:
                    seq_dict[location] = line.rstrip()
        for location in seq_dict:
            sequence = seq_dict[location]
            revised_seq = self.revise_rnabob_opt(sequence, topology_describe, ignore_list)
            dotbracket = self.convert(revised_seq, topology_describe, ignore_list)
            oplines.append(''.join(['>', '|'.join(location), '\n']))
            oplines.append(''.join([sequence.replace('|', ''), '\n']))
            oplines.append(''.join([dotbracket, '\n']))
        with open(optfile, 'w') as F:
            F.writelines(oplines)
    def calculate_gibbs_energy(self, input_file, opbed):
        echo_seq_program = subprocess.Popen(["cat", input_file], stdout=subprocess.PIPE)
        inter_rnaeval_file = input_file + '.rnaeval'
        with open(inter_rnaeval_file, 'w') as interfile:
            rnaeval_program = subprocess.Popen(["RNAeval"], stdin=echo_seq_program.stdout, stdout=interfile)
            rnaeval_program.wait()
        oplines = []
        with open(inter_rnaeval_file, 'r') as F:
            for line in F:
                if line.startswith('>'):
                    location = line.strip('>\n').split('|')
                    name = location[-1]
                    start, end, strand = [location[0], location[1], '+'] if int(location[0]) < int(location[1]) else [
                        location[1], location[0], '-']
                elif re.match('[a-zA-Z]', line):
                    sequence = line.rstrip()
                else:
                    splitlines = line.rstrip().split(' ')
                    dotbracket = splitlines[0]
                    energy_value = ''.join(splitlines[1:]).replace('(', '').replace(')', '').strip()
                    oplines.append([name, start, end, energy_value, sequence, strand, dotbracket])
        os.remove(inter_rnaeval_file)
        
        oplines=sorted(oplines, key=lambda x:[x[0], int(x[1])])
        with open(opbed, 'w') as F:
            for line in oplines:
                F.write('\t'.join(line))
                F.write('\n')
    def Search_with_Motif(self, opfile):
        with open(self.descr, 'r') as F:
            topology = F.readline().rstrip()
            ignor_name_list = [line.split(' ')[0] for line in F if 'CUGANGA' in line or 'GAA' in line]
        rnabob_opfile = ''.join([self.genome_file, '.rna_bob.txt'])
        self.hammerhead_search(rnabob_opfile)
        structure_file = ''.join([self.genome_file, '.structure.txt'])
        self.rnabob_opt_to_fa(rnabob_opfile, topology, ignor_name_list, structure_file)
        self.calculate_gibbs_energy(structure_file, opfile)
        os.remove(rnabob_opfile)
        os.remove(structure_file)

class Structure_analyze:
    def __init__(self, genome, genome_blastn_db, wkdir, process_num):
        self.process_num = int(process_num)
        self.wkdir = wkdir
        if not os.path.exists(self.wkdir):
            os.mkdir(self.wkdir)
            os.chdir(self.wkdir)
        else:
            os.chdir(self.wkdir)

        self.hh1descr = 'hh1.descr'
        self.hh1mdescr = 'hh1m.descr'

        with open(self.hh1descr, 'w') as F:
            F.write("""h1 s1 h2 s2 h3 s3 h3' s4 h4 s5 h4' s6 h2' s7 h1'\nh1 0:0 NNN[3]:[3]NNN\nh2 0:0 NNN[9]:[9]NNN\nh3 0:0 NNN[9]:[9]NNN\nh4 0:0 AYN[9]:[9]NRU\ns1 0 NNN[9]\ns2 0 CUGANGA[2]\ns3 0 NNN[12]\ns4 0 GAA\ns5 0 NNN[12]\ns6 0 H\ns7 0 NNN[9]\n""")
        with open(self.hh1mdescr, 'w') as F:
            F.write("""s1 h1 s2 h2 s3 h2' s4 h1' s5\nh1 0:0 NN[10]:[10]NN\nh2 0:0 NN[9]:[9]NN\ns1 0 NNN\ns2 0 CUGANGA[1]\ns3 0 NNN[15]\ns4 0 GAAAN[6]NUH\ns5 0 NNN\n""")

        if not os.path.exists('trf'):
            os.mkdir('trf')

        if not os.path.exists('HMM_cluster'):
            os.mkdir('HMM_cluster')

        if not os.path.exists('rtztbl'):
            os.mkdir('rtztbl')

        if os.path.exists('Clusters/'):
            shutil.rmtree('Clusters/')
            os.mkdir('Clusters/')
        else:
            os.mkdir('Clusters/')
            
        if not os.path.exists('genomes'):
            os.mkdir('genomes')

        self.genome=genome
        self.genomdb=genome_blastn_db
    def intersect(self, location1, location2, slip=0, lportion=0.0, rportion=0.0):
        location1 = sorted([int(location1[0]), int(location1[1])])
        location2 = sorted([int(location2[0]), int(location2[1])])
        if location1[0] - slip > location2[1] or location1[1] < location2[0] - slip:
            return False
        else:
            total_list = sorted([location1[0], location1[1], location2[0], location2[1]])
            portion1 = (total_list[2] - total_list[1] + 1) / (location1[1] - location1[0] + 1)  ##  proportion of the intersected sequence occupied on seq1
            portion2 = (total_list[2] - total_list[1] + 1) / (location2[1] - location2[0] + 1)  ##  proportion of the intersected sequence occupied on seq2
            if portion1 >= lportion or portion2 >= rportion:
                return True
            else:
                return False
    def run_trf(self, datfile):
        bnm = os.path.basename(datfile)
        trf_list = []
        with open(datfile, 'r') as F:
            header = F.readline()
            chrmid = header.split(' ')[0].strip('@\n')
            for line in F:
                splitlines = line.rstrip().split(' ')
                if float(splitlines[4]) < 50:  ## hammerhead motif usually longer than 50 nt
                    continue
                start, end = splitlines[0:2]
                repeat_time, concensus_size = splitlines[3:5]
                concensus_seq = splitlines[13]
                full_seq = splitlines[14]
                trf_list.append([chrmid, start, end, repeat_time, concensus_size, concensus_seq, full_seq])
        trf_list = sorted(trf_list, key=lambda x: [int(x[1]), int(x[2])])
        Final_lines = []
        if not trf_list:
            pass
        elif len(trf_list) == 2:
            first_line, second_line = trf_list
            if self.intersect(first_line[1:3], second_line[1:3], lportion=0.8, rportion=0.8):  ## The first and second one are the same. Select the one with the shorter consensus size.
                if float(first_line[4]) < float(second_line[4]):  ## The current line has a shorter consensus size.
                    Final_lines.append(first_line)
                else:
                    Final_lines.append(second_line)
            else:
                Final_lines.append(first_line)
                Final_lines.append(second_line)
        elif len(trf_list) == 1:
            Final_lines = trf_list
        else:
            FirstLine = trf_list[0]
            for trf_case in trf_list[1:]:
                if self.intersect(FirstLine[1:3], trf_case[1:3], lportion=0.8, rportion=0.8):  ## The first and second one are the same. Select the one with shorter consensus size.
                    if float(trf_case[4]) < float(FirstLine[4]):  ## The current line has a shorter consensus size.
                        FirstLine = trf_case
                else:  ## The first and second one are not the same, output the last one.
                    tobe_write = FirstLine
                    chrmid = tobe_write[0]
                    Final_lines.append(tobe_write)
                    FirstLine = trf_case
            Final_lines.append(FirstLine)
        trf_fa = ''.join(['trf/', bnm, '.trf.fa'])
        trf_location = {}
        trf_count = 1
        with open(trf_fa, 'w') as F:
            for trf_case in Final_lines:
                trf_name = '.'.join([trf_case[0], str(trf_count)])
                trf_count += 1
                trf_location[trf_name] = trf_case[:6]
                F.write(''.join(['>', trf_name, '\n']))
                #start, stop = int(trf_case[1]) - 1, int(trf_case[2])
                #F.write(str(self.genome_dict[trf_case[0]][start:stop]))
                F.write(trf_case[6])
                F.write('\n')
        return trf_location, trf_fa
    def Findclusters(self, bedfile, opbed, window=300):
        with open(opbed, 'w') as bedF:
            merge_bed_run = subprocess.Popen(['bedtools', 'cluster', '-i', bedfile, '-d', str(window), '-s'],
                                             stderr=subprocess.DEVNULL, stdout=bedF)
            merge_bed_run.wait()

        merge_dict = defaultdict(list)
        with open(opbed, 'r') as bedF:
            for line in bedF:
                line = line.strip().split('\t')
                cluster = line[-1]
                merge_dict[cluster].append(line[:-1])
        motif_count = 0
        active_count = 0
        rpt_times = []
        p_strand = 0
        n_strand = 0
        for cluster in merge_dict:
            coordlist = merge_dict[cluster]
            hhm_block_length=len(coordlist)
            strand = coordlist[0][5]
            if strand == '+':
                p_strand += hhm_block_length
            if strand == '-':
                n_strand += hhm_block_length
            rpt_times.append(hhm_block_length)
            motif_count+=hhm_block_length
            if hhm_block_length >=2:
                active_count+=1
        maximum_rpt=max(rpt_times)
        strand = '+' if p_strand >= n_strand else '-'
        return motif_count, active_count, maximum_rpt, strand
    def HHR_motif_decision(self, query_fa):
        hh1_bed = ''.join([query_fa, '.hh1'])
        hh1m_bed = ''.join([query_fa, '.hh1m'])
        Hmm_Identify(query_fa, self.hh1descr).Search_with_Motif(hh1_bed)
        Hmm_Identify(query_fa, self.hh1mdescr).Search_with_Motif(hh1m_bed)
        hh1_bed_dict = defaultdict(list)
        with open(hh1_bed, 'r') as F:
            for line in F:
                splitlines = line.rstrip().split('\t')
                hh1_bed_dict[splitlines[0]].append(line)

        hh1m_bed_dict = defaultdict(list)
        with open(hh1m_bed, 'r') as F:
            for line in F:
                splitlines = line.rstrip().split('\t')
                hh1m_bed_dict[splitlines[0]].append(line)

        motif_keys = list(hh1m_bed_dict.keys())
        motif_keys.extend(list(hh1_bed_dict.keys()))

        descr_lines = {}
        for name in set(motif_keys):
            hh1_list = hh1_bed_dict[name]
            hh1m_list = hh1m_bed_dict[name]
            if len(hh1_list) >= len(hh1m_list):
                hh1_list.append('hh1')
                descr_lines[name] = hh1_list
            else:
                hh1m_list.append('hh1m')
                descr_lines[name] = hh1m_list
        return descr_lines
    def trf_prog(self, hits_fa):
        bnm = os.path.basename(hits_fa)
        trf_opfile = ''.join(['trf/',bnm, '.dat'])
        with open(trf_opfile, 'w') as trF:
            TRF_task=subprocess.Popen(['trf', hits_fa ,'2', '7', '7', '80', '10', '50', '1000', '-h', '-ngs'], stderr=subprocess.DEVNULL, stdout=trF)
            TRF_task.wait()
        return trf_opfile
    def Main(self, datfile):
        bnm=os.path.basename(datfile).replace('.fa.dat', '')
        trf_location, trf_fa = self.run_trf(datfile)
        if not trf_location:
            with open('log.txt', 'a') as logF:
                logF.write(bnm)
                logF.write("\n")
            return 0

        descr_type_dict = self.HHR_motif_decision(trf_fa)

        TRF_num = 1
        retrozyme_candidates=[]
        for trfname in descr_type_dict:  ## To check trf region first
            chrmid, start, end, repeat_time, concensus_size, concensus_seq = trf_location[trfname]
            retrozyme_name = ''.join(['Rtz_', chrmid, '.', str(TRF_num)])
            TRF_num += 1
            hmm_lines = descr_type_dict[trfname]
            hhrtype = hmm_lines[-1]
            hmm_bed = ''.join(['HMM_cluster/', retrozyme_name, '.hmm.bed'])
            with open(hmm_bed, 'w') as F:
                F.writelines(hmm_lines[:-1])
            hmm_cluster = ''.join(['HMM_cluster/' , retrozyme_name, '.hmm.clust'])
            motif_count, active_count, rpt_times, strand = self.Findclusters(hmm_bed, hmm_cluster, window=int(concensus_size)+20)
            os.remove(hmm_bed)
            if int(active_count) >=1: #Only report cases with active retrozymes.
                retrozyme_candidates.append([chrmid, start, end, repeat_time, concensus_size, strand,'trf', retrozyme_name,
                                                str(motif_count), str(active_count), str(rpt_times), concensus_seq, hhrtype])
        retrozyme_table = ''.join(['rtztbl/', bnm, '.rtztbl'])
        if retrozyme_candidates:
            with open(retrozyme_table, 'w') as F:
                for line in retrozyme_candidates:
                    F.write('\t'.join(line))
                    F.write('\n')
        with open('log.txt', 'a') as logF:
            logF.write(bnm)
            logF.write("\n")
    def vsearch_clust(self, input_fa, opfile, marker, id=0.8):
        cons_name = ''.join([input_fa, '.reduce.temp'])
        cluster_file = ''.join([input_fa, '.clust.info'])
        ### vsearch can automaticaly output consencus sequence of clusters.
        run_cluster = subprocess.Popen(
            ['vsearch', '--cluster_fast', input_fa, '--id', str(id), '--strand', 'both', '--threads', str(self.process_num),
             '--clusterout_id', '--consout', cons_name, '--blast6out', cluster_file, '--quiet'],
            stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        run_cluster.wait()

        retrozyme_seq = []
        with open(cons_name, 'r') as F:
            for line in F:
                if line.startswith(">"):
                    series_num = line.split(';')[0].replace('>centroid=', '')
                    retrozyme_seq.append(series_num)

        cluster_info = defaultdict(list)
        with open(cluster_file, 'r') as F:
            for line in F:
                splitlines = line.rstrip().split('\t')
                cluster_name = splitlines[1]
                cluster_info[cluster_name].append(splitlines[0])

        cluster_dict = {}
        initio_num = 1
        for k in cluster_info:
            clust_list = cluster_info[k]
            clustname = '_'.join([marker, str(initio_num)])
            cluster_dict[k]=clustname
            for i in clust_list:
                cluster_dict[i] = clustname
            initio_num += 1
        monomoer_cluster_list = set(retrozyme_seq) - set(cluster_dict.keys())
        for k in monomoer_cluster_list:
            clustname = '_'.join([marker, str(initio_num)])
            cluster_dict[k] = clustname
            initio_num += 1
        opseq = ''
        with open(cons_name, 'r') as F:
            for line in F:
                if line.startswith(">"):
                    series_num = line.split(';')[0].replace('>centroid=', '')
                    cluster_name = cluster_dict[series_num]
                    opseq+=''.join(['>', cluster_name, '\n'])
                else:
                    opseq += line
        with open(opfile, 'w') as F:
            F.write(opseq)
        os.remove(cons_name)
        return cluster_dict
    def Copy_count(self, query_file, insertion_length_vary, distance):
        blastn_opt = ''.join([query_file, '.tbl'])
        blastn_pro = subprocess.Popen(
            ['blastn', '-db', self.genomdb, '-query', query_file, '-num_threads', str(self.process_num), '-max_target_seqs', '200000', '-evalue', '1e-5',
             '-outfmt', '6 qseqid sseqid pident qstart qend sstart send evalue qlen', '-out', blastn_opt])
        blastn_pro.wait()

        opdict = defaultdict(list)
        with open(blastn_opt, 'r') as F:
            for line in F:
                splitlines = line.rstrip().split('\t')
                query_name, chrm, identity, qstart, qend, sstart, send, evalue, qlen = splitlines[:9]
                hit_len = abs(int(qend) - int(qstart)) + 1
                min_len = min(insertion_length_vary[query_name])
                min_len = min([min_len, int(qlen)])
                coverage = hit_len/int(qlen)
                if float(identity) < 80 or hit_len < min_len * 0.8:
                    continue
                if int(sstart) < int(send):
                    START = sstart
                    END = send
                    strand = '+'
                else:
                    START = send
                    END = sstart
                    strand = '-'
                opdict[query_name].append([chrm, START, END, identity, ':'.join([str(coverage), str(min_len)]), strand])

        os.remove(blastn_opt)

        clust_count = {}
        for query_name in opdict:
            oplines = sorted(opdict[query_name], key=lambda x: [x[0], int(x[1])])
            blastn_bed = ''.join(['Clusters/', query_name, '.bed6'])
            with open(blastn_bed, 'w') as F:
                for line in oplines:
                    F.write('\t'.join(line))
                    F.write('\n')

            #### To cluster different full hits by distance
            clustblastn_bed = ''.join(['Clusters/', query_name, '.clust.bed6'])
            with open(clustblastn_bed, 'w') as F:
                cluster_bedtools = subprocess.Popen(['bedtools', 'cluster', '-i', blastn_bed, '-s', '-d', str(distance)],
                                                    stderr=subprocess.DEVNULL, stdout=F)
                cluster_bedtools.wait()

            os.remove(blastn_bed)
            with open(clustblastn_bed, 'r') as F:
                num_clusters = len({line.rstrip().split('\t')[-1] for line in F})
            clust_count[query_name]=str(num_clusters)
        return clust_count

    def MultipleMain(self):
        ## Run trf program firstly.
        if Args.trf:
            self.genome_dict = SeqIO.parse(self.genome, 'fasta')
            self.genome_dict = {k.id: k.seq.upper() for k in self.genome_dict}
            Total_genome_list = set(self.genome_dict.keys())
            trf_exist = set({i.replace('.fa.dat', '') for i in os.listdir('trf') if i.endswith('.fa.dat') and os.path.getsize(''.join(['trf/', i]))})
            if trf_exist:
                trf_left = Total_genome_list - trf_exist
            else:
                trf_left = Total_genome_list

            subgenome_list = []
            for chrm in trf_left:
                subgenome = ''.join(['genomes/', chrm, '.fa'])
                subgenome_list.append(subgenome)
                if os.path.exists(subgenome):
                    continue
                with open(subgenome, 'w') as F:
                    F.write(''.join(['>', chrm, '\n']))
                    F.write(str(self.genome_dict[chrm]))
                    F.write('\n')

            del self.genome_dict
            gc.collect()

            if subgenome_list:
                sys.stdout.write('%s candidates will be run by trf program!\n' % str(len(subgenome_list)))
                if len(subgenome_list) < self.process_num:
                    processnum = len(subgenome_list)
                else:
                    processnum = self.process_num
                planpool = ThreadPool(processnum)
                for subgenome in subgenome_list:
                    planpool.apply_async(self.trf_prog, args=(subgenome,))
                planpool.close()
                planpool.join()

        ## Run hammerhead detection program
        trf_datfile = [''.join(['trf/', i]) for i in os.listdir('trf') if i.endswith('.fa.dat') and os.path.getsize(''.join(['trf/', i]))]
        already_list = []
        if os.path.exists('log.txt'):
            with open('log.txt', 'r') as F:
                already_list = [''.join(['trf/', line.rstrip(), '.fa.dat']) for line in F]

        Tobe_run_list = set(trf_datfile) - set(already_list)

        if Tobe_run_list:
            if len(Tobe_run_list) < self.process_num:
                processnum = len(Tobe_run_list)
            else:
                processnum = self.process_num
            sys.stdout.write('%s candidates will be run by hammerhead search program!\n' % str(len(Tobe_run_list)))
            hammerheadpool = ThreadPool(processnum)
            for trfdat in Tobe_run_list:
                hammerheadpool.apply_async(self.Main, args=(trfdat,))
            hammerheadpool.close()
            hammerheadpool.join()

        retrozyme_list = [''.join(['rtztbl/', file]) for file in os.listdir('rtztbl')]
        if not retrozyme_list:
            sys.stdout.write('No retrozymes detected for %s!\n' % Args.opdir)
            return 0

        # To build consensus sequence.
        retrozyme_tbl = 'retrozyme.tbl'
        retrozyme_monomer = 'retrozyme.trf.monomer.fa'
        retrozyme_monomer_cons = 'retrozyme.trf.monomer.cons.fa'

        if not self.genomdb:
            self.genomdb = 'GenomeDB'
            mkblastdb = subprocess.Popen(['makeblastdb', '-dbtype', 'nucl', '-in', self.genome, '-out', self.genomdb],
                                         stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
            mkblastdb.wait()

        complement = {'A': "T", "T": "A", "G": "C", "C": "G",
                      "K": "M", "M": "K", "Y": "R", "R": "Y", "S": "S", "W": "W",
                      "B": "V", "V": "B", "H": "D", "D": "H", "N": "N", "X": "X"}

        insertion_length = {}
        with open(retrozyme_monomer, 'w') as wF:
            for rtztbl in retrozyme_list:
                with open(rtztbl, 'r') as rF:
                    for line in rF:
                        line = line.rstrip().split('\t')
                        if line[5] == '-': ## In case the consensus sequence is on reverse strand.
                            seq = list(line[11])
                            seq.reverse()
                            seq = ''.join([complement[i] for i in seq])
                        else:
                            seq = line[11]
                        insertion_length[line[7]] = line[4]
                        wF.write(''.join(['>', line[7],'\n']))
                        wF.write(''.join([seq, '\n']))
        ## Build consensus sequence
        subclass_dict = self.vsearch_clust(retrozyme_monomer, retrozyme_monomer_cons, 'subclass', id=0.8)
        ## To obtain the varied length of every class of retrozymes.
        insertion_length_vary=defaultdict(list)
        for i in subclass_dict:
            key = subclass_dict[i]
            length = insertion_length[i]
            insertion_length_vary[key].append(int(length))

        cluster_count_dict = self.Copy_count(retrozyme_monomer_cons, insertion_length_vary, 1000)

        ## To merge all information into a summary table
        with open(retrozyme_tbl, 'w') as wF:
            for rtztbl in retrozyme_list:
                with open(rtztbl, 'r') as rF:
                    for line in rF:
                        line = line.rstrip().split('\t')
                        if line[5] == '-': ## In case the consensus sequence is on reverse strand.
                            seq = list(line[11])
                            seq.reverse()
                            seq = ''.join([complement[i] for i in seq])
                            line[11]=seq
                        insertion_name = line[7]
                        subclass_name = subclass_dict[insertion_name]
                        if subclass_name in cluster_count_dict:
                            count = cluster_count_dict[subclass_name]
                        else:   ## N might exists in retrozymes
                            count = '0'
                        line.extend([subclass_name, count])
                        wF.write('\t'.join(line))
                        wF.write('\n')
                                    
if __name__ == "__main__":
    begin_time = time.time()
    parser = argparse.ArgumentParser(description="Hammerhead and retrozyme detection.")
    parser.add_argument("-g", "--genome", type=str, required=True, help="The reference genome.")
    parser.add_argument('-gdb', '--blastndb', type=str, required=False, default='',
                                  help='Blastn database name')
    parser.add_argument("-p", "--process", type=int, default=2, required=False,
                             help="Number of processes to be used.")
    parser.add_argument("-t", "--trf", type=int, default=0, choices = [0, 1], required=False, help="To run trf or not. default not run.")
    parser.add_argument("-o", "--opdir", type=str, required=True, help="The output file.")
    Args = parser.parse_args()

    if Args.blastndb:
        BLASTNDB=os.path.abspath(Args.blastndb)
        if not os.path.exists(''.join([BLASTNDB,'.ndb'])):
            BLASTNDB=''
    else:
        BLASTNDB=''

    Retrozyme_run=Structure_analyze(genome=os.path.abspath(Args.genome), genome_blastn_db=BLASTNDB,
                      wkdir=os.path.abspath(Args.opdir), process_num=Args.process)
    Retrozyme_run.MultipleMain()

    over_time = time.time()
    total_time = round((over_time - begin_time) / 3600, 4)
    print(total_time, ' hours used for retrozyme detection for %s.\n' % Args.opdir)
                 
