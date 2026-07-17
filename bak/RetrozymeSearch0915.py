import os, re, subprocess, argparse, time,sys
from Bio import SeqIO
from collections import defaultdict
from bs4 import BeautifulSoup
from multiprocessing.pool import ThreadPool

## Try to avoid pairing the constraint loops pf minimum hh1 (CUGANGA and GAAANNUH). revised at 30/08/2022
## Fix the problem of name system of hhm motifs. for example, in previous version, the motifs NC_030677.2-516402-516508-p should be in negative strand.
## To fix the problem that ltr type of retrozyme contain much of satellite sequence.  revised on 13/09/2022.
## To solve the problem that some retrozyme doesn't contain hammerhead at all.  1) Find active retrozyme firstly, and then find their homologous counterparts.

class Structure_analyze:
    def __init__(self, genome, hhmtbl):
        self.genome_dict = SeqIO.parse(genome, 'fasta')
        self.genome_dict = {k.id: k.seq.upper() for k in self.genome_dict}
        self.hhm_tbl=hhmtbl
        self.trfdistance=150
        self.ltrdistance=100
    def Iter_hhmtbl(self):
        with open(self.hhm_tbl,'r') as F:
            for line in F:
                splitlines=line.rstrip().split('\t')
                start=int(splitlines[1])-1000
                end=int(splitlines[2])+1000
                start=start if start>0 else 1
                strand=splitlines[-1]
                hhm_coord=[i.split('-') for i in splitlines[5].split(',')]
                yield (splitlines[0], start, end, strand, hhm_coord)

    def run_ltrharvest(self, hits_fa, dbname):
        db_tsk = subprocess.Popen(['gt', 'suffixerator', '-db', hits_fa, '-indexname', dbname, '-tis', '-suf', '-lcp', '-des', '-ssp', '-sds', '-dna'],
                                  stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        db_tsk.wait()
        ltr_task = subprocess.Popen(['gt', 'ltrharvest', '-index', dbname, '-mindistltr', '200'], stderr=subprocess.DEVNULL, stdout=subprocess.PIPE)
        ltr_task.wait()
        ltr_locations=[]
        for i in ltr_task.stdout.readlines():
            line=repr(i).replace("b'","").replace("\\n'","")
            if not line.startswith('#'):
                sub_locations=re.split('\s+',line.rstrip())
                lltr=[int(sub_locations[3]), int(sub_locations[4])]
                rltr=[int(sub_locations[6]), int(sub_locations[7])]
                ltr_locations.append([lltr,rltr])
        return ltr_locations

    def html_to_csv(self, htmlfile):
        path = htmlfile
        data = []
        # for getting the header from
        # the HTML file
        list_header = []
        soup = BeautifulSoup(open(path), 'html.parser')
        header = soup.find_all("table")[0].find("tr")

        for items in header:
            try:
                list_header.append(items.get_text())
            except:
                continue

        # for getting the data
        HTML_data = soup.find_all("table")[0].find_all("tr")[1:]

        for element in HTML_data:
            sub_data = []
            for sub_element in element:
                try:
                    sub_data.append(sub_element.get_text())
                except:
                    continue
            try:
                start, end = sub_data[0].split('--')
                data.append([start, end, sub_data[1], sub_data[2]])
            except:
                continue
        return data

    def run_trf(self, hits_fa, identifier):

        TRF_task=subprocess.Popen(['trf', hits_fa ,'2', '7', '7', '80', '10', '50', '500'],
                                  stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        TRF_task.wait()
        opt_list=[]
        for htmlfile in os.listdir():
            if htmlfile.startswith(identifier) and htmlfile.endswith('html') and not htmlfile.endswith('txt.html'):
                opt_list.extend(self.html_to_csv(htmlfile))
        opt_list=sorted(opt_list, key=lambda x:[x[0],int(x[1])])
        return opt_list

    def intersect(self, location1, location2, slip=0, lportion=0.0, rportion=0.0):
        location1 = sorted(location1)
        location2 = sorted(location2)
        if location1[0] - slip > location2[1] or location1[1] < location2[0] - slip:
            return False
        else:
            total_list = sorted([location1[0], location1[1], location2[0], location2[1]])
            portion1 = (total_list[2] - total_list[1] + 1) / (location1[1] - location1[0] + 1)  ## how much proportion the intersected sequence occupied on seq1
            portion2 = (total_list[2] - total_list[1] + 1) / (location2[1] - location2[0] + 1)  ## how much proportion the intersected sequence occupied on seq2
            if portion1 >= lportion and portion2 >= rportion:
                return True
            else:
                return False

    def MainRun(self, location_info):
        chrname, start, end, strand, hhm_coord=location_info
        HHM_Start, HHM_End = int(hhm_coord[0][0]), int(hhm_coord[-1][1])
        sequence=self.genome_dict[chrname][start-1:end]
        label='-'.join([chrname,str(start),str(end), strand])  ## To distinguish different files produced via multiple processing algorithm.
        sys.stdout.write('Processing block %s ...\n' % label)
        seq_file=''.join([label,'.fa'])
        with open(seq_file, 'w') as F:
            F.write(''.join(['>',label,'\n']))
            F.write(str(sequence))
        dbname=''.join([label,'.dbname'])

        retrozyme_list=[]

        ## run trf to detect tandem repeats
        trf_locations=self.run_trf(seq_file, label)
        ## make the trf as priority as ltr might be inside trf (if the ltr doesn't contain indertion refions)

        trf_status=0 ## set initial status as 0, meaning that there is no trf found in this region

        TRF_boundary=[]
        for trf in trf_locations:   ## A sequence might contain more than 1 trf blocks which can be overlapped.
            tstart,tend,pdsize,cpnum=trf
            tstart, tend, pdsize=start+int(tstart), start+int(tend), int(pdsize) ## adjust the tandem repeat coord to genome coord.
            TRF_boundary.extend([tstart, tend])  ## append trf unit location
            if pdsize<self.trfdistance:
                continue
            trf_loc_list=list(range(tstart, tend+1, pdsize))
            trf_loc_list.append(tend)   ## Remember to add the ending position
            trf_loc_list=[(trf_loc_list[i], trf_loc_list[i+1]) for i in range(len(trf_loc_list)-1)]  ## list every single unit in tandem repeat
            trf_num=0   ## the number of trf units that contain hammerhead motifs
            trf_label=[]

            count=-1 ### initialize as -1, then the first index will be 0
            for trf_loc in trf_loc_list:
                count+=1
                trf_label.append('n')
                for hhm_loc in hhm_coord:
                    if self.intersect(trf_loc, (int(hhm_loc[0]), int(hhm_loc[1])), rportion=0.95):  ## the hamerhead sequence should be inside trf unit
                        trf_num+=1
                        trf_label[-1]=str(count)
                        continue

            trf_label='-'.join(trf_label)
            trf_label=re.split('n', trf_label)
            for blc in trf_label:    ## To find the tandem repeat blocks: ---x---x---
                blc=blc.strip('-').split('-')
                if len(blc)>=2:
                    trf_start=trf_loc_list[int(blc[0])][0]
                    trf_end=trf_loc_list[int(blc[-1])][-1]
                    retrozyme_name='-'.join([chrname,str(trf_start), str(trf_end),'p']) if strand == '+' else '-'.join([chrname, str(trf_start), str(trf_end),'n'])
                    retrozyme_list.append([chrname, str(trf_start), str(trf_end), retrozyme_name, 'trf', str(pdsize), str(len(blc)), strand])

        if TRF_boundary:
            TRF_boundary=sorted(TRF_boundary)
            TRF_Start, TRF_End=TRF_boundary[0], TRF_boundary[-1]         ## To find the boundary of tandem region.
            if self.intersect([TRF_Start, TRF_End], [HHM_Start, HHM_End], rportion= 0.6):     ## if the whole trf region intersect with hammerhead region, that means this region is tandem repeats, then ltrharvest program will not run.
                trf_status = 1   ## set as 1 meaning that there are tandem repeats found

        if not trf_status:
            ####trfretrozyme_list=[sorted(trfretrozyme_list, key=lambda x:float(x[4]))[-1]]  ## To select the trf whose unit number is the largest
            ## run ltrharvest to detect if tandem repeats exist
            ltr_locations=self.run_ltrharvest(seq_file, dbname)
            for loc in ltr_locations:
                left_ltr, right_ltr=loc
                left_ltr=(left_ltr[0]+start, left_ltr[1]+start)  ## adjust the ltr coord to genome coord.
                right_ltr=(right_ltr[0]+start, right_ltr[1]+start)
                left_hhm,right_hhm=[],[]
                for hhm_loc in hhm_coord:
                    if self.intersect(left_ltr, (int(hhm_loc[0]), int(hhm_loc[1])), rportion=1):  ## the hamerhead sequence should be inside ltr repeats
                        left_hhm.append('-'.join(hhm_loc))
                    elif self.intersect(right_ltr, (int(hhm_loc[0]), int(hhm_loc[1])), rportion=1):
                        right_hhm.append('-'.join(hhm_loc))
                if left_hhm and right_hhm:
                    ## chrom id, LTR start, LTR end, distance between two LTR arms, number of hhm motifs in left LTR and right LTR, the location of hammerhead motifs (in left and right arms)
                    inter_seq_len=right_ltr[0]-left_ltr[1]
                    if inter_seq_len >= self.ltrdistance:
                        retrozyme_name = '-'.join([chrname, str(left_ltr[0]), str(right_ltr[1]), 'p']) if strand == '+' else \
                            '-'.join([chrname, str(left_ltr[0]), str(right_ltr[1]), 'n'])
                        retrozyme_list.append([chrname, str(left_ltr[0]), str(right_ltr[1]), retrozyme_name, 'ltr',
                                                  str(inter_seq_len), str(len(left_hhm)+len(right_hhm)), strand])

        os.system('rm -f %s*' % label)
        return retrozyme_list

    def getseq(self, location, strand):
        complement = {'A': "T", "T": "A", "G": "C", "C": "G",
                      "K": "M", "M": "K", "Y": "R", "R": "Y", "S": "S", "W": "W",
                      "B": "V", "V": "B", "H": "D", "D": "H", "N": "N", "X": "X"}
        start, end = location[1:3]
        ncnumber = location[0]
        seq = self.genome_dict[ncnumber][int(start)-1:int(end)]
        seq=str(seq).upper()
        if strand == '-':
            seq = list(seq)
            seq.reverse()
            seq = ''.join([complement[i] for i in seq])
        return seq

    def run_process(self):
        planpool = ThreadPool(int(Args.process))
        run_result=[]
        for location in self.Iter_hhmtbl():
            run_result.append(planpool.apply_async(self.MainRun, args=(location,)))
        planpool.close()
        planpool.join()
        opt=[]
        trf_seq_list=[]
        ltr_seq_list=[]
        for result in run_result:
            result_get=result.get()
            for subloc in result_get:
                opt.append('\t'.join(subloc))
                seq=self.getseq(subloc[:3], subloc[-1])
                if subloc[4]=='trf':
                    trf_seq_list.append(''.join(['>', subloc[3], '\n', seq, '\n']))
                else:
                    ltr_seq_list.append(''.join(['>', subloc[3], '\n', seq, '\n']))
        with open('retrozyme.tbl','w') as F:
            F.write('\n'.join(opt))
            F.write('\n')
        with open('ltr-retrozyme.fa','w') as F:
            F.writelines(ltr_seq_list)
        with open('trf-retrozyme.fa', 'w') as F:
            F.writelines(trf_seq_list)

def hammerhead_search(genome_file, hammerhead_description, optfile):
    with open(optfile, 'w') as rnabob_result:
        rnabob_program = subprocess.Popen(["rnabob", "-c", "-q", "-F", "-s", hammerhead_description, genome_file],
                                          stderr=subprocess.DEVNULL, stdout=rnabob_result)
        rnabob_program.wait()

def revise_rnabob_opt(rnabob_opt, topology_describe, ignore_list):
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
        lh_index = topology_list.index(lh)    ## the index for left helix
        rh_index = topology_list.index(lh + "'") ## the index for right helix

        left_drop, right_drop = [], []
        if lh_index + 2 == rh_index:   ### to find the single loop that between paired helixs
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
            if ls_name in ignore_list:   #### specifially for hh1.minimal
                ## To obtain the location of key motifs that should not be paired.
                constraint1_loc = [(i.start(), i.end()) for i in re.finditer('CTGA.GA', ls_seq)]
                constraint2_loc = [(i.start(), i.end()) for i in re.finditer('GAAA..{0,6}?.T[ATC]', rs_seq)]
                left_pair_length=constraint1_loc[0][0]
                right_pair_length=len(rs_seq)-constraint2_loc[0][1]
                left_seq_range = left_pair_length if left_pair_length < right_pair_length else right_pair_length
                ### only the 5' end of first constraint and 3'end second constraint could be paired. So only round 1 pair was considered.
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
            left_drop2, right_drop2=[],[]
            lh2_name = topology_list[lh_index + 2]  ## try to paire the single strand into another end of helix
            short_length2=short_length-len(left_drop)
            drop_length2 = 0
            if short_length2>4 and ls_name not in ignore_list: ### skip round 2 pairing if loop is constraint
                for i in range(short_length2 - 4):
                    left_nuc = ls_seq[-i-1]
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

def convert(rnabob_opt, topology_describe, ignore_list):
    rnabob_opt = rnabob_opt.strip('|').split('|')
    topology_list = topology_describe.strip().split(' ')

    full_loops = []
    for idx in range(len(topology_list) - 1):
        tpname=topology_list[idx]
        if tpname.startswith('s') and tpname not in ignore_list:
            #if topology_list[idx - 1] == topology_list[idx + 1].strip("'"):
            seq=rnabob_opt[idx]
            if len(seq)>=5:
                full_loops.append(''.join(['>',topology_list[idx],'\n',seq,'\n']))
    loop_name_dict = {}
    if full_loops:
        echo_seq_program = subprocess.Popen(["echo", '-e', "".join(full_loops)], stdout=subprocess.PIPE)
        echo_seq_program.wait()
        rnafold_program = subprocess.Popen(["RNAfold", '--noPS'], stdin=echo_seq_program.stdout, stdout=subprocess.PIPE)
        rnafold_program.wait()
        for i in rnafold_program.stdout.readlines():
            line=repr(i).replace("b'","").replace("\\n'","")
            if line.startswith('>'):
                loop_name=line.strip('>\n')
            elif re.findall('\d', line):
                loop_name_dict[loop_name]=re.split('\s+', line)[0]

    opt_topology = []
    for idx in range(len(rnabob_opt)):
        topology_name=topology_list[idx]
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

def unstrict_convert(rnabob_opt, topology_describe):
    topology_list = topology_describe.strip().split(' ')
    topology_dict = {i + 1: topology_list[i] for i in range(len(topology_list))}
    rnabob_opt = rnabob_opt.strip('|').split('|')
    opt_topology = []
    index_num = 0
    for element in rnabob_opt:
        length_element = len(element)
        if re.fullmatch('h\d+', topology_dict[index_num + 1]):
            opt_topology.append((index_num, '(' * length_element))
        elif topology_dict[index_num + 1].endswith("'"):
            opt_topology.append((index_num, ')' * length_element))
        else:
            opt_topology.append((index_num, '.' * length_element))
        index_num += 1

    opt_topology = [i[1] for i in sorted(opt_topology, key=lambda x: x[0])]
    return ''.join(opt_topology)

def rnabob_opt_to_fa(input_file, topology_describe, ignore_list, optfile):
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
        revised_seq = revise_rnabob_opt(sequence, topology_describe, ignore_list)
        if int(Args.adjust) == 1:
            dotbracket = convert(revised_seq, topology_describe, ignore_list)
        else:
            dotbracket = unstrict_convert(revised_seq, topology_describe)
        oplines.append(''.join(['>', '|'.join(location), '\n']))
        oplines.append(''.join([sequence.replace('|', ''), '\n']))
        oplines.append(''.join([dotbracket, '\n']))
    with open(optfile, 'w') as F:
        F.writelines(oplines)

def calculate_gibbs_energy(input_file):
    inter_rnaeval_file = input_file + '.rnaeval'
    with open(inter_rnaeval_file, 'w') as interfile:
        rnaeval_program = subprocess.Popen(["RNAeval", input_file],
                                           stderr=subprocess.DEVNULL, stdout=interfile)
        rnaeval_program.wait()
    oplines = []
    with open(inter_rnaeval_file, 'r') as F:
        for line in F:
            if line.startswith('>'):
                location = line.strip('>\n').split('|')
                name = location[-1]
                start, end, strand = [location[0], location[1], '+\n'] if int(location[0]) < int(location[1]) else [
                    location[1], location[0], '-\n']
            elif re.match('[a-zA-Z]', line):
                sequence = line.rstrip()
            else:
                splitlines = line.rstrip().split(' ')
                dotbracket = splitlines[0]
                energy_value = ''.join(splitlines[1:]).replace('(', '').replace(')', '').strip()
                hhm_nm = '-'.join([name.split('.')[0], start, end, 'n']) if strand == '-\n' else '-'.join(
                    [name, start, end, 'p'])     ## revised 30/08/2022
                oplines.append('\t'.join([name, hhm_nm, start, end, energy_value, sequence, dotbracket, strand]))
    return oplines

def Search_with_Motif(genome, description_hhm):
    marker = genome.replace('.fa', '').replace('.fasta', '').split('/')[-1]
    with open(description_hhm, 'r') as F:
        topology = F.readline().rstrip()
        ignor_name_list=[line.split(' ')[0] for line in F if 'CUGANGA' in line or 'GAA' in line]

    rnabob_dir = 'rnabob_opt'
    rnabob_opfile = '/'.join([rnabob_dir, '%s.rna_bob.txt' % marker])
    hammerhead_search(genome, description_hhm, rnabob_opfile)
    structure_file = '/'.join([rnabob_dir, '%s.structure.txt' % marker])
    rnabob_opt_to_fa(rnabob_opfile, topology,ignor_name_list, structure_file)
    return calculate_gibbs_energy(structure_file)

def cluster(hammerhead_tbl, opt, maxmim_gap=500):
    bedlines = []
    with open(hammerhead_tbl, 'r') as F:
        for line in F:
            splitlines = line.rstrip().split('\t')
            bedlines.append(['^_^'.join([splitlines[0], splitlines[7]]), splitlines[2], splitlines[3], splitlines[1], splitlines[4], splitlines[7]])
    bed_file_1 = hammerhead_tbl + '.bed1'
    merged_bed = hammerhead_tbl + '.merged.bed'

    bedlines = sorted(bedlines, key=lambda x: (x[0], int(x[1])))
    with open(bed_file_1, 'w') as F:
        F.write('\n'.join(['\t'.join(line) for line in bedlines]))

    with open(merged_bed, 'w') as F:
        bedtools_merge = subprocess.Popen(['bedtools', 'merge', '-i', bed_file_1, '-s'], stderr=subprocess.DEVNULL,
                                          stdout=F)
        bedtools_merge.wait()

    bed_lines = []
    with open(merged_bed, 'r') as F:
        for line in F:
            splitlines = line.rstrip().split('\t')
            subject_name, strand = splitlines[0].split('^_^')
            bed_lines.append([subject_name, splitlines[1], splitlines[2], '.', '.', strand])

    bed_lines = sorted(bed_lines, key=lambda x: (x[0], int(x[1])))

    bed_file = hammerhead_tbl + '.bed2'
    with open(bed_file, 'w') as F:
        F.write('\n'.join(['\t'.join(line) for line in bed_lines]))

    cluster_file = bed_file + '.cl'
    with open(cluster_file, 'w') as F:
        bedtools_cluster = subprocess.Popen(['bedtools', 'cluster', '-i', bed_file, '-s', '-d', str(maxmim_gap)],
                                            stderr=subprocess.DEVNULL, stdout=F)
        bedtools_cluster.wait()

    cluster_dict = defaultdict(list)
    with open(cluster_file, 'r') as F:
        for line in F:
            splilines = line.rstrip().split('\t')
            cluster_dict[splilines[-1]].append(splilines[:-1])

    OPLINES = []
    for k in cluster_dict:
        tandem_seris = cluster_dict[k]
        repeat_num = len(tandem_seris)
        if repeat_num >= 2:
            locations = sorted([(int(i[1]), int(i[2])) for i in tandem_seris], key=lambda x: x[0])
            idnumber = tandem_seris[0][0]
            strand = tandem_seris[0][-1]
            start = locations[0][0]
            end = locations[-1][-1]
            distance = ':'.join([str(locations[i + 1][0] - locations[i][1]) for i in range(repeat_num - 1)])
            location = ','.join(['-'.join([str(i[0]), str(i[1])]) for i in locations])
            OPLINES.append('\t'.join([idnumber, str(start), str(end), str(repeat_num), distance, location, strand]))
    if not OPLINES:
        sys.stdout.write('Cound not find hammerhead cluster! Exit...\n')
        sys.exit(0)
    with open(opt, 'w') as F:
        F.write('\n'.join(OPLINES))
        F.write('\n')

def sub_retro(Args):
    table_file = os.path.abspath(Args.table)
    os.chdir(WORKDIR)
    cluster(hammerhead_tbl=table_file, opt='hhm-block.tbl', maxmim_gap=int(Args.gap))
    SA = Structure_analyze(genome=GENOME, hhmtbl='hhm-block.tbl')
    SA.run_process()

def motif(Args):
    descriptor_file = os.path.abspath(Args.descriptor)
    os.chdir(WORKDIR)

    OPLINES = []
    genome_dir = 'seq'

    if not os.path.exists(genome_dir):
        os.mkdir(genome_dir)

    planpool = ThreadPool(int(Args.process))

    genome_dict = SeqIO.parse(GENOME, 'fasta')
    rnabob_dir = 'rnabob_opt'
    if not os.path.exists(rnabob_dir):
        os.mkdir(rnabob_dir)
    for record in genome_dict:
        fasta_name = record.id.replace('/', '--')
        filename = ''.join([genome_dir, '/', fasta_name, '.fa'])
        record.seq = record.seq.upper().replace('X', 'N')
        SeqIO.write([record], filename, 'fasta')
        OPLINES.append(planpool.apply_async(Search_with_Motif, args=(filename, descriptor_file,)))

    planpool.close()
    planpool.join()
    hhm_table = 'hammerhead.tbl'
    with open(hhm_table, 'w') as F:
        for opt in OPLINES:
            F.writelines(opt.get())
    return hhm_table

def retrozyme(Args):
    hhm_table = motif(Args)
    cluster(hammerhead_tbl=hhm_table, opt='hhm-block.tbl', maxmim_gap=int(Args.gap))
    SA = Structure_analyze(genome=GENOME, hhmtbl='hhm-block.tbl')
    SA.run_process()

if __name__ == "__main__":
    begin_time = time.time()

    parser = argparse.ArgumentParser(description="Hammerhead and retrozyme detection.")

    sub_parser = parser.add_subparsers(help='sub-command help')

    parser_main = sub_parser.add_parser('retrozyme', help='retrozyme search help. De-novo prediction.')
    parser_main.add_argument("-g", "--genome", type=str, required=True, help="The reference genome.")
    parser_main.add_argument("-d", "--descriptor", type=str, required=True, help="The descriptor file.")
    parser_main.add_argument("-p", "--process", type=int, default=2, required=False,
                             help="Number of processes to be used.")
    parser_main.add_argument("-a", "--adjust", type=int, default=1, required=False, help="Do adjustment to rnabob raw output? 1: yes; 0: no. (default: yes)")
    parser_main.add_argument("-r", "--rnafold", type=float, required=False, default=-10.0,
                             help="The maximum energy value.")
    parser_main.add_argument("-l", "--gap", type=int, required=False, default=500,
                             help="The maximum gap between intervals.")
    parser_main.add_argument("-o", "--opdir", type=str, required=True, help="The output file.")
    parser_main.set_defaults(func=retrozyme)

    parser_search = sub_parser.add_parser('motif', help='hammerhead motif search help.')
    parser_search.add_argument("-g", "--genome", type=str, required=True, help="The reference genome.")
    parser_search.add_argument("-d", "--descriptor", type=str, required=True, help="The descriptor file.")
    parser_search.add_argument("-p", "--process", type=int, default=2, required=False,
                               help="Number of processes to be used.")
    parser_search.add_argument("-a", "--adjust", type=int, default=1, required=False, help="Do adjustment to rnabob raw output? 1: yes; 0: no. (default: yes)")
    parser_search.add_argument("-o", "--opdir", type=str, required=True, help="The output file.")
    parser_search.set_defaults(func=motif)

    parser_retrozyme = sub_parser.add_parser('sub_retro',
                                             help="retrozyme search help. Starts after the ending of 'motif' sub-program.")
    parser_retrozyme.add_argument("-t", "--table", type=str, required=True,
                                  help="The hammerhead motif table (output of 'motif' sub-program).")
    parser_retrozyme.add_argument("-g", "--genome", type=str, required=True, help="The reference genome.")
    parser_retrozyme.add_argument("-r", "--rnafold", type=float, required=False, default=-10.0,
                                  help="The maximum energy value.")
    parser_retrozyme.add_argument("-l", "--gap", type=int, required=False, default=500,
                                  help="The maximum gap between intervals.")
    parser_retrozyme.add_argument("-p", "--process", type=int, default=2, required=False,
                                  help="Number of processes to be used.")
    parser_retrozyme.add_argument("-o", "--opdir", type=str, required=True, help="The output file.")

    parser_retrozyme.set_defaults(func=sub_retro)

    Args = parser.parse_args()

    GENOME = os.path.abspath(Args.genome)
    WORKDIR = os.path.abspath(Args.opdir)
    if not os.path.exists(WORKDIR):
        os.mkdir(WORKDIR)

    Args.func(Args)
    over_time = time.time()
    total_time = round((over_time - begin_time) / 3600, 4)
    print(total_time, ' hours used for hammerhead searching and/or retrozyme detection.\n')