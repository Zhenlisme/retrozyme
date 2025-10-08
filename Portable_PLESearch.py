import os, re, subprocess, sys, argparse, shutil
from Bio import SeqIO
from multiprocessing.pool import ThreadPool
from collections import defaultdict

"""
This program is supposed to find and distinguish different subtypes of PLE: Neptune, Poseidon, Nematis. 
The RT gene of PLE should be in upstream (within 1000 bp) of GIY-YIG gene. (positive strand).
"""

class Homologous_search:
    def __init__(self, RT_GIY_hmm, genome, wkdir, process_num):
        self.RT_GIY_hmm = RT_GIY_hmm
        self.genome = genome
        self.genome_dict = SeqIO.parse(genome, 'fasta')
        self.genome_dict = {k.id: k.seq.upper() for k in self.genome_dict}
        self.process_num=int(process_num)
        self.wkdir = wkdir
        if not os.path.exists(self.wkdir):
            os.mkdir(self.wkdir)
            os.chdir(self.wkdir)
        else:
            os.chdir(self.wkdir)
    def hmmsearch(self, subgenome):
        orf_file = ''.join([subgenome, '.orf'])
        hmm_opt = ''.join([subgenome, '.hmmsearch.out'])
        """The index of getorf output starts from 1, not 0"""
        get_orf = subprocess.Popen(['getorf', '-sequence', subgenome, '-outseq', orf_file, '-minsize', '200', '-maxsize', '30000'],
                                   stderr=subprocess.DEVNULL)
        get_orf.wait()

        RT_dict, GIY_dict = defaultdict(list), defaultdict(list)
        RT_opline, GIY_opline = [], []
        if os.path.getsize(orf_file):   ## In case the orf file is empty.
            run_hmmsearch = subprocess.Popen(['hmmsearch', '--domtblout', hmm_opt, '--noali', '-E', '1e-5', self.RT_GIY_hmm, orf_file],
                                        stdout=subprocess.DEVNULL)
            run_hmmsearch.wait()
            os.remove(orf_file)
        else:
            return RT_opline, GIY_opline
            
        if not os.path.exists(hmm_opt):
            return RT_opline, GIY_opline
        with open(hmm_opt, 'r') as F:
            for line in F:
                if line.startswith('#'):
                    continue
                splitlines = re.split('\s+', line.rstrip())
                sub_class = splitlines[3]
                subchrname = "_".join(splitlines[0].split('_')[:-1])
                chrm_name, START = subchrname.split('startat')
                start, end = re.findall('\[(\d+)\s+-\s+(\d+)\]', line)[0]
                start = str(int(start) + int(START))
                end = str(int(end) + int(START))

                aa_start,aa_end=splitlines[19:21]
                score = splitlines[7]
                c_evalue,i_evalue = splitlines[11:13]
                if float(c_evalue) > 1e-5 or float(i_evalue) > 1e-5:
                    continue

                ## To transform amino acide coord to nucleotide coord
                if int(end) > int(start):
                    strand = '+'
                    nuc_start = int(aa_start)*3-3+int(start)
                    nuc_end = int(aa_end)*3+int(start)-1
                    orf_loc = '-'.join([start, end])
                else:
                    nuc_end = int(start) - 3*int(aa_start)+3
                    nuc_start = int(start) - 3*int(aa_end) +1
                    strand = '-'
                    orf_loc = '-'.join([end, start])
                if sub_class.startswith('RT'):
                    RT_dict[splitlines[0]].append([chrm_name, str(nuc_start), str(nuc_end), sub_class, score, strand, orf_loc])
                else:
                    GIY_dict[splitlines[0]].append([chrm_name, str(nuc_start), str(nuc_end), sub_class, score, strand, orf_loc])
        for key in RT_dict:
            rt_candidate = sorted(RT_dict[key], key=lambda x:float(x[4]))[-1] ## Select the case with highest score.
            RT_opline.append(rt_candidate)
        for key in GIY_dict:
            giy_candidate = sorted(GIY_dict[key], key=lambda x: float(x[4]))[-1]  ## Select the case with highest score.
            GIY_opline.append(giy_candidate)
        return RT_opline, GIY_opline
    def intersect(self, location1, location2, slip=0, lportion=0.0, rportion=0.0):
        location1 = sorted(location1)
        location2 = sorted(location2)
        if location1[0] - slip > location2[1] or location1[1] < location2[0] - slip:
            return False
        else:
            total_list = sorted([location1[0], location1[1], location2[0], location2[1]])
            portion1 = (total_list[2] - total_list[1] + 1) / (
                        location1[1] - location1[0] + 1)  ## how much proportion the intersected sequence occupiedonseq1
            portion2 = (total_list[2] - total_list[1] + 1) / (
                        location2[1] - location2[0] + 1)  ## how much proportion the intersected sequence occupiedonseq2
            if portion1 >= lportion and portion2 >= rportion:
                return True
            else:
                return False
    def ORF_boundary(self, sequencefile):
        orf_file = ''.join([sequencefile, '.orf'])
        # """The index of getorf output starts from 1, not 0"""
        get_orf = subprocess.Popen(['getorf', '-sequence', sequencefile, '-outseq', orf_file, '-minsize', '400'],
                                   stderr=subprocess.DEVNULL)
        get_orf.wait()

        ORF_OPT = []
        with open(orf_file, 'r') as F:
            for line in F:
                if line.startswith('>'):
                    chrm_name = "_".join(line.strip('\n').split(' ')[0].split('_')[:-1])
                    start, end = re.findall('\[(\d+)\s+-\s+(\d+)\]', line)[0]
                    ORF_OPT.append(start)
                    ORF_OPT.append(end)
        ORF_OPT = sorted(ORF_OPT, key=lambda x: int(x))
        os.remove(orf_file)
        if ORF_OPT:
            return [ORF_OPT[0], ORF_OPT[-1]]
        else:
            return ['0', '0']
    def parser_hmmsearch(self, subgenome):
        RT_opline, GIY_opline = self.hmmsearch(subgenome)
        if not RT_opline or not GIY_opline:
            return []
        RT_opline=sorted(RT_opline, key=lambda x:[x[0], int(x[1])])
        GIY_opline=sorted(GIY_opline, key=lambda x:[x[0], int(x[1])])

        GIY_bed = ''.join([subgenome, '.giy.bed'])
        RT_bed = ''.join([subgenome, '.rt.bed'])

        with open(GIY_bed, 'w') as F:
            GIY_opline='\n'.join(['\t'.join(line) for line in GIY_opline])
            F.write(GIY_opline)
        with open(RT_bed, 'w') as F:
            RT_opline = '\n'.join(['\t'.join(line) for line in RT_opline])
            F.write(RT_opline)

        def merge_rtgiy(bedfile, opfile):
            cluster_file = ''.join([bedfile, '.cluster'])
            with open(cluster_file, 'w') as F:
                merge_bed_run=subprocess.Popen(['bedtools', 'cluster', '-i', bedfile, '-d', '450', '-s'], stdout=F)
                merge_bed_run.wait()

            merge_dict=defaultdict(list)
            merge_list=[]
            with open(cluster_file, 'r') as F:
                for line in F:
                    line=line.strip().split('\t')
                    cluster=line[7]
                    merge_dict[cluster].append(line[:7])
            for cluster in merge_dict:
                ## To record the domain location.
                coordlist=merge_dict[cluster]
                coord_set = [int(i[1]) for i in coordlist]
                coord_set.extend([int(i[2]) for i in coordlist])
                coord_set = sorted(coord_set)
                start = coord_set[0]
                stop = coord_set[-1]
                ## To determain sub_class, select the case with highest score
                sub_class=sorted(coordlist, key=lambda x:float(x[4]))[-1][3]
                ## To merge the ORF location
                orf_list=sorted([int(b) for i in coordlist for b in i[-1].split('-')])
                orf_coord = '-'.join([str(orf_list[0]), str(orf_list[-1])])
                merge_list.append([coordlist[0][0], str(start), str(stop),sub_class, orf_coord, coordlist[0][5]])
            merge_list=sorted(merge_list, key=lambda x:[x[0], int(x[1])])
            with open(opfile, 'w') as F:
                merge_list = ['\t'.join(i) for i in merge_list]
                F.write('\n'.join(merge_list))
                F.write('\n')
            if merge_list:
                signal = 1
            else:
                signal = 0
            os.remove(cluster_file)
            return signal

        merge_giy = ''.join([subgenome, '.giy.merge.bed'])
        merge_rt = ''.join([subgenome, '.rt.merge.bed'])

        giy_signal = merge_rtgiy(GIY_bed, merge_giy)
        rt_signal = merge_rtgiy(RT_bed, merge_rt)
        if not giy_signal or not rt_signal:  ## Either hel or rep data is null
            return []

        window_opt =  ''.join([subgenome, '.rt_giy.window.bed'])
        with open(window_opt, 'w') as window_F:
            bedtools_intersect1=subprocess.Popen(['bedtools', 'window','-a', merge_rt, '-b', merge_giy, '-w', '1500','-sm','-sw'], stdout=window_F)
            bedtools_intersect1.wait()

        #os.remove(merge_giy)
        #os.remove(merge_rt)
        #os.remove(GIY_bed)
        #os.remove(RT_bed)

        opseq=[]
        bedlist=[]
        with open(window_opt, 'r') as F:
            for line in F:
                splitlines=line.rstrip().split('\t')
                if splitlines[3]=='RT_TERT':  ## ignore TERT domain
                    continue
                strand = splitlines[5]
                chrm=splitlines[0]
                RT_start, RT_end = splitlines[1:3]
                GIY_start, GIY_end = splitlines[7:9]

                if self.intersect([int(RT_start), int(RT_end)], [int(GIY_start), int(GIY_end)],
                                  lportion=0.2, rportion=0.2):   ### If the helicase and rep domain have a intersection, skip
                    continue
                RT_orf = splitlines[4].split('-')
                GIY_orf = splitlines[10].split('-')
                loc=sorted([RT_orf[0], RT_orf[1], GIY_orf[0], GIY_orf[1]], key=lambda x:int(x))
                start,end=loc[0],loc[-1]
                maximum_length=len(self.genome_dict[chrm])
                bedlist.append([chrm, str(start), str(end), '-'.join([RT_start, RT_end]), '-'.join([GIY_start, GIY_end]),
                                strand, splitlines[3], splitlines[9]])
        #bedlist=sorted(bedlist, key=lambda x:[x[0], int(x[1])])
        #os.remove(window_opt)
        sys.stdout.write('Find %s RT-GIY blocks in %s\n' % (str(len(bedlist)), os.path.basename(subgenome).replace('.fa','')))
        return bedlist
    def split_genome(self, chunk_size=200000000, flanking_size=20000, num_groups = 2):
        if not os.path.exists('genomes'):
            os.mkdir('genomes')
        subgenome_list = []
        for chrm in self.genome_dict:
            seq_len = len(self.genome_dict[chrm])
            if seq_len < 1000:  ##Skip chrms whose length is shorter than 1000 bp
                sys.stdout.write(
                    "Chrm %s will not be used to detect autonomous Helitron/Helentron as its length is shorter than 1000 bp\n" % chrm)
                continue
            subgenome_list.append((chrm, seq_len))

        num_groups = num_groups if num_groups <= len(subgenome_list) else len(subgenome_list)
        ###  To split the genomes into several files.
        # Calculate target sum for each group
        total_sum = sum([i[1] for i in subgenome_list])
        target_sum = total_sum / num_groups
        # Sort numbers in descending order
        numbers = sorted(subgenome_list, key=lambda x: -x[1])
        # Split numbers into groups with similar sums
        groups = [[] for i in range(num_groups)]
        group_sums = [0] * num_groups
        for number in numbers:
            # Find the group with the smallest current sum and add the number to it
            min_sum_index = group_sums.index(min(group_sums))
            groups[min_sum_index].append(number)
            group_sums[min_sum_index] += number[1]

        ## To split big chrms into smaller chunks.
        subgenome_list = []
        init_num = 1
        for subgroup in groups:
            subgenome = ''.join(['genomes/subgenome', str(init_num), '.fa'])
            init_num += 1
            with open(subgenome, 'w') as F:
                for chrminfo in subgroup:
                    chrid = chrminfo[0]
                    seq = self.genome_dict[chrid]
                    seq_len = len(seq)
                    num = seq_len // chunk_size
                    for i in range(num + 1):
                        start, stop = i * chunk_size, (i + 1) * chunk_size + flanking_size
                        if start >= seq_len:
                            continue
                        if stop > seq_len:
                            stop = seq_len
                        subchrm = 'startat'.join([chrid, str(start)])
                        chunk_seq = str(self.genome_dict[chrid][start:stop])
                        F.write(''.join(['>', subchrm, '\n']))
                        F.write(chunk_seq)
                        F.write('\n')
            subgenome_list.append(subgenome)
        return subgenome_list
    def run_multiple_threads(self):
        subgenome_list = self.split_genome(chunk_size=200000000, flanking_size=20000, num_groups=200)
        if len(subgenome_list) < self.process_num:
            processnum = len(subgenome_list)
        else:
            processnum = self.process_num

        planpool = ThreadPool(processnum)
        run_result = []
        for subgenome in subgenome_list:
            run_result.append(planpool.apply_async(self.parser_hmmsearch, args=(subgenome,)))
        planpool.close()
        planpool.join()
        PLE_list=[]
        for result in run_result:
            result_get = result.get()
            PLE_list.extend(result_get)
        PLE_list=sorted(PLE_list, key=lambda x:[x[0], int(x[1])])
        with open('PLE.tbl', 'w') as F:
            F.write('\n'.join(['\t'.join(i) for i in PLE_list]))
            F.write('\n')
        #shutil.rmtree('genomes/')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PLE-like detection.")
    parser.add_argument("-hmm", "--hmm", type=str, required=True, help="The PLE protein hmm model file.")
    parser.add_argument("-g", "--genome", type=str, required=True, help="The reference genome.")
    parser.add_argument("-o", "--opdir", type=str, required=True, help="The output directory.")
    parser.add_argument("-p", "--process", type=int, default=2, required=False, help="Number of threads to be used.")
    Args = parser.parse_args()
    HomoSearch = Homologous_search(os.path.abspath(Args.hmm),
                                   os.path.abspath(Args.genome),
                                   os.path.abspath(Args.opdir),
                                   Args.process)
    HomoSearch.run_multiple_threads()
    
