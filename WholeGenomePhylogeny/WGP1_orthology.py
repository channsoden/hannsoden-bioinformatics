#!/global/scratch/hannsode/pkgs/anaconda2/bin/python
# Standard modules
import sys, os, shutil, pickle
from bisect import bisect_right

# Nonstandard modules
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet

# My modules
import fasta_tools, mummer_tools
from processing_tools import mapPool
from SLURM_tools import submit, job_wait

def orthology(args):
    # Find orthologous sequences between all genomes and write them to fasta files.

    # Align any queries to the reference that have not yet been aligned.
    deltafiles = runmummer(args)

    job = '{} {} {}'.format(sys.executable, __file__, args.output+'.args.pickle')
    #ID = submit(job,
    #            partition = 'savio_bigmem',
    #            account = 'co_rosalind',
    #            qos = 'savio_lowprio',
    #            time='12:0:0',
    #            job_name = 'py_ortho',
    #            cpus_per_task = 20,
    #            mem_per_cpu = '25G',
    #            modules = ['biopython/1.64-2.7.8',
    #                       'pandas/0.14.1'])
    ID = submit(job,
                partition = 'savio',
                account = 'co_rosalind',
                qos = 'rosalind_savio_normal',
                time='12:0:0',
                job_name = 'py_ortho',
                cpus_per_task = 20,
                mem_per_cpu = '3000',
                modules = ['biopython/1.64-2.7.8',
                           'pandas/0.14.1'])
    job_wait(ID)
    out_file = 'py_ortho_'+str(ID)+'.out'
    err_file = 'py_ortho_'+str(ID)+'.err'

    if not os.path.isfile('1_orthology/'+args.output+'.OS.txt'):
        exit_message = 'Orthology search failed.\n\n{}:\n{}'
        exit_message = exit_message.format(err_file, open(err_file, 'r').read())
        sys.stderr.write(exit_message)
        cleanup([out_file, err_file])
        return None

    cleanup([out_file, err_file])
    return get_segments(args)

def runmummer(args):
    basedir = os.getcwd()
    try:
        os.mkdir('1_orthology')
    except OSError:
        pass
    os.chdir('1_orthology')

    queries = [genome for genome in args.genomes if not genome.split('/')[-1] == args.reference.split('/')[-1]]
    mummer_runs = [[nucmer, [args.reference, genome]] for genome in queries]
    delta_files = mapPool(len(mummer_runs), mummer_runs)

    os.chdir(basedir)
    return delta_files
    
def nucmer(reference, query):
    prefix = os.path.splitext(os.path.basename(reference))[0] + '.' + os.path.splitext(os.path.basename(query))[0]
    deltafile = prefix + '.delta'
    tempfile = prefix + '.mgaps'
    filterfile = deltafile + '.filtered'
    while (not os.path.isfile(filterfile) or
           os.stat(filterfile).st_size == 0 or
           os.path.isfile(tempfile)):
        command_chain = ('nucmer -g 2000 --prefix={0} {1} {2}; '+
                         'delta-filter -1 {3} -o 0 > {4};')
        command_chain = command_chain.format(prefix, reference, query, deltafile, filterfile)
        ID = submit(command_chain,
                    partition='savio',
                    account='co_rosalind',
                    qos='rosalind_savio_normal',
                    time='2:0:0',
                    job_name = 'mummer',
                    mem_per_cpu = '2G')
        job_wait(ID)
        outfile = 'mummer_'+str(ID)+'.out'
        errfile = 'mummer_'+str(ID)+'.err'
        if not os.path.isdir('logs'):
            os.mkdir('logs')
        os.rename(outfile, 'logs/'+outfile)
        os.rename(errfile, 'logs/'+errfile)
        
    if os.path.isfile(deltafile):
        os.remove(deltafile)
        
    return filterfile

def merge_endpoints(endpointsA, endpointsB):
    for scaf, ends in endpointsB.items():
        try:
            endpointsA[scaf].extend(ends)
        except KeyError:
            endpointsA[scaf] = ends
    return endpointsA

def get_endpoints(df):
    with open(df, 'r') as fh:
        query = fh.readline().split()[1].strip()

    endpoints = {}
    for aln in mummer_tools.parse_delta(df):
        start = min(aln.rStart, aln.rEnd)
        end = -max(aln.rStart, aln.rEnd) # ends are arbitrarily coded negative
        try:
            endpoints[aln.rSeq].append(start)
            endpoints[aln.rSeq].append(end) 
        except KeyError:
            endpoints[aln.rSeq] = [start, end]
    
    return endpoints

def full_coverage(endpoints, max_coverage):
    segments = {}
    for scaf, ends in endpoints.items():
        ends.sort(key=abs)
        covered = 0
        segments[scaf] = []
        for i, end in enumerate(ends):
            covered += -2 * (end<0) + 1
            if covered == max_coverage:
                segments[scaf].append((end, abs(ends[i+1])))
                
    # filter to segments that are longer than 60 bp
    for scaf, coords in segments.items():
        segments[scaf] = [coord for coord in coords if coord[1] - coord[0] >= 60]

    return segments

def extract_orthologous_sequence(reference, delta_files, uni_shared_ref):
    try:
        os.mkdir(args.output+'.shared_alignment')
    except OSError:
        pass

    template = '{}_{}-{}' # scaffold_startPos-endPos
    results = open(args.output+'.OS.txt', 'w')
    
    # Write the reference segments to files
    ref = fasta_to_dict(reference)
    for scaf, seglist in uni_shared_ref.items():
        for segStart, segEnd in seglist:
            segID = template.format(scaf, segStart, segEnd)
            faFH = open(args.output+'.shared_alignment/{}.{}.fa'.format(args.output, segID), 'w')
            rec = SeqRecord(ref[scaf].seq[segStart-1:segEnd],
                            id=segID,
                            description=args.reference+':'+segID)
            SeqIO.write(rec, faFH, 'fasta')
            results.write('{}.{}.fa\n'.format(args.output, segID))
            faFH.close()

    for df in delta_files:
        # Load the query genome
        qerFile = open(df, 'r').readline().split()[1].strip()
        qer = fasta_to_dict(qerFile)

        for aln in mummer_tools.parse_delta(df):
            # Find all uni segments that are covered by this alignment
            seg = (aln.rStart, aln.rEnd)
            starts_after = bisect_right(uni_shared_ref[aln.rSeq], seg)
            if starts_after > 0 and uni_shared_ref[aln.rSeq][starts_after-1][1] > seg[0]:
                # if there are no uni segments to this rSeq, then starts_after == 0
                starts_after -= 1
            rseg = (seg[1], 'end')
            ends_before = bisect_right(uni_shared_ref[aln.rSeq], rseg)
            overlapping = uni_shared_ref[aln.rSeq][starts_after:ends_before]
            
            for rStart, rEnd in overlapping:
                # Trim the uni segments on the ends down to only the part covered by this alignment
                tStart = max(seg[0], rStart)
                tEnd = min(seg[1], rEnd)
                # Map the reference coordinates to query coordinates
                qStart, qEnd = aln_adjust(tStart, tEnd, aln)

                # Create the sequence record
                segID = template.format(aln.rSeq, tStart, tEnd)
                desc = ('{}:'+template).format(qerFile, aln.qSeq, qStart, qEnd)
                seq = qer[aln.qSeq].seq[qStart-1:qEnd]
                if not aln.orientation:
                    seq = seq.reverse_complement()
                rec = SeqRecord(seq,
                                id=segID,
                                description=desc)
                
                faFile = args.output+'.shared_alignment/'+args.output+'.'+template.format(aln.rSeq, rStart, rEnd)+'.fa'

                if tStart == rStart and tEnd == rEnd:
                    faFH = open(faFile, 'a')
                    SeqIO.write(rec, faFH, 'fasta')
                else:
                    # This alignment only covers part of the segments, and must be concatenated with another alignment.
                    # Most likely due to a rearangment in the segment.
                    records = merge_records(rec, qerFile, qStart, qEnd, faFile, tStart, tEnd, seq, aln, template)
                    faFH = open(faFile, 'w')
                    SeqIO.write(records, faFH, 'fasta')
                faFH.close()
    
    return args.output+'.shared_alignment/'+args.output+'.'+template.format(*['*']*3)+'.fa'

def fasta_to_dict(fastafile):
    with open(fastafile, 'r') as fh:
        return SeqIO.to_dict(SeqIO.parse(fh, "fasta"))

def aln_adjust(start, end, aln):
    # Maps a segment in the reference to a segment in the query using the alignment

    # Parse the delta coordinates to find the query based start and end positions.
    a = [aln.rStart] # alignment positon. Note that these are 1 based positions from MUMmer, not python indeces.
    i = [0] # insertions up to corresponding index of alignment position
    d = [0] # deletions up to corresponding index of alignemnt position
    for dig in aln.delta:
        a.append(a[-1] + abs(dig))
        if dig > 0:
            # dig is an insertion
            i.append(i[-1] + 1)
            d.append(d[-1])
        else:
            # dig is a deletion
            i.append(i[-1])
            d.append(d[-1] + 1)
                    
    startIdx = bisect_right(a, start) -1
    endIdx = bisect_right(a, end) -1

    if aln.orientation:
        # If start == aln.rStart, then (start - aln.rStart + d[] - i[]) will = 0. Start at the beggining of the alignment.
        # If start > aln.rStart, then the segment starts partway through the alignment. Add this distance adjusted for any insertions or deletions in this distance.
        qStart = aln.qStart + (start - aln.rStart + d[startIdx] - i[startIdx])
        # If end == aln.rEnd, then aln.rEnd - aln.rStart = rLen. Adjust for all indels and you get qLen.
        qEnd = aln.qStart + (end - aln.rStart +d[endIdx] - i[endIdx])
    else:
        # The alignment is to the reverse strand on the query, so work backwards.
        qEnd = aln.qStart - (start - aln.rStart + d[startIdx] - i[startIdx])
        qStart = aln.qStart - (end - aln.rStart +d[endIdx] - i[endIdx])
    
    return min(qStart, qEnd), max(qStart, qEnd)

def merge_records(rec, qerFile, qStart, qEnd, faFile, tStart, tEnd, seq, aln, template):
    records = [rec]
    for recB in SeqIO.parse(faFile, 'fasta'):
        sub_rSegs, sub_qSegs = recB.description.split()
        recQer, sub_qSegs = sub_qSegs.split(':')
        if recQer == qerFile:
            # parse the record description of previously written sub-segments
            sub_rScafs, sub_rSegs = map(list, zip(*[ss.split('_') for ss in sub_rSegs.split('|')]))
            sub_rSegs = [map(int, ss.split('-')) for ss in sub_rSegs]
            sub_qScafs, sub_qSegs = map(list, zip(*[ss.split('_') for ss in sub_qSegs.split('|')]))
            sub_qSegs = [map(int, ss.split('-')) for ss in sub_qSegs]
            i = bisect_right(sub_rSegs, [tStart, tEnd])
            
            # insert the sequence in the right place in recB
            before = sum([1+e-s for s,e in sub_qSegs[:i]])
            seq = recB.seq[:before] + seq + recB.seq[before:]
            sub_rScafs = sub_rScafs[:i] + [aln.rSeq] + sub_rScafs[i:]
            sub_rSegs = sub_rSegs[:i] + [[tStart, tEnd]] + sub_rSegs[i:]
            sub_qScafs = sub_qScafs[:i] + [aln.qSeq] + sub_qScafs[i:]
            sub_qSegs = sub_qSegs[:i] + [[qStart, qEnd]] + sub_qSegs[i:]
            
            # reform the new description strings
            segID = '|'.join([template.format(l, *v) for l,v in zip(sub_rScafs, sub_rSegs)])
            desc = qerFile+':' + '|'.join([template.format(l, *v) for l,v in zip(sub_qScafs, sub_qSegs)])
            
            # overwrite the old record with the new one
            records[0]  = SeqRecord(seq,
                                    id=segID,
                                    description=desc)
        else:
            records.append(recB)
    return records

def seg_stats(seg_list):
    lengths = [end-(start-1) for scaf_list in seg_list.values() for start, end in scaf_list]
    n = len(lengths)
    if not n:
        return 0, 0, 0
    total = sum(lengths)
    mean = float(total) / n
    return n, mean, total

def get_segments(args):
    fh = open('1_orthology/'+args.output+'.OS.txt', 'r')
    segment_files = [line.strip() for line in fh if line.strip()]
    return segment_files

def cleanup(logs):
    if not os.path.isdir('1_orthology/logs'):
        os.mkdir('1_orthology/logs')
    for log in logs:
        os.rename(log, '1_orthology/logs/'+log)

if __name__ == '__main__':
    with open(sys.argv[1], 'rb') as fh:
        args = pickle.load(fh)

    # Get the names of the delta files
    delta_files = runmummer(args)

    basedir = os.getcwd()
    os.chdir('1_orthology')

    print 'Finding segments of the reference that have orthologous sequences in all queries (uni_shared_ref). . .'
    # uni_shared_ref = {ref_scaf: [(start, end), (start, end), . . .], ref_scaf: [. . .], . . .}
    jar = args.output+'.uni_shared_ref.pickle'
    if not os.path.isfile(jar) or not os.path.isfile(args.output+'.OS.txt'):
        ref_lens = fasta_tools.get_scaffold_lengths(args.reference)
        # Make a list of endpoints of orthologous segments for each ref scaffold.
        # The queries of endpoints are not differentiated, only starts and ends (negative).
        # endpoints = [refstart, -refend, refstart, refstart, -refend, . . .]
        endpoints = {key: [1, -value] for key, value in ref_lens.items()}
        for df in delta_files:
            endpoints = merge_endpoints(endpoints, get_endpoints(df))
        # Sort endpoints, find segments where coverage == # of queries
        uni_shared_ref = full_coverage(endpoints, len(delta_files)+1)
        pickle.dump(uni_shared_ref, open(jar, 'wb'))
    else:
        uni_shared_ref = pickle.load(open(jar, 'rb'))    

    num_segs, mean_size, total_orth = seg_stats(uni_shared_ref)
    print 'Found {} orthologous segments totalling {} bp (mean {} bp).'.format(num_segs, total_orth, mean_size)
    if not num_segs:
        sys.exit('Aborting whole genome phylogeny: no orthologous segments found.')
    elif total_orth < 10000:
        sys.exit('Aborting whole genome phylogeny: too little orthology found.')
    elif mean_size < 100:
        sys.exit('Aborting whole genome phylogeny: orthologous segments too short.')
    else:
        pass
        
    print 'Extracting and writing orthologous segments. . .'
    # Go back to the alignment files, and for each uni segment use the alignments to write a fasta file
    # containing the orthologous segments from each genome.
    extract_orthologous_sequence(args.reference, delta_files, uni_shared_ref)

    if os.path.isdir(basedir+'/2_alignment'):
        [shutil.move(args.output+'.shared_alignment/'+f, basedir+'/2_alignment/'+f) for f in os.listdir(args.output+'.shared_alignment/')]
        shutil.rmtree(args.output+'.shared_alignment/')
    else:
        os.rename(args.output+'.shared_alignment', basedir+'/2_alignment')
    
    os.chdir(basedir)
