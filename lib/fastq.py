# -*- coding: utf-8 -*-

from Bio import SeqIO
import subprocess
import _common
import operator

ILLUMINA_V14 = "A"
ILLUMINA_V18 = "B"

FASTQ_PARAMS = {ILLUMINA_V14: {"callback": "_illumina14",
                               "qual_offset": 64},
                ILLUMINA_V18: {"callback": "_illumina18",
                               "qual_offset": 33}}

def validate_reads(input_fastq):
    # figure out phred type and read length for each file
    fastq_ver = []
    phred_ver = []
    readlen = []

    for fname in input_fastq:
        this_ver = fastq_version(fname)
        this_readlen = fastq_readlen(fname)
        
        fastq_ver.append(this_ver)
        readlen.append(this_readlen)
        phred_ver.append(fastq_ver_to_phred(this_ver))

    for param_list in [fastq_ver, phred_ver, readlen]:
        if len(set(param_list)) <> 1:
            raise ValueError("Passed reads files have parameter types: %s" % (param_list,))

    return fastq_ver[0], phred_ver[0], readlen[0]

def _illumina14(description):
    meta_split = description.split(":")

    try:
        # IMPORTANT: 1 = KEEP THE READ, 0 = DISCARD THE READ
        # INTERNALLY, Y = DISCARD THE READ, N = KEEP THE READ
        if meta_split[7] == "1":
            meta_split[7] = "N"
        elif meta_split[7] == "0":
            meta_split[7] = "Y"
        else:
            raise ValueError("Filtered field must be 1/0")
    except IndexError:
        raise ValueError("Illumina v1.4 metadata format invalid")
    
    return meta_split[6], meta_split[7], ":".join(meta_split[:6])
    
def _illumina18(description):
    try:
        main_split = description.split(" ")
        meta_split = main_split[1].split(":")
    except IndexError:
        raise ValueError("Illumina v1.8+ metadata format invalid")
        
    if len(meta_split) <> 4:
        raise ValueError("Illumina v1.8+ metadata format invalid")
    
    return meta_split[0], meta_split[1], main_split[0]

def fastq_readlen(fname):
    readlens = dict(zip(range(1000), [0] * 1000))
    reads = 0
    
    for seq_rec in fast_fastq(open(fname, "r")):
        readlens[len(seq_rec.sequence)] += 1
        reads += 1
        
        if reads == 10000:
            break

    readlens = sorted(readlens.iteritems(), key=operator.itemgetter(1), reverse=True)

    return readlens[0][0]

def fastq_version(fname):
    fastq_vers = []
    
    for seq_rec in fast_fastq(open(fname, "r")):
        fastq_vers.append(fastq_desc_to_ver(seq_rec.id))

        if len(fastq_vers) == 10000:
            break

    fastq_vers = list(set(fastq_vers))
    
    if len(fastq_vers) <> 1:
        raise ValueError("Multiple FASTQ types in one file")
        
    return fastq_vers[0]

def fastq_desc_to_ver(desc):
    desc_split = desc.split(" ")
    
    if len(desc_split) == 2:
        # Illumina v1.8
        read_ver = ILLUMINA_V18
        
        meta_split = desc_split[1].split(":")
        
        if len(meta_split) <> 4:
            raise ValueError("Illumina v1.8+ metadata format invalid")
    elif len(desc_split) == 1:
        # Illumina v1.4
        read_ver = ILLUMINA_V14
        
        meta_split = desc_split[0].split(":")
        
        if len(meta_split) <> 8:
            raise ValueError("Illumina v1.4 metadata format invalid")
    else:
        raise ValueError("Could not detect FASTQ pipeline version")
        
    return read_ver

def fastq_ver_to_callback(ver):
    return globals()[FASTQ_PARAMS[ver]["callback"]]

def fastq_ver_to_phred(ver):
    try:
        offset = FASTQ_PARAMS[ver]["qual_offset"]
    except KeyError:
        raise ValueError("Specified FASTQ type has no quality offset")

    return offset

def fast_fastq(fp_in):
    recs = 0
    buf = []
    
    for line in fp_in:
        buf.append(line)
        
        if len(buf) == 4:
            recs += 1
            if recs % 1000 == 0:
                print recs, "records"
            yield fastq_record(buf)
            buf = []
    
class fastq_record():
    def __init__(self, lines):
        lines = [x.strip() for x in lines]

        self.id = lines[0][1:]
        self.sequence = lines[1]
        self.quals = lines[3]

    def raw(self):
        return "\n".join(["@%s" % (self.id,), self.sequence, "+", self.quals, ""])

def fastq_filter_trim(fp_in, qual_offset, min_qual, min_fraction, min_length):
    fastq_filter = subprocess.Popen([_common.FASTQ_FILTER,
                                     "-Q", str(qual_offset),
                                     "-q", str(min_qual),
                                     "-p", str(min_fraction)],
                                    stdin=fp_in,
                                    stdout=subprocess.PIPE)
    fastq_trimmer = subprocess.Popen([_common.FASTQ_TRIMMER,
                                      "-Q", str(qual_offset),
                                      "-t", str(min_qual),
                                      "-l", str(min_length)],
                                     stdin=fastq_filter.stdout,
                                     stdout=subprocess.PIPE)

    for seq_rec in fast_fastq(fastq_trimmer.stdout):
        yield seq_rec

def single_parser(fp_in, fp_out, callback, qual_offset, min_qual, min_len):
    qual_filtered = fastq_filter_trim(fp_in, qual_offset, min_qual,
                                      50, min_len)
    
    for seq_rec in qual_filtered:
        mate_pair, filtered, _ = callback(seq_rec.id)
        
        if mate_pair != "1":
            raise ValueError("Found mate_pair = 2 in single-end library")
        
        if filtered == "N":
            fp_out.write(seq_rec.raw())
        elif filtered == "Y":
            pass
        else:
            raise ValueError("Filtered field must be Y/N")

def paired_parser(fp_in, fp_out_left, fp_out_right, fp_out_orphans,
                  callback, qual_offset, min_qual, min_len):
    # filter reads, track readtags
    left_reads = []
    right_reads = []

    qual_filtered = fastq_filter_trim(fp_in, qual_offset, min_qual, 50,
                                      min_len)
    
    for seq_rec in qual_filtered:
        mate_pair, filtered, readtag = callback(seq_rec.id)

        if filtered == "N":
            if mate_pair == "1":
                left_reads.append(readtag)
            elif mate_pair == "2":
                right_reads.append(readtag)
            else:
                raise ValueError("Paired end field must be 1/2")
        elif filtered == "Y":
            pass
        else:
            raise ValueError("Filtered field must be Y/N")

    # determine pairs and orphans
    left_set = set(left_reads)
    right_set = set(right_reads)

    paired = left_set.intersection(right_set)
    orphans = left_set.symmetric_difference(right_set)

    # track order reads are written
    left_reads = []
    right_reads = []

    # re-filter reads
    fp_in.seek(0)
    qual_filtered = fastq_filter_trim(fp_in, qual_offset, min_qual, 50,
                                      min_len)

    for seq_rec in qual_filtered:
        mate_pair, filtered, readtag = callback(seq_rec.id)
        
        # strip paired-end information from read IDs
        seq_rec.id = readtag

        if readtag in paired:
            if mate_pair == "1":
                left_reads.append(readtag)
                fp_out_left.write(seq_rec.raw())
            elif mate_pair == "2":
                right_reads.append(readtag)
                fp_out_right.write(seq_rec.raw())
            else:
                raise ValueError("Paired end field must be 1/2")
        elif seq_rec.id in orphans:
            fp_out_orphans.write(seq_rec.raw())
            
    # verify that reads are sorted
    for left, right in zip(left_reads, right_reads):
        if left != right:
            raise ValueError("Reads are not sorted")
