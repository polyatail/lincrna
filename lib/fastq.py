# -*- coding: utf-8 -*-

from Bio import SeqIO
import subprocess
import _common
import operator
import os

ILLUMINA_V14 = "A"
ILLUMINA_V18 = "B"

FASTQ_PARAMS = {ILLUMINA_V14: {"callback": "_illumina14",
                               "qual_offset": 64},
                ILLUMINA_V18: {"callback": "_illumina18",
                               "qual_offset": 33}}

def _illumina14(description):
    meta_split = description.split(":")

    try:
        # IMPORTANT: 1 = KEEP THE READ, 0 = DISCARD THE READ
        # INTERNALLY, Y = DISCARD THE READ, N = KEEP THE READ
        if meta_split[7] == "1":
            meta_split[7] = False
        elif meta_split[7] == "0":
            meta_split[7] = True
        else:
            raise ValueError("Filtered field must be 1/0")
    except IndexError:
        raise ValueError("Illumina v1.4 metadata format invalid")
    
    return int(meta_split[6]), meta_split[7], ":".join(meta_split[:6])
    
def _illumina18(description):
    try:
        main_split = description.split(" ")
        meta_split = main_split[1].split(":")
    except IndexError:
        raise ValueError("Illumina v1.8+ metadata format invalid")
        
    if len(meta_split) <> 4:
        raise ValueError("Illumina v1.8+ metadata format invalid")

    if meta_split[1] == "Y":
        meta_split[1] = True
    elif meta_split[1] == "N":
        meta_split[1] = False
    else:
        raise ValueError("Filtered field must be Y/N")
    
    return int(meta_split[0]), meta_split[1], main_split[0]

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

def fastq_phred(fname):
    fastq_phreds = {33: 0,
                    59: 0,
                    64: 0}

    reads = 0
    
    for seq_rec in fast_fastq(open(fname, "r")):
        fastq_phreds[guess_quality_offset(seq_rec.quals)] += 1
        reads += 1

        if reads == 10000:
            break

    fastq_phreds = sorted(fastq_phreds.iteritems(), key=operator.itemgetter(1), reverse=True)

    return fastq_phreds[0][0]

def guess_quality_offset(quals):
    result = {33: 0,
              59: 0,
              64: 0}

    quals = [ord(x) for x in quals]

    for i in quals:
        if i >= 33 and i <= 74:
            result[33] += 1
            
        if i >= 59 and i <= 104:
            result[59] += 1
            
        if i >= 64 and i <= 104:
            result[64] += 1

    result = sorted(result.iteritems(), key=operator.itemgetter(1), reverse=True)

    return result[0][0]

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

def fast_fastq_offset(fp_in):
    recs = 0
    buf = []
    offset = fp_in.tell()
    
    while True:
        line = fp_in.readline()
        
        if not line:
            raise StopIteration

        buf.append(line)
        
        if len(buf) == 4:
            recs += 1
            if recs % 1000 == 0:
                print recs, "records"
            yield fastq_record(buf, offset)
            buf = []     
            offset = fp_in.tell()

def fast_fastq(fp_in):
    buf = []
     
    for line in fp_in:
        buf.append(line)
         
        if len(buf) == 4:
            yield fastq_record(buf, 0)
            buf = []

class fastq_record():
    def __init__(self, lines, offset):
        lines = [x.strip() for x in lines]

        self.id = lines[0][1:]
        self.sequence = lines[1]
        self.quals = lines[3]
        self.offset = offset

    def raw(self):
        return "\n".join(["@%s" % (self.id,), self.sequence, "+", self.quals, ""])

def fastq_filter_trim_native(seq_rec, qual_offset, min_qual, min_fraction, min_length):
    quals = [ord(x) - qual_offset for x in seq_rec.quals]
    
    perc_above_cutoff = float(sum([1 for x in quals if x >= min_qual])) / \
                        float(len(quals))

    if perc_above_cutoff * 100 < min_fraction:
        return False

    for pos, (nt, qual) in enumerate(reversed(zip(seq_rec.sequence, quals))):
        if qual < min_qual:
            continue
        else:
            break
    
    if pos > 0:
        seq_rec.sequence = seq_rec.sequence[:-pos]
        seq_rec.quals = seq_rec.quals[:-pos]

    if len(seq_rec.sequence) < min_length:
        return False
    
    return seq_rec
        
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
        
        if mate_pair != 1:
            raise ValueError("Found mate_pair = 2 in single-end library")
        
        if filtered == False:
            fp_out.write(seq_rec.raw())
        elif filtered == True:
            pass
        else:
            raise ValueError("Filtered field must be boolean")

def paired_parser(fp_in, fp_out_left, fp_out_right, fp_out_orphans,
                  callback, qual_offset, min_qual, min_len):
    # filter reads, track readtags
    left_reads = []
    right_reads = []

    qual_filtered = fastq_filter_trim(fp_in, qual_offset, min_qual, 50,
                                      min_len)
    
    for seq_rec in qual_filtered:
        mate_pair, filtered, readtag = callback(seq_rec.id)

        if filtered == False:
            if mate_pair == 1:
                left_reads.append(hash(readtag))
            elif mate_pair == 2:
                right_reads.append(hash(readtag))
            else:
                raise ValueError("Paired end field must be 1/2")
        elif filtered == True:
            pass
        else:
            raise ValueError("Filtered field must be boolean")

    # determine pairs and orphans
    left_set = set(left_reads)
    del left_reads
    right_set = set(right_reads)
    del right_reads

    paired = left_set.intersection(right_set)
    orphans = left_set.symmetric_difference(right_set)

    del left_set, right_set

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

        if hash(readtag) in paired:
            if mate_pair == 1:
                left_reads.append(readtag)
                fp_out_left.write(seq_rec.raw())
            elif mate_pair == 2:
                right_reads.append(readtag)
                fp_out_right.write(seq_rec.raw())
            else:
                raise ValueError("Paired end field must be 1/2")
        elif hash(readtag) in orphans:
            fp_out_orphans.write(seq_rec.raw())
            
    # verify that reads are sorted
    for left, right in zip(left_reads, right_reads):
        if left != right:
            raise ValueError("Reads are not sorted")

def split_paired_parser(fp_in, fp_out_left, fp_out_right, fp_out_orphans,
                        callback, qual_offset, min_qual, min_len):
    fp_fastq = fast_fastq_offset(fp_in)

    # guess that mate-pairs start at 50% through the file
    fp_in.seek(0)
    read1 = fp_fastq.next()

    fp_in.seek(os.stat(fp_in.name).st_size / 2)
    read2 = fp_fastq.next()

    if callback(read1.id)[2] == callback(read2.id)[2]:
        pair1_offset = 0
        pair2_offset = read2.offset
    else:
        # find out where the mate-pairs starts
        pair1_offset = 0
        pair2_offset = 0
    
        for seq_rec in fp_fastq:
            mate_pair, filtered, readtag = callback(seq_rec.id)
            
            if mate_pair == 2:
                pair2_offset = seq_rec.offset
                break
        
    while True:
        fp_in.seek(pair1_offset)
        try:
            pair1_read = fp_fastq.next()
        except StopIteration:
            break
        pair1_offset = pair1_read.offset + len(pair1_read.raw())
        
        mate_pair1, filtered1, readtag1 = callback(pair1_read.id)
        
        fp_in.seek(pair2_offset)
        try:
            pair2_read = fp_fastq.next()
        except StopIteration:
            break
        pair2_offset = pair2_read.offset + len(pair2_read.raw())
        
        mate_pair2, filtered2, readtag2 = callback(pair2_read.id)

        assert readtag1 == readtag2

        pair1_read.id = readtag1
        pair2_read.id = readtag2
        
        if filtered1 == False:
            filtered_read1 = fastq_filter_trim_native(pair1_read, qual_offset,
                                                      min_qual, 50, min_len)
            
            if filtered_read1 == False:
                filtered1 = True
            else:
                pair1_read = filtered_read1
                
        if filtered2 == False:
            filtered_read2 = fastq_filter_trim_native(pair2_read, qual_offset,
                                                      min_qual, 50, min_len)
            
            if filtered_read2 == False:
                filtered2 = True
            else:
                pair2_read = filtered_read2
                
        if filtered1 == True and filtered2 == True:
            # do nothing
            pass
        elif filtered1 == False and filtered2 == True:
            # left is an orphan
            fp_out_orphans.write(pair1_read.raw())
        elif filtered1 == True and filtered2 == False:
            # right is an orphan
            fp_out_orphans.write(pair2_read.raw())
        elif filtered1 == False and filtered2 == False:
            # they form a good pair            
            fp_out_left.write(pair1_read.raw())
            fp_out_right.write(pair2_read.raw())
        