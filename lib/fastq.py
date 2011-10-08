# -*- coding: utf-8 -*-

from Bio import SeqIO

ILLUMINA_V14 = 0x01
ILLUMINA_V18 = 0x02

FASTQ_PARAMS = {ILLUMINA_V14: {"callback": "_illumina14",
                               "quals": "phred64"},
                ILLUMINA_V18: {"callback": "_illumina18",
                               "quals": "phred33"}}

def validate_reads(input_fastq):
    # figure out phred type and read length for each file
    fastq_ver = []
    phred_ver = []
    readlen = []

    for fname in input_fastq:
        fastq_ver.append(fastq_version(fname))
        phred_ver.append(phred_version(fastq_version(fname)))
        readlen.append(fastq_readlen(fname))

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

def phred_version(fastq_type):
    try:
        return FASTQ_PARAMS[fastq_type]["quals"]
    except KeyError:
        raise ValueError("Specified FASTQ type has no quality value")

def fastq_readlen(fname):
    first_record = SeqIO.parse(fname, "fastq").next()

    return len(first_record.seq)

def fastq_version(fname):
    first_record = SeqIO.parse(fname, "fastq").next()

    return fastq_desc_to_ver(first_record.description)

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

def fastq_callback(fname):
    return globals()[FASTQ_PARAMS[fastq_version(fname)]["callback"]]
    
def single_parser(fp_in, fp_out, callback):
    for seq_rec in SeqIO.parse(fp_in, "fastq"):
        mate_pair, filtered, _ = callback(seq_rec.description)
        
        if mate_pair != "1":
            raise ValueError("Found mate_pair = 2 in single-end library")
        
        if filtered == "N":
            SeqIO.write(seq_rec, fp_out, "fastq")
        elif filtered == "Y":
            pass
        else:
            raise ValueError("Filtered field must be Y/N")
            
def paired_parser(fp_in, fp_out_left, fp_out_right, callback):      
    left_count = 0
    right_count = 0
    
    for seq_rec in SeqIO.parse(fp_in, "fastq"):
        mate_pair, filtered, readtag = callback(seq_rec.description)

        seq_rec.id = readtag
        seq_rec.name = readtag
        seq_rec.description = ""

        if filtered == "N":
            if mate_pair == "1":
                left_count += 1
                SeqIO.write(seq_rec, fp_out_left, "fastq")
            elif mate_pair == "2":
                right_count += 1
                SeqIO.write(seq_rec, fp_out_right, "fastq")
            else:
                raise ValueError("Paired end field must be 1/2")
        elif filtered == "Y":
            pass
        else:
            raise ValueError("Filtered field must be Y/N")
            
    fp_out_left.seek(0)
    fp_out_right.seek(0)
    
    left_parser = SeqIO.parse(fp_out_left, "fastq")
    right_parser = SeqIO.parse(fp_out_right, "fastq")

    if left_count <> right_count:
        raise ValueError("Left read count (%s) != right read count (%s)" % \
                         (left_count, right_count))

    for _ in range(left_count):
        left_rec = left_parser.next()
        right_rec = right_parser.next()

        if callback(left_rec.description)[2] != callback(right_rec.description)[2]:
            raise ValueError("Output reads are not in order")
