# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/20/2011"
__version__ = 0.0

from optparse import OptionParser
import os
import subprocess
import time
import _common
    
def _index_bam(bamfile):
    if not os.path.exists(bamfile):
        raise ValueError("BAM file not found")
    elif os.path.exists(bamfile + ".bai"):
        return True
        
    samtools_proc = subprocess.Popen(["samtools",
                                      "index",
                                      bamfile])
                                      
    while samtools_proc.poll() == None:
        time.sleep(1)
        
    if os.path.exists(bamfile + ".bai"):
        return True
    else:
        return False
    
def _bed_to_gtf(infile, outfile):
    out_fp = open(outfile, "w")    
    
    with open(infile, "r") as in_fp:
        for line_num, line in enumerate(in_fp):
            if not line.strip(): continue

            line_split = line.strip().split("\t")
            
            # 0 = chromosome
            # 1 = start
            # 2 = end
            # 9 = number of exons
            # 10 = end sizes
            # 11 = end starts            
            
            tx_out = [line_split[0],
                      "Scripture",
                      "transcript",
                      line_split[1],
                      line_split[2],
                      "1000",
                      ".",
                      ".",
                      "gene_id \"SCRG%06d\"; transcript_id \"SCRT%06d\";" % \
                      (line_num + 1, line_num + 1)]
                      
            out_fp.write("\t".join(map(str, tx_out)) + "\n")
                      
            for exon_num, exon_start, exon_size in zip(range(int(line_split[9])),
                                                       line_split[11].split(","),
                                                       line_split[10].split(",")):
                exon_out = [line_split[0],
                            "Scripture",
                            "exon",
                            int(line_split[1]) + int(exon_start),
                            int(line_split[1]) + int(exon_start) + int(exon_size),
                            "1000",
                            ".",
                            ".",
                            "gene_id \"SCRG%06d\"; transcript_id \"SCRT%06d\";" % \
                            (line_num + 1, line_num + 1)]
                          
                out_fp.write("\t".join(map(str, exon_out)) + "\n")
                
    out_fp.close()

def run_cufflinks(prefix, cufflinks_options, sample_reads, sample_name):
    outdir = os.path.join(options.output_dir, prefix + sample_name)

    if os.path.exists(outdir):
        raise ValueError("Output directory %s already exists!" % (outdir,))

    os.mkdir(outdir)
    
    with open(os.path.join(outdir, "cufflinks.log"), "w") as cl_log:
        cufflinks_cmd = ["cufflinks", "-q",
                         "-p", str(options.num_threads),
                         "-o", outdir,
                         "-m", str(options.fragment_mean),
                         "-s", str(options.fragment_sd)] + \
                         cufflinks_options + \
                         [sample_reads]

        cl_log.write(" ".join(cufflinks_cmd) + "\n\n")                         
                         
        cufflinks_proc = subprocess.Popen(cufflinks_cmd, stderr=cl_log)
                                        
    while cufflinks_proc.poll() == None:
        time.sleep(1)
        
def run_scripture(prefix, scripture_options, sample_reads, sample_name):
    outdir = os.path.join(options.output_dir, prefix + sample_name)

    if os.path.exists(outdir):
        raise ValueError("Output directory %s already exists!" % (outdir,))

    os.mkdir(outdir)

    for chrom in options.chroms:
        with open(os.path.join(outdir, "scripture.log"), "a") as sc_log:
            scripture_cmd = ["scripture",
                             "-alignment", sample_reads,
                             "-out", os.path.join(outdir, chrom + ".segments"),
                             "-sizeFile", os.path.join(args[0], options.assembly + ".sizes"),
                             "-chr", chrom,
                             "-chrSequence", os.path.join(args[0], chrom + ".fa")] + \
                             scripture_options
                             
            sc_log.write(" ".join(scripture_cmd) + "\n\n")    
                             
            scripture_proc = subprocess.Popen(scripture_cmd, stderr=sc_log)
                                        
        while scripture_proc.poll() == None:
            time.sleep(1)

    with open(os.path.join(outdir, "transcripts.bed"), "w") as final_fp:            
        for chrom in options.chroms:
            if not os.path.exists(os.path.join(outdir, chrom + ".segments")):
                continue
            
            with open(os.path.join(outdir, chrom + ".segments")) as chrom_fp:
                final_fp.write(chrom_fp.read())
                
    _bed_to_gtf(os.path.join(outdir, "transcripts.bed"),
                os.path.join(outdir, "transcripts.gtf"))

def run_cuffcompare(cuffcomp_options, input_gtfs):
    outdir = os.path.join(options.output_dir, "cuffcompare")
    
    if os.path.exists(outdir):
        raise ValueError("Output directory %s already exists!" % (outdir,))
        
    os.mkdir(outdir)
    
    with open(os.path.join(outdir, "cuffcompare.log"), "w") as cc_log:
        cuffcomp_cmd = ["cuffcompare", "-V",
                        "-o", os.path.join(outdir, "cuffcompare"),
                        "-s", args[0]] + \
                        ([] if not options.reference else ["-r", options.reference]) + \
                        cuffcomp_options + \
                        input_gtfs
        
        cc_log.write(" ".join(cuffcomp_cmd) + "\n\n")
        
        cuffcomp_proc = subprocess.Popen(cuffcomp_cmd, stderr=cc_log)
                                          
    while cuffcomp_proc.poll() == None:
        time.sleep(1)

def main():
    if not os.path.exists(options.output_dir):
        os.mkdir(options.output_dir)

    # run each sample separately in cufflinks transcript discovery mode
    for cond_name, cond_reads in zip(options.labels, args[1:]):
        run_cufflinks("cufflinks_", ["--min-frags-per-transfrag", "0"], cond_reads, cond_name)
        
    # make sure there's a BAM index on every input file
    for cond_reads in args[1:]:
        _index_bam(cond_reads)
        
    # run each sample separately through scripture's slow, shitty algorithm
    # and convert final BED output to GTF
    for cond_name, cond_reads in zip(options.labels, args[1:]):
        run_scripture("scripture_", [], cond_reads, cond_name)

    # run cuffcompare to generate the unique intersection of all output GTFs
    input_gtfs = []
    
    for cond_name in options.labels:
        input_gtfs.append(os.path.join(options.output_dir,
                                       "cufflinks_" + cond_name,
                                       "transcripts.gtf"))
        input_gtfs.append(os.path.join(options.output_dir,
                                       "scripture_" + cond_name,
                                       "transcripts.gtf"))
                                       
    run_cuffcompare([], input_gtfs)
    
    # re-run each sample separately in abundance calculation mode
    for cond_name, cond_reads in zip(options.labels, args[1:]):
        run_cufflinks("FPKM_", ["--min-frags-per-transfrag", "0",
                                "-G", os.path.join(options.output_dir,
                                                   "cuffcompare",
                                                   "cuffcompare.combined.gtf")] + \
                               ([] if not options.mask else ["-M", options.mask]),
                      cond_reads, cond_name)
                        
if __name__ == "__main__":
    parser = OptionParser(usage="%prog [options] <assembly> <sample1.sam> [...sampleN.sam]",
                          version="%prog " + str(__version__))
                          
    parser.add_option("-o",
                      dest="output_dir",
                      metavar="[./txome_out]",
                      default="./txome_out",
                      help="write output files to this directory")
                      
    parser.add_option("-p",
                      dest="num_threads",
                      type="int",
                      metavar="[1]",
                      default=1,
                      help="number of threads used during analysis")
                      
    parser.add_option("-L",
                      dest="labels",
                      metavar="sample1,...sampleN",
                      default=None,
                      help="comma-separated list of condition labels")
                      
    parser.add_option("-r",
                      dest="reference",
                      metavar="reference.gtf",
                      default=None,
                      help="reference annotation (i.e. refGene) GTF file")
                      
    parser.add_option("-M",
                      dest="mask",
                      metavar="mask.gtf",
                      default=None,
                      help="ignore all alignment within transcripts in this file")
                      
    parser.add_option("-m",
                      dest="fragment_mean",
                      type="int",
                      metavar="[200]",
                      default=200,
                      help="average fragment length (unpaired reads only)")
                      
    parser.add_option("-s",
                      dest="fragment_sd",
                      type="int",
                      metavar="[80]",
                      default=80,
                      help="fragment length std deviation (unpaired reads only)")
    
    options, args = parser.parse_args()
    
    if len(args) < 2:
        raise ValueError("Not enough arguments")
        
    if options.labels != None:
        options.labels = options.labels.split(",")
        
        if len(options.labels) <> len(args) - 1:
            raise ValueError("When using -L, must specify a label for every condition")
    else:
        options.labels = ["sample%s" % (x,) for x in range(len(args))]
        
    options.assembly, options.chroms = _common._find_genome_files(_common.genome_fasta[args[0]])
        
    main()
