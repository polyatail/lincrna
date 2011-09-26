# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "9/19/2011"
__version__ = 0.0

from multiprocessing import Pool, TimeoutError
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import _common
import os
import tempfile
import time
import subprocess

species_to_name = {"mm9": "Mouse",
                   "rn4": "Rat",
                   "cavPor2": "Guinea_Pig",
                   "oryCun1": "Rabbit",
                   "hg18": "Human",
                   "panTro2": "Chimp",
                   "rheMac2": "Rhesus",
                   "otoGar1": "Bushbaby",
                   "tupBel1": "TreeShrew",
                   "sorAra1": "Shrew",
                   "eriEur1": "Hedgehog",
                   "canFam2": "Dog",
                   "felCat3": "Cat",
                   "equCab1": "Horse",
                   "bosTau3": "Cow",
                   "dasNov1": "Armadillo",
                   "loxAfr1": "Elephant",
                   "echTel1": "Tenrec"}
    
def _score_block(work_block, assembly):
    # each process loads its own MafIO instances
    chrom_maf = _common._load_maf(assembly)
    
    # for each transcript, pull the multiple alignment using MafIO and write it
    # to at tempfile. keep track of tempfiles for every transcript_id
    temps_to_tx_id = {}
    
    for transcript in work_block:
        # pull the multiple alignment of this transcript using MafIO. each
        # transcript is a list of this format:
        #     [transcript_id, chrom, exonStarts, exonEnds, strand]
        multi_align = chrom_maf[transcript[1]].get_spliced(transcript[2],
                                                           transcript[3],
                                                           "+")

        # convert species codes from UCSC's multiz alignments into species
        # names that are compatible with PhyloCSF. make sure that our species
        # of interest (given by 'assembly') appears at the top of the file
        multi_seq = []
        assembly_rec = None
                                                           
        for record in multi_align:
            species_code = record.id.split(".")[0]
            
            try:
                record.id = species_to_name[species_code]
                record.name = species_to_name[species_code]
                record.description = ""
            except KeyError:
                continue
            
            if species_code == assembly:
                assembly_rec = record
            else:
                multi_seq.append(record)
                
        multi_seq.insert(0, assembly_rec)            
        multi_align = MultipleSeqAlignment(multi_seq)
        
        # write the alignment to temporary file        
        maf_temp = tempfile.NamedTemporaryFile()        
        AlignIO.write(multi_align, maf_temp, "fasta")
        maf_temp.flush()
        
        # keep track of which temporary file correponds with each transcript
        temps_to_tx_id[maf_temp.name] = transcript[0]

    # write a list of all temporary files (containing multiple alignments) to
    # yet another temporary file. this will be passed to PhyloCSF     
    files = tempfile.NamedTemporaryFile("w", delete=False)
    files.write("\n".join(temps_to_tx_id.keys()))
    files.flush()

    phylocsf_proc = subprocess.Popen([_common.PHYLOCSF_PATH,
                                      "29mammals",
                                      "--files", files.name,
                                      "--frames=6",
                                      "--orf=StopStop3",
                                      "--minCodons=29",
                                      "--removeRefGaps"],
                                      stdout=subprocess.PIPE)

    phylocsf_out = {}

    while phylocsf_proc.poll() == None:
        time.sleep(1)
        
    for line in phylocsf_proc.stdout:
        if not line.strip():
            continue

        line_split = line.strip().split("\t")
        
        tx_id = temps_to_tx_id[line_split[0]]
        
        if line_split[1] == "abort":
            raise ValueError("Abort returned for %s: %s" % (tx_id, line_split[2]))
        elif line_split[1] == "failure":
            phylocsf_out[tx_id] = line_split[2]
        else:
            phylocsf_out[tx_id] = float(line_split[2])
        
    return phylocsf_out
    
def score_transcripts(transcripts, assembly, num_threads):
    pool = Pool(processes=num_threads)
    
    results = []
    
    # fill the queue
    work_block = []
    
    for tx_id in transcripts.keys():
        work_block.append([tx_id,
                           transcripts[tx_id]["chrom"],
                           transcripts[tx_id]["exonStarts"].values(),
                           transcripts[tx_id]["exonEnds"].values()])
                           
        if len(work_block) >= 1:
            results.append(pool.apply_async(_score_block, [work_block, assembly]))
            work_block = []
    else:
        if len(work_block) > 0:
            results.append(pool.apply_async(_score_block, [work_block, assembly]))

    # get results
    result_data = {}

    while True:
        for num, data in enumerate(results):
            try:
                this_result = data.get(timeout=0)
                
                for tx_id, csf_score in this_result.items():
#                    print tx_id, csf_score
                    result_data[tx_id] = csf_score

                del results[num]
            except TimeoutError:
                pass
            
        if len(results) == 0:
            break
            
#        print "Waiting for", len(results), "work blocks to return"
        
        time.sleep(1)
        
    return result_data
