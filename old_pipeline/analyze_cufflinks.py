mysql_host = "littman60"
mysql_user = "root"
mysql_pass = "aids"

_cuffcomp_tx_table = "cuffcomp_transcripts"
_cuffcomp_exon_table = "cuffcomp_exons"
_cuffcomp_loci_table = "cuffcomp_loci"

# sequence file access options
_maf_species = "mm9"
_maf_dir = "/comp_sync/data/foreign/ucsc/20100921_multiz30way"
_seq_dir = "/comp_sync/data/foreign/ucsc/20100816_mm9_sequence"
_phast_dir = "/comp_node0/andrew/work/biolibrary_data/binary_phastcons_elements"

# multiprocessing options
_max_children = 32
_batch_size = 1000

import sys
import os
import time
import itertools
import warnings

import MySQLdb

from Andrew import _fetch_sequence, _find_orfs

# returns a MySQL dict cursor
def _mysql_db_connect(mysql_db):
    mysql_conn = MySQLdb.connect(host = mysql_host,
                                 user = mysql_user,
                                 passwd = mysql_pass,
                                 db = mysql_db)

    return mysql_conn.cursor(MySQLdb.cursors.DictCursor)

# split task into batches and feed 'em to a process pool
def _tx_multiprocess(mysql_db, test_func, args = []):
    from multiprocessing import JoinableQueue, Process
    import Queue
        
    all_my_work = JoinableQueue(_max_children)

    # launch children
    sys.stderr.write("Launching children... ")
    
    all_my_children = []
    
    for i in xrange(0, _max_children):
        child = Process(target = test_func, args = [all_my_work] + args)
        child.daemon = True
        child.start()
        
        all_my_children.append(child)

    sys.stderr.write("done!\n")
        
    # fetch transcripts
    sys.stderr.write("Fetching transcripts... ")
    
    db_conn = _mysql_db_connect(mysql_db)
    db_conn.execute("SELECT transcript_id, chrom, strand, exonStarts, exonEnds FROM %s" % (_cuffcomp_tx_table,))

    fetched = db_conn.fetchall()
    
    sys.stderr.write("done!\n")

    # fill the queue
    queued_recs = 0

    for i in xrange(0, len(fetched), _batch_size):
        batch = list(fetched[i:i + _batch_size])

        all_my_work.put(batch)
        queued_recs += len(batch)

        perc_complete = (float(queued_recs) / float(len(fetched))) * 100
            
        sys.stderr.write("\r%-45s [%-25s] %3d%%" % ("Populating queue...",
                                                    "*" * int(perc_complete / 4),
                                                    perc_complete))

    # wait for work to be finished
    sys.stderr.write("\rWaiting for workers to finish... ")

    all_my_work.close()
    all_my_work.join()

    sys.stderr.write("done!\n")
 
# CSF Score
def calc_csf(mysql_db, win_len, win_overlap):
    _tx_multiprocess(mysql_db, _calc_csf, [win_len, win_overlap])
    
def _calc_csf(queue, win_len, win_overlap):
    from Bio.SeqUtils.CSFScore import CSFMultizMaf

    csf = CSFMultizMaf(_maf_dir, _maf_species)
    win_len = int(win_len)
    win_overlap = int(win_overlap)
    csf.load_matrix("/tmp/matrix.pickle")
        
    while True:
        batch = queue.get()

        for tx in batch:
            try:
                score = csf.sliding_window(tx["chrom"],
                                   tx["strand"],
                                   0,
                                   0,
                                   map(int, tx["exonStarts"].split(",")[:-1]),
                                   map(int, tx["exonEnds"].split(",")[:-1]),
                                   win_len,
                                   win_overlap)
                                   
                print tx["transcript_id"], score
            except ValueError:
                warnings.warn("Caught ValueError processing %s" % tx["transcript_id"]) 

        queue.task_done()     

# RFC Score
def calc_rfc(mysql_db, win_len, win_overlap):
    _tx_multiprocess(mysql_db, _calc_rfc, [win_len, win_overlap])
    
def _calc_rfc(queue, win_len, win_overlap):
    from Bio.SeqUtils.RFCScore import RFCMultizMaf

    rfc = RFCMultizMaf(_maf_dir, _maf_species)
    win_len = int(win_len)
    win_overlap = int(win_overlap)

    while True:
        batch = queue.get()
        
        for tx in batch:
            print tx["transcript_id"]
            try:
                rfc.sliding_window(tx["chrom"],
                                   tx["strand"],
                                   0,
                                   0,
                                   map(int, tx["exonStarts"].split(",")[:-1]),
                                   map(int, tx["exonEnds"].split(",")[:-1]),
                                   win_len,
                                   win_overlap)
            except ValueError:
                warnings.warn("Caught ValueError processing %s" % tx["transcript_id"])   

        queue.task_done()

# find longest orf in transcript
def longest_orf(mysql_db):
    _tx_multiprocess(mysql_db, _longest_orf)

def _longest_orf(queue):
    while True:
        batch = queue.get()

        for tx in batch:
            try:
                seq = _fetch_sequence(tx["chrom"],
                                      map(int, tx["exonStarts"].split(",")[:-1]),
                                      map(int, tx["exonEnds"].split(",")[:-1]),
                                      _seq_dir)
                                     
                _find_orfs(seq, True)
            except ValueError:
                warnings.warn("Caught ValueError processing %s" % tx["transcript_id"])            
            
        queue.task_done()
        
# phastcons conservation
def phastcons(mysql_db):
    _tx_multiprocess(mysql_db, _phastcons)
    
def _phastcons(queue):
    from Andrew.PhastCons import PhastConsMmap
    
    phast = PhastConsMmap(_phast_dir, 2500)
    
    while True:
        batch = queue.get()
        
        for tx in batch:
            try:
                score = phast.zscore(tx["chrom"],
                                      map(int, tx["exonStarts"].split(",")[:-1]),
                                      map(int, tx["exonEnds"].split(",")[:-1]))
            except ValueError:
                warnings.warn("Caught ValueError processing %s" % tx["transcript_id"])
            
        queue.task_done()

if __name__ == "__main__":
    #TODO replace with options parser
    if sys.argv[1].startswith("_"):
        raise ValueError("I won't let you touch my privates")

    try:
        run_class = globals()[sys.argv[1]] (*sys.argv[2:])
    except KeyError:
        print "usage"
