cuffcomp_tables = {"tx": "cuffcomp_transcripts",
                   "exons": "cuffcomp_exons",
                   "loci": "cuffcomp_loci"}

cuffcomp_schemas = \
"""
CREATE TABLE IF NOT EXISTS `%(exons)s` (
  `id` bigint(20) NOT NULL auto_increment,
  `transcript_id` char(20) NOT NULL,
  `chr` char(15) NOT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  `strand` char(1) NOT NULL,
  `exon_num` int(11) NOT NULL,
  `class_code` char(1) NOT NULL,
  `nearest_ref` text,
  `gene_name` text,
  `siphy_score` float default NULL,
  PRIMARY KEY  (`id`),
  KEY `transcript_id` (`transcript_id`)
) ENGINE=MyISAM AUTO_INCREMENT=411634 DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `%(loci)s` (
  `locus_id` char(20) NOT NULL,
  `chr` char(15) NOT NULL,
  `start` int(11) NOT NULL,
  `end` int(11) NOT NULL,
  `strand` char(1) NOT NULL,
  `nearest_ref` text,
  `gene_name` text,
  PRIMARY KEY  (`locus_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `%(tx)s` (
  `transcript_id` char(20) NOT NULL,
  `locus_id` char(20) NOT NULL,
  `class_code` char(1) NOT NULL,
  `nearest_ref_1kb` text,
  `manual_annotation` enum('yes','no') NOT NULL default 'no',
  `nearest_ref` text,
  `gene_name` char(40) default NULL,
  `known_noncoding` enum('yes','no','unknown','not_in_annodb','unannotated','ambiguous') NULL,
  `num_exons` int(11) default NULL,
  `transcript_length` int(11) default NULL,
  `presplice_length` int(11) default NULL,
  `u_exon` enum('yes','no') NOT NULL,
  `antisense` enum('yes','no','strand_unknown','not_in_annodb','unannotated','ambiguous') NULL,
  `longest_orf_length` int(11) default NULL,
  `longest_orf_coords` varchar(255) default NULL,
  `phastcons_count` int(11) default NULL,
  `phastcons_zscore` float default NULL,
  `phastcons_pvalue` float default NULL,
  `rfc_score` int(11) default NULL,
  `rfc_voters` int(11) default NULL,
  `rfc_norm` float default NULL,
  `csf_score` float default NULL,
  `csf_score_20placental` float default NULL,
  PRIMARY KEY  (`transcript_id`),
  KEY `gene_name` (`gene_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE IF NOT EXISTS `rfc_longest_orf` (
  `transcript_id` varchar(40) NOT NULL,
  `mm9` float default NULL,
  `rn4` float default NULL,
  `cavPor2` float default NULL,
  `oryCun1` float default NULL,
  `hg18` float default NULL,
  `panTro2` float default NULL,
  `ponAbe2` float default NULL,
  `rheMac2` float default NULL,
  `calJac1` float default NULL,
  `otoGar1` float default NULL,
  `tupBel1` float default NULL,
  `sorAra1` float default NULL,
  `eriEur1` float default NULL,
  `canFam2` float default NULL,
  `felCat3` float default NULL,
  `equCab1` float default NULL,
  `bosTau3` float default NULL,
  `dasNov1` float default NULL,
  `loxAfr1` float default NULL,
  `echTel1` float default NULL,
  `monDom4` float default NULL,
  `ornAna1` float default NULL,
  `galGal3` float default NULL,
  `anoCar1` float default NULL,
  `xenTro2` float default NULL,
  `tetNig1` float default NULL,
  `fr2` float default NULL,
  `gasAcu1` float default NULL,
  `oryLat1` float default NULL,
  `danRer5` float default NULL,
  PRIMARY KEY  (`transcript_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
""" % cuffcomp_tables

__DEBUG_LEVEL = 5

if __debug__:
    def debug_print (msg, msg_level):
        if __DEBUG_LEVEL >= msg_level:
            print "debug%d\t%s" % (msg_level, msg)

import sys

sys.path.append ("/comp_node0/andrew/work/biolibrary")

def create_comp_tables (database):
    import functions

    db_conn = functions.mysql_db_connect (database)
    db_conn.execute (cuffcomp_schemas)

def _check_num_rows (db_conn, table, field, count):
    # check that number of rows in table = number of queries generated
    db_conn.execute ("SELECT COUNT(%s) FROM %s" % (field, table))
    rows_in_table = db_conn.fetchone()["COUNT(%s)" % (field,)]

    if count != rows_in_table:
        print "WARNING: Number of rows does not match number of queries!"
        print "\trows=%d lines=%d" % (rows_in_table, count)
    else:
        print "Check OK -- number of rows matches number of queries"

def get_transcript_sequence (assembly, database, transcript, sep_exons = "", stdout = True):
    """
    a wrapper function to call _get_transcript_seq from command line
    """

    import functions, bio_functions

    # input validation
    if assembly not in bio_functions.valid_chr:
        raise SystemError ("process_cuffcompare::get_transcript_sequence (%s, %s, %s, %s) assembly not valid" %
                           (assembly, database, transcript, sep_exons))

    db_conn = functions.mysql_db_connect (database)

    ##FIXME##
    # temporary fix to let mac access the old version's db style
    if database == "lincrna":
        table_name = "exons"
    else:
        table_name = cuffcomp_tables["exons"]

    sequence = _get_transcript_seq (db_conn, assembly, transcript, table_name, sep_exons, True)

    if sequence == False:
        raise SystemError ("process_cuffcompare::get_transcript_sequence (%s, %s, %s, %s) transcript not found" %
                           (assembly, database, transcript, sep_exons))

    if stdout == True:
        print ">%s length=%s exons=%s\n%s" % (transcript, len (sequence[1].replace(sep_exons, "")), sequence[0], sequence[1])
    elif stdout == False:
        return sequence[1]

def _get_transcript_seq (db_conn, assembly, transcript, table_name = cuffcomp_tables["exons"], sep_exons = "", num_exons = False):
    import bio_functions

    # get all exons
    db_conn.execute ("SELECT * FROM %s WHERE transcript_id = '%s' ORDER BY exon_num ASC" % (table_name, transcript))

    if db_conn.rowcount == 0:
        return False

    exon_result = db_conn.fetchall()

    transcript_seq = []

    # loop through all exons
    for exon_data in exon_result:
        # cat together exon sequences
        transcript_seq.append (bio_functions.get_sequence (assembly,
                                                           exon_data["chr"],
                                                           exon_data["start"],
                                                           exon_data["end"] + 1))

    if num_exons == True:
        return (db_conn.rowcount, sep_exons.join (transcript_seq))
    else:
        return sep_exons.join (transcript_seq)

def import_cuffcomp (database, directory, prefix, subtask = None):
    """
    a wrapper function to call various import tasks
    """

    import os

    valid_subtasks = {"exons": (".combined.gtf", cuffcomp_tables["exons"], _import_cuffcomp_GTF),
                      "transcripts": (".tracking", cuffcomp_tables["tx"], _import_cuffcomp_tracking),
                      "loci": (".loci", cuffcomp_tables["loci"], _import_cuffcomp_loci)}

    # input sanitization
    directory = os.path.realpath (directory)
    prefix = os.path.basename (os.path.realpath (prefix))

    # figure out what tasks to run
    tasks_to_do = []

    if subtask == None:
        tasks_to_do.extend (valid_subtasks.keys())
    else:
        tasks_to_do.extend (subtask.split (","))

    # input validation
    filenames = {}

    for i in tasks_to_do:
        if i not in valid_subtasks.keys():
            raise SystemError ("process_cuffcompare::import_cuffcomp (%s, %s, %s, %s) invalid subtask" % \
                               (database, directory, prefix, subtask))

        filenames[i] = "%s/%s%s" % (directory, prefix, valid_subtasks[i][0])

        if not os.path.exists (filenames[i]):
            raise SystemError ("process_cuffcompare::import_cuffcomp (%s, %s, %s, %s) %s does not exist" % \
                               (database, directory, prefix, subtask, filenames[i]))

    # do it to it, brother
    for i in tasks_to_do:
        valid_subtasks[i][2] (filenames[i], database, valid_subtasks[i][1])
        print

def _import_cuffcomp_GTF (filename, database, table_name):
    """
    import data from '<prefix>.combined.gtf' into db.cuffcomp_exons

    id, transcript_id, chr, start, end, strand, exon_num, class_code,
    nearest_ref, gene_name
    """

    import functions, bio_functions

    cuffcomp_gtf = bio_functions.ParseGTF (filename)

    db_conn = functions.mysql_db_connect (database)

    queries_generated = 0

    print "Importing data into %s.%s..." % (database, table_name)

    for data in cuffcomp_gtf.iteritems():
        if queries_generated % 300 == 0:
            functions.write_status_line("Processing %s exon %s" % (data["transcript_id"], data["exon_number"]),
                                        int(100 * (float(queries_generated) / float(cuffcomp_gtf.recordcount))))

        build_query = {}

        build_query["transcript_id"] = data["transcript_id"]
        build_query["chr"] = data["chr"]
        build_query["start"] = data["start"]
        build_query["end"] = data["end"]
        build_query["strand"] = data["strand"]
        build_query["exon_num"] = data["exon_number"]
        build_query["class_code"] = data["class_code"]

        # these two won't always be in the dict, check if they exist first
        if "nearest_ref" in data:
            build_query["nearest_ref"] = data["nearest_ref"]

        if "gene_name" in data:
            build_query["gene_name"] = data["gene_name"]

        query = functions.mysql_build_query ("INSERT INTO %s SET " % (table_name,), build_query)

        db_conn.execute (query)

        queries_generated += 1

    functions.write_status_line ("Processing complete!", 100, True)

    # check if insert was OK
    _check_num_rows (db_conn, table_name, "id", queries_generated)

    # cleanup
    db_conn.close()

def _import_cuffcomp_tracking (filename, database, table_name):
    """
    import data from '<prefix>.tracking' into db.cuffcomp_transcripts

    transcript_id, locus_id, class_code, nearest_ref, gene_name
    """

    import functions

    lines_in_file = functions.count_lines_in_file (filename)
    cuffcomp_tracking = open(filename, "r")

    db_conn = functions.mysql_db_connect (database)

    queries_generated = 0

    print "Importing data into %s.%s..." % (database, table_name)

    for line in cuffcomp_tracking:
        # should be tab-delimited: split it into a list
        data = line.split("\t")

        if queries_generated % 100 == 0:
            functions.write_status_line("Processing %s" % (data[0],),
                                        int(100 * (float(queries_generated) / float(lines_in_file))))

        build_query = {}

        build_query["transcript_id"] = data[0]
        build_query["locus_id"] = data[1]
        build_query["class_code"] = data[3]

        # can be a true nearest_ref or '-', check first
        if data[2].rstrip() != "-":
            # expect data[2] to be 'gene_name|nearest_ref', split it up
            ref_split = data[2].split("|")

            build_query["gene_name"] = ref_split[0]
            build_query["nearest_ref"] = ref_split[1]

        query = functions.mysql_build_query ("INSERT INTO %s SET " % (table_name,), build_query)

        db_conn.execute (query)

        queries_generated += 1

    functions.write_status_line ("Processing complete!", 100, True)

    # check if insert was OK
    _check_num_rows (db_conn, table_name, "transcript_id", queries_generated)

    # cleanup
    db_conn.close()
    cuffcomp_tracking.close()

def _import_cuffcomp_loci (filename, database, table_name):
    """
    import data from '<prefix>.loci' into db.cuffcomp_loci

    locus_id, chr, start, end, strand, nearest_ref, gene_name
    """

    import functions, re

    lines_in_file = functions.count_lines_in_file (filename)
    cuffcomp_loci = open(filename, "r")

    db_conn = functions.mysql_db_connect (database)

    queries_generated = 0

    print "Importing data into %s.%s..." % (database, table_name)

    for line in cuffcomp_loci:
        # should be tab-delimited, split it into a list
        data = line.split("\t")

        if queries_generated % 100 == 0:
            functions.write_status_line ("Processing %s " % (data[0],),
                                         int(100 * (float(queries_generated) / float(lines_in_file))))

        build_query = {}

        build_query["locus_id"] = data[0]

        # expect data[1] to be chr[strand]start-end
        coord_list = re.sub("[\[\]]", "\t", data[1]).split("\t")

        build_query["chr"] = coord_list[0]
        build_query["strand"] = coord_list[1]
        build_query["start"] = coord_list[2].split("-")[0]
        build_query["end"] = coord_list[2].split("-")[1]

        gene_names = []
        nearest_refs = []

        # data[2] can be 'gene_name|nearest_ref, ...' or '-', check first
        if data[2].rstrip() != "-":
            ref_split = data[2].split(",")

            for ref in ref_split:
                ref = ref.split("|")

                gene_names.append(ref[0])
                nearest_refs.append(ref[1])

            build_query["gene_name"] = ",".join(set(gene_names))
            build_query["nearest_ref"] = ",".join(set(nearest_refs))

        query = functions.mysql_build_query ("INSERT INTO %s SET " % (table_name,), build_query)

        db_conn.execute (query)

        queries_generated += 1

    functions.write_status_line ("Processing complete!", 100, True)

    # check if insert was OK
    _check_num_rows (db_conn, table_name, "locus_id", queries_generated)

    # cleanup
    db_conn.close()
    cuffcomp_loci.close()

def manual_annotation (assembly, database, force = "no"):
    """
    attempt to manually annotate transcripts based on our own logic
    """

    def find_best_match (needle_start, needle_end, haystack_start, haystack_end):
        smallest_dist = None
        best_match = None
        
        for hs_start, hs_end in zip (haystack_start, haystack_end):
            this_dist = abs (hs_start - needle_start) + abs (hs_end - needle_end)
            
            if smallest_dist == None or this_dist < smallest_dist:
                smallest_dist = this_dist
                best_match = (hs_start, hs_end, this_dist)
                
        return best_match
    
    db_priority = ("ensembl", "ucsc", "refgene", "mir")

    import functions, bio_functions, annodb

    # input validation
    if assembly not in bio_functions.valid_chr:
        raise SystemError ("process_cuffcompare::annotate_unknowns (%s, %s) assembly not valid" %
                           (assembly, database))
                           
    annodb_obj = annodb.AnnoDB (assembly)

    db_conn = functions.mysql_db_connect (database)

    db_conn.execute ("SELECT transcript_id FROM %s" % (cuffcomp_tables["tx"],))
    total_num_transcripts = db_conn.rowcount
    transcript_result = db_conn.fetchall()

    transcripts_processed = 0

    print "Manually annotating transcripts..."

    for transcript_data in transcript_result:
        if transcripts_processed % 1 == 0:
            functions.write_status_line ("Processing %s" % (transcript_data["transcript_id"],),
                                         int(100 * (float(transcripts_processed) / float(total_num_transcripts))))

        db_conn.execute ("SELECT chr, MIN(start) AS start, MAX(end) AS end FROM %s WHERE transcript_id = '%s'" % (cuffcomp_tables["exons"], transcript_data["transcript_id"]))
        exon_result = db_conn.fetchone()

        anno_info = annodb_obj.query_all_db (annodb_obj.QT_ALL_OVERLAPS, combine = True, chrom = exon_result["chr"], start = exon_result["start"], end = exon_result["end"])
        
        if len (anno_info) == 0:
            manual_nearest_ref = "NULL"
            manual_gene_name = "NULL"
        elif len (anno_info) == 1:
            manual_nearest_ref = "'" + anno_info[0]["name"] + "'"
            manual_gene_name = "'" + anno_info[0]["name2"] + "'"
        elif len (anno_info) > 1:
            # get all exons for this transcript
            db_conn.execute ("SELECT * FROM %s WHERE transcript_id = '%s' ORDER BY exon_num ASC" % (cuffcomp_tables["exons"], transcript_data["transcript_id"]))
            all_exons_result = db_conn.fetchall()
            
            exon_starts = []
            exon_ends = []
            
            for exon_data in all_exons_result:
                exon_starts.append (int (exon_data["start"]))
                exon_ends.append (int (exon_data["end"]))
                
            if __debug__: debug_print ("manual_annotation transcript_id=%s" % (transcript_data["transcript_id"],), 3)
            if __debug__: debug_print ("manual_annotation exon_starts=%s" % (exon_starts,), 3)
            if __debug__: debug_print ("manual_annotation exon_ends=%s" % (exon_ends,), 3)  

            lowest_score = None
            lowest_score_data = None         
            
            if __debug__: debug_print ("manual_annotation looking at %s annotations" % (len (anno_info),), 3)
                
            # assign a distance to each overlapping annotation
            for anno in anno_info:
                try:
                    anno_exon_starts = map (int, anno["exonStarts"].rstrip(",").split(","))
                    anno_exon_ends = map (int, anno["exonEnds"].rstrip(",").split(","))
                except KeyError:
                    anno_exon_starts = [int (anno["start"])]
                    anno_exon_ends = [int (anno["end"])]
                
                if __debug__: debug_print ("manual_annotation comparing to %s" % (anno["name"],), 3)
                if __debug__: debug_print ("manual_annotation anno_exon_starts=%s" % (anno_exon_starts,), 3)
                if __debug__: debug_print ("manual_annotation anno_exon_ends=%s" % (anno_exon_ends,), 3)
                
                unmatchable_exons = abs (len (anno_exon_starts) - len (exon_starts))
                
                if __debug__: debug_print ("manual_annotation unmatchable_exons=%s" % (unmatchable_exons,), 3)

                this_tx_score = 0    
                unmatched_exons = 0
                
                for cl_exon_start, cl_exon_end in zip (exon_starts, exon_ends):
                    if __debug__: debug_print ("manual_annotation looking at cl exon (%s, %s)" % (cl_exon_start, cl_exon_end), 3)
                    
                    cl_to_anno_best_match = find_best_match (cl_exon_start, cl_exon_end, anno_exon_starts, anno_exon_ends)
                    anno_to_cl_best_match = find_best_match (cl_to_anno_best_match[0], cl_to_anno_best_match[1], exon_starts, exon_ends)

                    if __debug__: debug_print ("manual_annotation cl_to_anno_best_match=%s" % (cl_to_anno_best_match,), 3)
                    if __debug__: debug_print ("manual_annotation anno_to_cl_best_match=%s" % (anno_to_cl_best_match,), 3)

                    if (anno_to_cl_best_match[0], anno_to_cl_best_match[1]) == (cl_exon_start, cl_exon_end):
                        # the cl exon's best match's best match is the cl exon
                        # add it to the score

                        this_tx_score += anno_to_cl_best_match[2]
                        recip_match = "yes"
                    else:
                        unmatched_exons += 1
                        recip_match = "no"
                        
                    if __debug__: debug_print ("manual_annotation recip_match=%s" % (recip_match,), 3)
                        
                if __debug__: debug_print ("manual_annotation this_tx_score=%s" % (this_tx_score,), 3)
                if __debug__: debug_print ("manual_annotation unmatched_exons=%s" % (unmatched_exons,), 3)

                if lowest_score == None or this_tx_score < lowest_score:
                    if __debug__: debug_print ("manual_annotation %s < %s, current winner" % (this_tx_score, lowest_score), 3)
                    lowest_score_data = [anno]
                    lowest_score = this_tx_score
                elif this_tx_score == lowest_score:
                    if __debug__: debug_print ("manual_annotation %s = %s, adding to winner list" % (this_tx_score, lowest_score), 3)
                    lowest_score_data.append (anno)
            
            if __debug__: debug_print ("manual_annotation lowest_score_data=%s" % (lowest_score_data,), 3)                
                
            if len (lowest_score_data) > 1:
                break_flag = False
                
                for db in db_priority:
                    if break_flag == True: break
                
                    for score_anno in lowest_score_data:
                        if score_anno["db"] == db:
                            lowest_score_data = score_anno
                            break_flag = True
                            break
            elif len (lowest_score_data) == 1:
                lowest_score_data = lowest_score_data[0]
                
            if __debug__: debug_print ("manual_annotation picked %s" % (lowest_score_data["name"],), 3)
                
            manual_nearest_ref = "'" + lowest_score_data["name"] + "'"
            manual_gene_name = "'" + lowest_score_data["name2"] + "'"
            
        db_conn.execute ("UPDATE cuffcomp_transcripts SET manual_nearest_ref = %s, manual_gene_name = %s WHERE transcript_id = '%s'" % (manual_nearest_ref, manual_gene_name, transcript_data["transcript_id"]))

        transcripts_processed += 1

    functions.write_status_line ("Processing complete!", 100, True)
    print
    
def nearest_annotation (assembly, database, force = "no"):
    """
    Calculate distance to nearest annotation for each transcript
    """

    import functions, bio_functions, annodb

    # input validation
    if assembly not in bio_functions.valid_chr:
        raise SystemError ("process_cuffcompare::analysis_first_pass (%s, %s) assembly not valid" %
                           (assembly, database))
                           
    annodb_inst = annodb.AnnoDB (assembly)

    db_conn = functions.mysql_db_connect (database)

    if force == "yes":
        dont_overwrite = ""
    elif force == "no":
        dont_overwrite = " WHERE 5p_dist IS NULL"

    db_conn.execute ("SELECT transcript_id FROM %s%s" % (cuffcomp_tables["tx"], dont_overwrite))
    total_num_transcripts = db_conn.rowcount
    transcript_result = db_conn.fetchall()

    transcripts_processed = 0

    print "Finding nearest annotations..."

    for transcript_data in transcript_result:
        if transcripts_processed % 1 == 0:
            functions.write_status_line ("Processing %s" % (transcript_data["transcript_id"],),
                                         int(100 * (float(transcripts_processed) / float(total_num_transcripts))))

        db_conn.execute ("SELECT chr, MIN(start), MAX(end) FROM %s WHERE transcript_id = '%s'" % (cuffcomp_tables["exons"], transcript_data["transcript_id"]))
        exon_result = db_conn.fetchone()
        
        anno_info_5p = annodb_inst.query_all_db (annodb_inst.QT_NEAREST_5P, combine = True, chrom = exon_result["chr"], start = exon_result["MIN(start)"], end = exon_result["MAX(end)"])
        anno_info_3p = annodb_inst.query_all_db (annodb_inst.QT_NEAREST_3P, combine = True, chrom = exon_result["chr"], start = exon_result["MIN(start)"], end = exon_result["MAX(end)"])

        build_query = {}
        
        if len (anno_info_5p) > 0:
            build_query["5p_dist"] = "'" + str (int (exon_result["MIN(start)"]) - int (anno_info_5p[0]["end"])) + "'"
            build_query["5p_nearest_ref"] = "'" + anno_info_5p[0]["name"] + "'"
        else:
            build_query["5p_dist"] = -1
            build_query["5p_nearest_ref"] = "NULL"
            
        if len (anno_info_3p) > 0:
            build_query["3p_dist"] = "'" + str (int (anno_info_3p[0]["start"]) - int (exon_result["MAX(end)"])) + "'"
            build_query["3p_nearest_ref"] = "'" + anno_info_3p[0]["name"] + "'"
        else:
            build_query["3p_dist"] = -1
            build_query["3p_nearest_ref"] = "NULL"

        db_conn.execute (functions.mysql_build_query ("UPDATE %s SET " % (cuffcomp_tables["tx"],), build_query, " WHERE transcript_id = '%s'" % (transcript_data["transcript_id"],), enclose = ""))
        
        transcripts_processed += 1

    functions.write_status_line ("Processing complete!", 100, True)
    print

def analysis_first_pass (assembly, database, force = "no"):
    """
    populates additional fields in transcript table

    num_exons, transcript_length, presplice_length, u_exon, antisense
    """

    import functions, bio_functions, annodb

    # input validation
    if assembly not in bio_functions.valid_chr:
        raise SystemError ("process_cuffcompare::analysis_first_pass (%s, %s) assembly not valid" %
                           (assembly, database))

    annodb_obj = annodb.AnnoDB (assembly)

    db_conn = functions.mysql_db_connect (database)

    if force == "no":
        dont_overwrite = " WHERE num_exons IS NULL"
    elif force == "yes":
        dont_overwrite = ""

    db_conn.execute ("SELECT manual_nearest_ref, manual_gene_name, nearest_ref, gene_name, transcript_id FROM %s%s" % (cuffcomp_tables["tx"], dont_overwrite))
    total_num_transcripts = db_conn.rowcount
    transcript_result = db_conn.fetchall()

    transcripts_processed = 0

    print "First pass transcript analysis..."

    for transcript_data in transcript_result:
        if transcripts_processed % 1 == 0:
            functions.write_status_line ("Processing %s" % (transcript_data["transcript_id"],),
                                         int(100 * (float(transcripts_processed) / float(total_num_transcripts))))

        # determine transcript length, pre-splice length, number of exons,
        # and if there's a class_code "u" exon
        db_conn.execute ("SELECT * FROM %s WHERE transcript_id = '%s'" % (cuffcomp_tables["exons"], transcript_data["transcript_id"]))
        num_exons = db_conn.rowcount
        exon_result = db_conn.fetchall()

        # defaults
        exon_starts = []
        exon_ends = []
        transcript_length = 0
        u_exon = "no"

        for exon_data in exon_result:
            # set the transcript strand to the exon strand
            # in theory every exon in the transcript should have the same strand
            if "strand" not in transcript_data:
                transcript_data["strand"] = exon_data["strand"]

            # is the class_code "u"?  if so, flag as having a "u" exon
            if exon_data["class_code"] == "u":
                u_exon = "yes"

            # for finding the pre-splice length
            exon_starts.append (exon_data["start"])
            exon_ends.append (exon_data["end"])

            transcript_length += exon_data["end"] - exon_data["start"] + 1;

        presplice_length = max(exon_ends) - min(exon_starts) + 1
        
        # use manual annotation if there's no cufflinks annotation
        if transcript_data["nearest_ref"] == None:
            transcript_data["nearest_ref"] = transcript_data["manual_nearest_ref"]
            
        if transcript_data["gene_name"] == None:
            transcript_data["gene_name"] = transcript_data["manual_gene_name"]

        # does the transcript have a reference annotation?
        if transcript_data["nearest_ref"] != None:
            # use the annotation information that corresponds to the nearest_ref ID
            # ENS = ensembl, uc = UCSC, mmu = mir, NM_/NR_ = refGene
            ref_prefix = {"ENS": "ensembl",
                          "uc": "ucsc",
                          "mmu": "mir",
                          "NM_": "refgene",
                          "NR_": "refgene"}

            # default, use all dbs
            use_this_db = "all"

            for k, v in ref_prefix.items():
                if transcript_data["nearest_ref"][:len (k)] == k:
                    use_this_db = v
                    break

            if use_this_db == "all":
                anno_info = annodb_obj.query_all_db (annodb_obj.QT_NAME_COLL, name = transcript_data["nearest_ref"])
            else:
                anno_info = annodb_obj.query_single_db (annodb_obj.QT_NAME_COLL, use_this_db, name = transcript_data["nearest_ref"])

            # maybe the annotation isn't in our database (shouldn't happen!)
            if len (anno_info) == 0:
                antisense_to_ref = "not_in_annodb"
                known_noncoding = "not_in_annodb"
            else:
                # does the transcript have a strand?
                if transcript_data["strand"] == ".":
                    antisense_to_ref = "strand_unknown"
                else:
                    # compare annotation strand to that of the transcript
                    last_comp = ""
                    antisense_to_ref = ""

                    for rec in anno_info:
                        if transcript_data["strand"] == rec["strand"]:
                            antisense_to_ref = "no"
                        else:
                            antisense_to_ref = "yes"

                        if last_comp != "" and last_comp != antisense_to_ref:
                            antisense_to_ref = "ambiguous"
                            break

                        last_comp = antisense_to_ref

                # is this known to be noncoding?
                last_nc_result = ""
                known_noncoding = ""

                for rec in anno_info:
                    try:
                        if rec["cdsEnd"] - rec["cdsStart"] == 0:
                            known_noncoding = "yes"
                        else:
                            known_noncoding = "no"
                            break
                    except KeyError:
                        known_noncoding = "unknown"

        else:
            antisense_to_ref = "unannotated"
            known_noncoding = "unannotated"

        build_query = {}

        build_query["num_exons"] = num_exons
        build_query["transcript_length"] = transcript_length
        build_query["presplice_length"] = presplice_length
        build_query["antisense"] = antisense_to_ref
        build_query["u_exon"] = u_exon
        build_query["known_noncoding"] = known_noncoding

        query = functions.mysql_build_query ("UPDATE %s SET " % (cuffcomp_tables["tx"],),
                                             build_query,
                                             " WHERE transcript_id = '%s'" % (transcript_data["transcript_id"],))

        db_conn.execute (query)

        transcripts_processed += 1

    functions.write_status_line ("Processing complete!", 100, True)
    print

def analysis_longest_orf (assembly, database, force = "no"):
    """
    determines the longest orf length for all spliced transcripts
    """

    import functions, bio_functions

    db_conn = functions.mysql_db_connect (database)

    if force == "no":
        dont_overwrite = " WHERE longest_orf_length IS NULL"
    elif force == "yes":
        dont_overwrite = ""

    db_conn.execute ("SELECT gene_name, transcript_id FROM %s%s" % (cuffcomp_tables["tx"], dont_overwrite))

    total_num_transcripts = db_conn.rowcount
    transcript_result = db_conn.fetchall()

    transcripts_processed = 0

    print "Finding longest ORFs..."

    for transcript_data in transcript_result:
        if transcripts_processed % 50 == 0:
            functions.write_status_line ("Processing %s" % (transcript_data["transcript_id"],),
                                         int(100 * (float(transcripts_processed) / float(total_num_transcripts))))

        transcript_seq = _get_transcript_seq (db_conn, assembly, transcript_data["transcript_id"], cuffcomp_tables["exons"], "", False)

        orfs = bio_functions.find_orfs (transcript_seq)

        longest_orf = 0
        longest_start = 0
        longest_end = 0

        for (length, start, end, strand, frame) in orfs:
            if length > longest_orf:
                longest_orf = length
                longest_start = start
                longest_end = end

        build_query = {}

        build_query["longest_orf_length"] = longest_orf
        build_query["longest_orf_coords"] = "%d,%d" % (longest_start, longest_end)

        query = functions.mysql_build_query ("UPDATE %s SET " % (cuffcomp_tables["tx"],),
                                             build_query,
                                             " WHERE transcript_id = '%s'" % (transcript_data["transcript_id"],))

        db_conn.execute (query)

        transcripts_processed += 1

    functions.write_status_line ("Processing complete!", 100, True)
    print

def analysis_rfc_score (assembly, database, chrom = None):
    """
    calculates reading frame conservation score for all transcripts
    """

    import functions, bio_functions, rfc
    from statlib import stats

    rfc = rfc.RFCScore (assembly)

    db_conn = functions.mysql_db_connect (database)

    if chrom == None:
        db_conn.execute ("SELECT transcript_id, longest_orf_coords FROM %s WHERE rfc_score IS NULL" % (cuffcomp_tables["tx"],))
    else:
        db_conn.execute ("SELECT cuffcomp_transcripts.transcript_id, cuffcomp_transcripts.longest_orf_coords FROM %s, cuffcomp_exons WHERE rfc_score IS NULL AND cuffcomp_transcripts.transcript_id = cuffcomp_exons.transcript_id AND cuffcomp_exons.chr = '%s' GROUP BY cuffcomp_transcripts.transcript_id" % (cuffcomp_tables["tx"], chrom))

    total_num_transcripts = db_conn.rowcount
    transcript_result = db_conn.fetchall()

    transcripts_processed = 0

    print "Calculating RFC score of longest ORF in all transcripts..."

    for transcript_data in transcript_result:
        functions.write_status_line ("Processing %s" % (transcript_data["transcript_id"],),
                                     int(100 * (float(transcripts_processed) / float(total_num_transcripts))))

        db_conn.execute ("SELECT chr, start, end FROM %s WHERE transcript_id = '%s'" % (cuffcomp_tables["exons"], transcript_data["transcript_id"]))

        exon_rows = db_conn.fetchall()

        # make exon schema for this transcript
        all_exons = []

        for exon_data in exon_rows:
            all_exons.append ((exon_data["chr"], int(exon_data["start"]), int(exon_data["end"]) + 1))

        orf_coords = transcript_data["longest_orf_coords"].split(",")

        longest_orf_score = rfc.get_longest_orf_score (all_exons, (int(orf_coords[0]), int(orf_coords[1])))

        build_query = {}
        build_query["transcript_id"] = transcript_data["transcript_id"]

        for asm, val in longest_orf_score.iteritems():
            build_query[asm] = val

        query = functions.mysql_build_query ("INSERT INTO rfc_longest_orf SET ", build_query)
        db_conn.execute (query)
        
        db_conn.execute ("UPDATE cuffcomp_transcripts SET rfc_score = '0' WHERE transcript_id = '%s'" % (transcript_data["transcript_id"],))

        transcripts_processed += 1

    functions.write_status_line ("Processing complete!", 100, True)
    print
    
    print "Finding 95th-percentile for complete, known protein-coding transcripts..."
    
    db_conn.execute ("SHOW FIELDS IN rfc_longest_orf")
    
    total_assemblies = db_conn.rowcount
    assemblies_processed = 0
    
    assemblies_result = db_conn.fetchall()
    
    asm_95p = {}
    asm_track = []
    
    for assembly_data in assemblies_result:
        asm = assembly_data["Field"]
         
        if asm == assembly or asm == "transcript_id": continue
        
        if not asm in asm_track:
            asm_track.append (asm)
        
        functions.write_status_line ("Processing %s" % (asm,),
                                     int(100 * (float(assemblies_processed) / float(total_assemblies))))

        db_conn.execute ("SELECT IF(%s IS NULL, -1, %s) AS score FROM rfc_longest_orf, cuffcomp_transcripts WHERE cuffcomp_transcripts.transcript_id = rfc_longest_orf.transcript_id AND cuffcomp_transcripts.known_noncoding = 'no' AND cuffcomp_transcripts.class_code = '='" % (asm, asm))

        score_result = db_conn.fetchall()
        
        scores = []        
        
        for score_data in score_result:
            if float (score_data["score"]) > -1:
                scores.append (float (score_data["score"]))

        asm_95p[asm] = stats.scoreatpercentile (scores, 0.05)
        
        assemblies_processed += 1
        
    functions.write_status_line ("Processing complete!", 100, True)
    print

    if chrom == None:
        db_conn.execute ("SELECT rfc_longest_orf.* FROM rfc_longest_orf") 
    else:
        db_conn.execute ("SELECT rfc_longest_orf.* FROM rfc_longest_orf, cuffcomp_exons WHERE rfc_longest_orf.transcript_id = cuffcomp_exons.transcript_id AND cuffcomp_exons.chr = '%s' GROUP BY rfc_longest_orf.transcript_id" % (chrom,))

    total_num_transcripts = db_conn.rowcount
    transcript_result = db_conn.fetchall()

    transcripts_processed = 0

    print "Tallying votes for protein-coding potential by organism..."

    for transcript_data in transcript_result:
        if transcripts_processed % 100 == 0:
            functions.write_status_line ("Processing %s" % (transcript_data["transcript_id"],),
                                         int(100 * (float(transcripts_processed) / float(total_num_transcripts))))
                                     
        tx_protein_score = 0
        voters = 0
                                     
        for asm in asm_track:
            if transcript_data[asm] == None: continue
            if float (transcript_data[asm]) < 0: continue
        
            if transcript_data[asm] >= asm_95p[asm]:
                tx_protein_score += 1
            else:
                tx_protein_score -= 1
                
            voters += 1
                
        db_conn.execute ("UPDATE cuffcomp_transcripts SET rfc_score = '%d', rfc_voters = '%d', rfc_norm = IF(rfc_score = 0, 0, rfc_score / rfc_voters) WHERE transcript_id = '%s'" % (tx_protein_score, voters, transcript_data["transcript_id"]))
                
        transcripts_processed += 1
        
    functions.write_status_line ("Processing complete!", 100, True)
    print
    
def analysis_csf_score (assembly, database, lookup, placental = True, forked = False):    
    import functions, bio_functions, csf, os, time
    
    if lookup == "fork":
        print "Calculating CSF Score..."
        
        functions.__LOG_OUTPUT = True
        
        forked_pids = []        
        
        for chrom in bio_functions.valid_chr[assembly]:
            fork_pid = os.fork()
            
            if fork_pid == 0:
                analysis_csf_score (assembly, database, chrom, placental, True)
                
                sys.exit (0)
            else:
                forked_pids.append (fork_pid)
        
        while len (forked_pids) > 0:
            for pid in forked_pids:
                waitpid_result = os.waitpid (pid, os.WNOHANG)

                if __debug__: debug_print ("analysis_csf_score waitpid_result (%d)=%s" % (pid, waitpid_result), 6)

                if waitpid_result[0] == pid:
                    if __debug__: debug_print ("analysis_csf_score finished, removing %d" % (pid,), 5)
                    forked_pids.remove (pid)

            time.sleep (0.1)
        
        sys.exit (0)

    csf_obj = csf.CSFScore (assembly)
    csf_obj.load_matrices()

    if placental == "True" or placental == True:
        csf_field = "csf_score_20placental"
    elif placental == "False" or placental == False:
        csf_obj.analysis_assemblies = []
        csf_field = "csf_score"
    else:
        raise SystemError ("csf::get_tx_scores placental must be true or false")

    db_conn = functions.mysql_db_connect (database)

    if lookup not in bio_functions.valid_chr[assembly]:
        db_conn.execute ("SELECT * FROM %s WHERE transcript_id = '%s'" % (cuffcomp_tables["tx"], lookup))
    else:
        db_conn.execute ("SELECT cuffcomp_transcripts.* FROM cuffcomp_transcripts, cuffcomp_exons WHERE cuffcomp_transcripts.transcript_id = cuffcomp_exons.transcript_id AND cuffcomp_exons.chr = '%s' AND cuffcomp_transcripts.csf_score IS NULL GROUP BY cuffcomp_transcripts.transcript_id ORDER BY cuffcomp_exons.start ASC" % (lookup,))

    transcript_result = db_conn.fetchall()
    total_num_transcripts = db_conn.rowcount

    transcripts_processed = 0
    
    if forked == False: print "Calculating CSF Score..."

    for transcript_data in transcript_result:
        if forked == True:
            functions.write_status_line ("[%s]\tProcessing %s" % (lookup, transcript_data["transcript_id"]),
                                         int(100 * (float(transcripts_processed) / float(total_num_transcripts))))
        else:
            functions.write_status_line ("Processing %s" % (transcript_data["transcript_id"],),
                                         int(100 * (float(transcripts_processed) / float(total_num_transcripts))))

        if __debug__: print

        db_conn.execute ("SELECT chr, start, end, strand FROM cuffcomp_exons WHERE transcript_id = '%s' ORDER BY exon_num ASC" % (transcript_data["transcript_id"],))

        exon_rows = db_conn.fetchall()

        # make exon schema for this transcript
        all_exons = []

        for exon_data in exon_rows:
            all_exons.append ((exon_data["chr"], int(exon_data["start"]), int(exon_data["end"]) + 1))

        if exon_data["strand"] == "-": antisense = True
        else: antisense = False

        csf_score = csf_obj.sliding_window (all_exons, antisense, None)

        db_conn.execute ("UPDATE cuffcomp_transcripts SET %s = '%s' WHERE transcript_id = '%s'" % (csf_field, csf_score, transcript_data["transcript_id"]))

        transcripts_processed += 1

    if forked == True:
        functions.write_status_line ("[%s]\tProcessing complete!" % (lookup,), 100, True)
    else:
        functions.write_status_line ("Processing complete!", 100, True)
        
    if forked == False: print

def analysis_phastcons_score (assembly, database):
    """
    calculates the phastcons count, z-score, and p-value for all transcripts
    """

    import functions, bio_functions, phastcons
    from statlib import stats

    db_conn = functions.mysql_db_connect (database)
#    db_conn.execute ("SELECT transcript_id FROM %s WHERE phastcons_zscore IS NULL" % (cuffcomp_tables["tx"]))
    db_conn.execute ("SELECT transcript_id FROM %s WHERE phastcons_zscore IS NULL AND nearest_ref IS NULL AND manual_nearest_ref IS NULL" % (cuffcomp_tables["tx"]))

    total_num_transcripts = db_conn.rowcount
    transcript_result = db_conn.fetchall()

    transcripts_processed = 0

    print "Calculating phastConsElements counts and z-scores..."

    for transcript_data in transcript_result:
        functions.write_status_line ("Processing %s" % (transcript_data["transcript_id"],),
                                     int(100 * (float(transcripts_processed) / float(total_num_transcripts))))

        db_conn.execute ("SELECT chr, start, end FROM %s WHERE transcript_id = '%s'" % (cuffcomp_tables["exons"], transcript_data["transcript_id"]))

        exon_rows = db_conn.fetchall()

        all_exons = []
        
        for exon_data in exon_rows:
            all_exons.append ((exon_data["chr"], exon_data["start"], exon_data["end"] + 1))
            
        (phastcons_count, phastcons_zscore, phastcons_pvalue) = phastcons.get_zscore (assembly, all_exons)

        db_conn.execute ("UPDATE %s SET phastcons_count = '%d', phastcons_zscore = '%f', phastcons_pvalue = '%f' WHERE transcript_id = '%s';" % (cuffcomp_tables["tx"], phastcons_count, phastcons_zscore, phastcons_pvalue, transcript_data["transcript_id"]))

        transcripts_processed += 1

    functions.write_status_line ("Processing complete!", 100, True)
    print

def analysis_siphy_score (assembly, database):
    """
    calculates the max 12bp sliding window score in siphy for all exons
    """

    import functions, bio_functions, siphy

    db_conn = functions.mysql_db_connect (database)
    db_conn.execute ("SELECT id, chr, start, end FROM %s WHERE siphy_score IS NULL" % (cuffcomp_tables["exons"]))

    total_num_exons = db_conn.rowcount
    exon_result = db_conn.fetchall()

    exons_processed = 0

    print "Calculating max 12 bp SiPhy score for all exons..."

    siphy_score = 0

    for exon_data in exon_result:
        functions.write_status_line ("Processing exon %s (last score=%d)" % (exon_data["id"], siphy_score),
                                     int(100 * (float(exons_processed) / float(total_num_exons))))

        siphy_score = siphy.max_sliding_window (assembly, exon_data["chr"], exon_data["start"], exon_data["end"] + 1, 12)

        db_conn.execute ("UPDATE %s SET siphy_score = '%f' WHERE id = '%s'" % (cuffcomp_tables["exons"], siphy_score, exon_data["id"]))

        exons_processed += 1

    functions.write_status_line ("Processing complete!", 100, True)
    print

def _usage (help_cmd = None):
    if help_cmd == "import_cuffcomp":
        print \
"""Usage: python %(cmd)s import_cuff_comp <database> <directory> <prefix> [subtask]

 Inserts cuffcompare output into three database tables: cuffcomp_loci,
 cuffcomp_transcripts, and cuffcomp_exons.  Requires a valid path to
 cuffcompare output files (<prefix>.combined.gtf, <prefix>.loci, <prefix>.tracking)
 and a valid database name.

Options:
[subtask]  direct the script to perform one or more parts of the overall task
           by specifying a comma-delimited list of one or more of the following:
               genes
               exons
               transcripts
""" % \
{"cmd": sys.argv[0]}

    if help_cmd == "analysis_first_pass":
        print \
"""Usage: python %(cmd)s analysis_first_pass <assembly> <database> <tx_table> <exon_table>

 Counts the number of exons in each transcript, notes the spliced/pre-splice transcript
 length, if a transcript has any class_code "u" exons, and if it is antisense
 to the reference annotation in RefGene.

Options:
assembly  the version of the genome assembly in use
          (ex: mm9 for mouse, hg18 for human)
""" % \
{"cmd": sys.argv[0]}

    if help_cmd == "analysis_phastcons_score":
        print \
"""Usage: python %(cmd)s analysis_phastcons_score <assembly> <database> <tx_table> <exon_table>

 Counts the number of nucleotides in a spliced transcript that appear in the
 phastConsElements track from UCSC.  Takes a random sample of size- and schema-
 matched transcripts from the given assembly, and determines a z-score and
 p-value for the transcript based on the random distribution.

Options:
assembly  the version of the genome assembly in use
          (ex: mm9 for mouse, hg18 for human)
""" % \
{"cmd": sys.argv[0]}

    if help_cmd == "analysis_siphy_score":
        print \
"""Usage: python %(cmd)s analysis_siphy_score <assembly> <database> <exon_table>

 Slides 12bp windows along each exon individually and calculates the maximum score
 found.

Options:
assembly  the version of the genome assembly in use
          (ex: mm9 for mouse, hg18 for human)
""" % \
{"cmd": sys.argv[0]}

    if help_cmd == None:
        print \
"""Usage: python %(cmd)s <task> [options]
   Or: python %(cmd)s --help <task>

Tasks: %(tasks)s
""" % \
{"cmd": sys.argv[0], "tasks": "\n       ".join (valid_tasks)}

    sys.exit(0)

if __name__ == "__main__":
    valid_tasks = ("full_run",
                   "import_cuffcomp",
                   "manual_annotation",
                   "nearest_annotation",
                   "analysis_first_pass",
                   "analysis_longest_orf",
                   "analysis_phastcons_score",
                   "analysis_siphy_score",
                   "analysis_csf_score",
                   "analysis_rfc_score",
                   "get_transcript_sequence",
                   "create_comp_tables")

    # if ran with no parameters, print usage
    try:
        sys.argv[1]
    except IndexError:
        print _usage()

    # redirect to the appropriate help usage
    if sys.argv[1] == "--help":
        try:
            sys.argv[2]
        except IndexError:
            print _usage()

        _usage (sys.argv[2])

    # if the task is valid, pass the arguments to the function
    if sys.argv[1] in valid_tasks:
        globals()[sys.argv[1]] (*sys.argv[2:])
    else:
        print _usage()