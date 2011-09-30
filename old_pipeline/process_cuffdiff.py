#                 {type: (table name, key name, comment)}
cuffdiff_tables = {"tx": ("cuffdiff_transcripts",
                          "transcript_id",
                          "FPKM by transcript_id"),
                   "loci": ("cuffdiff_loci",
                            "locus_id",
                            "FPKM by locus_id"),
                   "genes": ("cuffdiff_gene_name",
                             "gene_name",
                             "Summed FPKM of ''='' or ''c'' transcripts that share a gene name")}

cuffdiff_views = {"known_noncoding": "SELECT cuffcomp_loci.locus_id, cuffcomp_transcripts.nearest_ref, cuffcomp_transcripts.gene_name, cuffcomp_loci.chr, cuffcomp_loci.start, cuffcomp_loci.end, cuffcomp_loci.strand, %s FROM cuffcomp_transcripts, cuffdiff_loci, cuffcomp_loci WHERE cuffdiff_loci.locus_id = cuffcomp_transcripts.locus_id AND cuffdiff_loci.locus_id = cuffcomp_loci.locus_id AND known_noncoding = 'yes' GROUP BY cuffcomp_loci.locus_id",
                  "unanno_lt300orf_gt1exon": "SELECT cuffcomp_transcripts.transcript_id, cuffcomp_transcripts.locus_id, cuffcomp_loci.chr, cuffcomp_loci.start, cuffcomp_loci.end, cuffcomp_transcripts.transcript_length, cuffcomp_transcripts.presplice_length, cuffcomp_transcripts.longest_orf_length, cuffcomp_transcripts.phastcons_zscore, %s FROM cuffcomp_transcripts, cuffcomp_loci, cuffdiff_loci WHERE cuffcomp_transcripts.locus_id = cuffdiff_loci.locus_id AND cuffcomp_loci.locus_id = cuffcomp_transcripts.locus_id AND cuffcomp_transcripts.num_exons > 1 AND cuffcomp_transcripts.nearest_ref IS NULL AND cuffcomp_transcripts.longest_orf_length <= 300"}

__DEBUG_LEVEL = 5

if __debug__:
    def debug_print (msg, msg_level):
        if __DEBUG_LEVEL >= msg_level:
            print "debug%d\t%s" % (msg_level, msg)
            
import sys

sys.path.append ("/comp_node0/andrew/work/biolibrary")

def create_diff_tables (database, field_name_list):
    import functions

    db_conn = functions.mysql_db_connect (database)

    field_names = field_name_list.split(",")

    if len (field_names) < 2:
        raise SystemError ("process_cuffdiff::create_diff_tables (%s, %s) invalid field name list" % \
                           (database, field_name_list))

    for table_type, table_info in cuffdiff_tables.items():
        fold_change_fields = []
        final_fields = ["`%s` char(40) NOT NULL" % (table_info[1],)]

        for k, v in enumerate (field_names):
            final_fields.append ("`%s` float NOT NULL" % (v,))

            for k1, v1 in enumerate (field_names[k + 1:]):
                fold_change_fields.append ("`%s_%s` float NULL" % (v, v1))

        final_fields.extend (fold_change_fields)
        final_fields.append ("PRIMARY KEY (`%s`)" % (table_info[1],))

        query = "CREATE TABLE `%s` (%s) ENGINE=MyISAM DEFAULT CHARSET=latin1 COMMENT='%s';" % (table_info[0], ", ".join (final_fields), table_info[2])

        db_conn.execute (query)

def create_views (database, field_name_list):
    import functions

    db_conn = functions.mysql_db_connect (database)

    field_names = field_name_list.split(",")

    loci_diff_fields = []
    fold_change_fields = []

    for k, v in enumerate (field_names):
        loci_diff_fields.append ("cuffdiff_loci.%s" % (v,))

        for k1, v1 in enumerate (field_names[k + 1:]):
            fold_change_fields.append ("cuffdiff_loci.%s_%s" % (v, v1))

    loci_diff_fields.extend (fold_change_fields)

    for view_name, view_query in cuffdiff_views.items():
        select_statement = view_query % (", ".join (loci_diff_fields))

        query = "CREATE VIEW %s AS %s" % (view_name, select_statement)

        db_conn.execute (query)

def _check_num_rows (db_conn, table, field, count):
    # check that number of rows in table = number of queries generated
    db_conn.execute ("SELECT COUNT(%s) FROM %s" % (field, table))
    rows_in_table = db_conn.fetchone()["COUNT(%s)" % (field,)]

    if count != rows_in_table:
        print "WARNING: Number of rows does not match number of queries!"
        print "\trows=%d lines=%d" % (rows_in_table, count)
    else:
        print "Check OK -- number of rows matches number of queries"

def _fold_change_query (table_name, field_names):
    import functions

    build_query = {}

    for k, v in enumerate (field_names):
        for k1, v1 in enumerate (field_names[k + 1:]):
            build_query["%s_%s" % (v, v1)] = "IF (%(num)s = 0 AND %(denom)s = 0, 0, IF (%(num)s > %(denom)s, %(num)s, -%(denom)s) / IF (%(num)s > %(denom)s, IF(%(denom)s < 1, 1, %(denom)s), IF(%(num)s < 1, 1, %(num)s)))" % {"num": v, "denom": v1}

    return functions.mysql_build_query ("UPDATE %s SET " % (table_name,), build_query, "", "")

def import_cuffdiff_nonovel (database, directory, desired_fields, field_names):
    def _multi_query_insert (table, fields, data):
        return "INSERT INTO %s (%s) VALUES %s;" % (table, ", ".join (fields), ", ".join (data))
        
    import os, functions
    
    # input sanitization
    directory = os.path.realpath (directory)
    
    desired_fields = desired_fields.split (",")
    field_names = field_names.split (",")

    db_conn = functions.mysql_db_connect (database)    
    
    with open ("tmp2", "r") as fp:
        for line in fp:
            db_conn.execute ("SELECT * FROM genes_fpkm WHERE gene_name = '%s'" % (line.strip(),))
            
            result = db_conn.fetchall()
            
            for row in result:
                print row
                
    sys.exit(0)
        
    
    # import fpkm by transcripts
    with open (directory + "/isoforms.fpkm_tracking", "r") as fp:
        # read the first line, figure out which fields are where
        line_split = fp.readline().split()
        
        field_to_fieldnum = []
        
        for filefield in desired_fields:
            for fieldnum, fieldval in enumerate (line_split):
                if filefield == fieldval:
                    field_to_fieldnum.append (fieldnum)

        assert len (field_to_fieldnum) == len (field_names) == len (desired_fields)
        
        nice_list = zip (field_names, field_to_fieldnum)
        
        gene_totals = {}
        
        # process the rest of the lines
        tmp_inserts = []
        
        for line in fp:
            line_split = line.split()
            
            # if gene name is blank, use the transcript name
            if line_split[3] == "-": line_split[3] = line_split[0]
            
            try:
                gene_totals[line_split[3]]
            except KeyError:
                gene_totals[line_split[3]] = {}
            
            this_insert = [line_split[0], line_split[3]]
            
            for sqlfield, fieldnum in nice_list:                
                try:
                    gene_totals[line_split[3]][sqlfield]
                except KeyError:
                    gene_totals[line_split[3]][sqlfield] = 0
                    
                gene_totals[line_split[3]][sqlfield] += float (line_split[fieldnum])
                
                this_insert.append (line_split[fieldnum])
                
            tmp_inserts.append ("('%s')" % ("', '".join (this_insert),))
            
            if len (tmp_inserts) % 1000 == 0:
                db_conn.execute (_multi_query_insert ("isoforms_fpkm", ["isoform_name", "gene_name"] + field_names, tmp_inserts))
                
                tmp_inserts = []
                
        if len (tmp_inserts) > 0:
            db_conn.execute (_multi_query_insert ("isoforms_fpkm", ["isoform_name", "gene_name"] + field_names, tmp_inserts))

    # import fpkm by genes, summing by gene name
    tmp_inserts = []
    
    for gene in gene_totals.iterkeys():
        this_insert = [gene]
        
        for sqlfield, _ in nice_list:
            this_insert.append (str (gene_totals[gene][sqlfield]))
            
        tmp_inserts.append ("('%s')" % ("', '".join (this_insert),))
        
        if len (tmp_inserts) % 1000 == 0:
            db_conn.execute (_multi_query_insert ("genes_fpkm", ["gene_name"] + field_names, tmp_inserts))
            
            tmp_inserts = []
            
    if len (tmp_inserts) > 0:
        db_conn.execute (_multi_query_insert ("genes_fpkm", ["gene_name"] + field_names, tmp_inserts))
        
    db_conn.execute (_fold_change_query ("genes_fpkm", field_names))
    db_conn.execute (_fold_change_query ("isoforms_fpkm", field_names))

def import_cuffdiff (database, directory, field_name_list, subtask = None):
    """
    a wrapper function to call various import tasks
    """

    import os

    valid_subtasks = {"loci": ("genes.fpkm_tracking", cuffdiff_tables["loci"][0], cuffdiff_tables["loci"][1]),
                      "transcripts": ("isoforms.fpkm_tracking", cuffdiff_tables["tx"][0], cuffdiff_tables["tx"][1])}

    # input sanitization
    directory = os.path.realpath (directory)

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
            raise SystemError ("process_cuffdiff::import_cuffdiff (%s, %s, %s, %s) invalid subtask" % \
                               (database, directory, field_name_list, subtask))

        filenames[i] = "%s/%s" % (directory, valid_subtasks[i][0])

        if not os.path.exists (filenames[i]):
            raise SystemError ("process_cuffdiff::import_cuffdiff (%s, %s, %s, %s) %s does not exist" % \
                               (database, directory, field_name_list, subtask, filenames[i]))

    field_names = field_name_list.split(",")

    if len (field_names) < 2:
        raise SystemError ("process_cuffdiff::import_cuffdiff (%s, %s, %s, %s) invalid field name list" % \
                           (database, directory, field_name_list, subtask))

    # do it to it
    for i in tasks_to_do:
        _import_tracking_file (filenames[i], database, valid_subtasks[i][1], [valid_subtasks[i][2]] + field_names)
        print

def _import_tracking_file (filename, database, table_name, field_names):
    """
    import a cuffdiff tracking file into database.table_name
    """

    import functions

    lines_in_file = functions.count_lines_in_file (filename)

    cuffdiff_tracking = open(filename, "r")

    db_conn = functions.mysql_db_connect (database)

    queries_generated = 0

    print "Importing data into %s.%s..." % (database, table_name)

    # skip the first (header) line
    for line in cuffdiff_tracking.readlines()[1:]:
        # should be tab-delimited: split it into a list
        data = line.split("\t")

        # output the status bar, style stolen from cufflinks
        if queries_generated % 100 == 0:
            functions.write_status_line("Processing %s" % (data[0],),
                                        int(100 * (float(queries_generated) / float(lines_in_file))))

        build_query = {}

        # populate build_query with select fields from this line
        build_query[field_names[0]] = data[0]

        # the fpkm fields start at 6 and are every 3rd field
        for k, v in enumerate (field_names[1:]):
            build_query[v] = data[6 + (3 * k)]

        query = functions.mysql_build_query ("INSERT INTO %s SET " % (table_name,), build_query)

        db_conn.execute (query)

        queries_generated += 1

    functions.write_status_line ("Processing complete!", 100, True)

    # check if insert was OK
    _check_num_rows (db_conn, table_name, field_names[0], queries_generated)

    # add fold change information
    db_conn.execute (_fold_change_query (table_name, field_names[1:]))

    # cleanup
    db_conn.close()
    cuffdiff_tracking.close()

def populate_gene_table (database, field_name_list):
    """
    wrapper to call _populate_gene_table from the command line
    """

    import process_cuffcompare

    field_names = field_name_list.split(",")

    if len (field_names) < 2:
        raise SystemError ("process_cuffdiff::populate_gene_table (%s, %s) invalid field name list" % \
                           (database, field_name_list))

    # do it to it
    _populate_gene_table (database, cuffdiff_tables["tx"][0], process_cuffcompare.cuffcomp_tables["tx"], cuffdiff_tables["genes"][0], field_names)
    print

def _populate_gene_table (database, diff_tx_table, comp_tx_table, new_table, field_names):
    """
    calculates fpkm by annotated gene name from cuffcomp/cuffdiff analyses
    """

    import functions

    db_conn = functions.mysql_db_connect (database)

    # get all possible not-null gene_names from cuffcomp_transcripts
    db_conn.execute ("SELECT gene_name FROM %s WHERE gene_name IS NOT NULL GROUP BY gene_name" % (comp_tx_table,))

    all_genes = db_conn.fetchall()

    gene_count = db_conn.rowcount

    # build list of fpkm field names
    select_fields = ", ".join(field_names)

    queries_generated = 0

    print "Populating %s.%s..." % (database, new_table)

    for gene_row in all_genes:
        # output the status bar, style stolen from cufflinks
        if queries_generated % 100 == 0:
            functions.write_status_line("Processing %s" % (gene_row["gene_name"],),
                                        int(100 * (float(queries_generated) / float(gene_count))))

        # get fpkm values for all transcripts with this gene_name
        db_conn.execute ("SELECT %(select_fields)s FROM %(comp_table)s, %(diff_table)s WHERE %(comp_table)s.transcript_id = %(diff_table)s.transcript_id AND %(comp_table)s.gene_name = '%(gene)s' AND (%(comp_table)s.class_code = '=' OR %(comp_table)s.class_code = 'c')" % \
                        {"select_fields": select_fields,
                         "comp_table": comp_tx_table,
                         "diff_table": diff_tx_table,
                         "gene": gene_row["gene_name"]})

        if db_conn.rowcount > 0:
            all_gene_fpkm = db_conn.fetchall()

            gene_sum = {}

            for fpkm in all_gene_fpkm:
                for field in field_names:
                    try:
                        gene_sum[field] += fpkm[field]
                    except KeyError:
                        gene_sum[field] = fpkm[field]

            build_query = {}

            build_query["gene_name"] = gene_row["gene_name"]

            for field in field_names:
                try:
                    build_query[field] = gene_sum[field]
                except KeyError:
                    build_query[field] = 0

            query = functions.mysql_build_query ("INSERT INTO %s SET " % (new_table,), build_query)

            db_conn.execute (query)

            queries_generated += 1

    functions.write_status_line ("Processing complete!", 100, True)

    # check if insert was OK
    _check_num_rows (db_conn, new_table, field_names[0], queries_generated)

    # add fold change information
    db_conn.execute (_fold_change_query (new_table, field_names))

    # cleanup
    db_conn.close()

def _usage (help_cmd = None):
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
                   "import_cuffdiff",
                   "import_cuffdiff_nonovel",
                   "populate_gene_table",
                   "create_views",
                   "create_diff_tables")

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