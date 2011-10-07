# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "10/7/2011"
__version__ = 0.0

import os
import subprocess
import _common
import time
import urllib
import gzip
import cStringIO

def _ucsc_path(assembly, table):
    return _common.UCSC_HTTP_PATH % {"assembly": assembly,
                             "table": table}
                             
def _http_gzip_fetch(url):
    url = urllib.urlopen(url)    
    url_raw = cStringIO.StringIO(url.read())
    url_gzip = gzip.GzipFile(fileobj = url_raw)
    
    try:
        url_gzip.readline()
        url_gzip.seek(0)
    except IOError:
        url_gzip = []
    
    return list(url_gzip)
    
def fetch_ucsc_table(assembly, table):
    return _http_gzip_fetch(_ucsc_path(assembly, table))
    
def fetch_ucsc_gtf(assembly, table):
    table = fetch_ucsc_table(assembly, table)
    
    # remove last two columns, since they fuck it up
    table_stripped = []

    for line in table:
        line_split = line.strip().split("\t")
        
        table_stripped.append("\t".join(line_split[:-2]).strip())

    table_stripped = "\n".join(table_stripped)
    
    # check if it's valid
    check_proc = subprocess.Popen([os.path.join(_common.UCSC_KENT_PATH, "genePredCheck"),
                                   "stdin"],
                                  stdin=subprocess.PIPE,
                                  stderr=subprocess.PIPE)    
    
    check_proc.stdin.write(table_stripped)
    check_proc.stdin.close()
    
    while check_proc.poll() == None:
        time.sleep(0.1)
        
    result = check_proc.stderr.readline().strip()
    
    if not result.endswith("failed: 0"):
        raise ValueError("Assembly: %s Table: %s failed genePredCheck" % (assembly,
                                                                          table))

    # convert it to a GTF
    convert_proc = subprocess.Popen([os.path.join(_common.UCSC_KENT_PATH, "genePredToGtf"),
                                     "file",
                                     "stdin",
                                     "stdout"],
                                     stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)

    convert_proc.stdin.write(table_stripped)
    convert_proc.stdin.close()
    
    output = []
    
    while convert_proc.poll() == None:
        data = convert_proc.stdout.read(4096)
        
        if line:
            output.append(data)
        
    result = convert_proc.stderr.readline().strip()
    
    if result != "":
        raise ValueError("Assembly: %s Table: %s failed genePredToGtf" % (assembly,
                                                                          table))

    return "".join(output).split("\n")