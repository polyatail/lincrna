# -*- coding: utf-8 -*-

__author__ = "Andrew Sczesnak"
__email__ = "andrew.sczesnak@med.nyu.edu"
__date__ = "11/10/2011"
__version__ = 0.0

import _common
import urllib
import urllib2
from MultipartPostHandler import MultipartPostHandler as mpph

class UCSCImage():
    def __init__(self, assembly):
        self.sessid = self.get_cookie()
        self.clade, self.organism, self.db = _common.UCSC_PARAMS[assembly]
        
        self.custom_tracks = []
    
    def get_cookie(self):
        req = urllib2.Request(_common.UCSC_HGCUSTOM)
        response = urllib2.urlopen(req)
        cookie = response.headers.get("set-cookie")
    
        for item in cookie.split(";"):
            if item.startswith("hguid"):
                sessid = item[6:]
    
        return sessid

    def add_track(self, filename):
        header = open(filename, "r").readline()

        if header.startswith("track"):
            split_header = header.split(" ")

            for item in split_header:
                if item.startswith("name"):
                    trackname = item[5:].strip().strip("\"").strip("\'")
        
        fields = {"hgsid": self.sessid,
                  "clade": self.clade,
                  "org": self.organism,
                  "db": self.db,
                  "Submit": "submit",
                  "hgct_customText": "", 
                  "hgct_docText": ""} 
        files = {"hgt.customFile": open(filename, "r"),
                 "hgt.docFile": open("/dev/null", "rb")}
    
        payload = mpph.multipart_encode(fields.items(), files.items())
    
        req = urllib2.Request(_common.UCSC_HGCUSTOM, payload[1])
        req.add_header("content-type", "multipart/form-data; boundary=" + payload[0])
        req.add_header("cookie", "hguid=" + self.sessid)
        response = urllib2.urlopen(req).read()

        trackname_pos = response.find("ct_" + trackname)
        
        self.custom_tracks.append(response[trackname_pos:trackname_pos+len(trackname)+8])
        
    def fetch_image(self, position, pix = 800, params = {}):
        params["hgsid"] = self.sessid
        params["db"] = self.db
        params["position"] = position
        params["pix"] = pix
        
        for trackname in self.custom_tracks:
            params[trackname] = "full"
        
        query_string = urllib.urlencode(params.items())
        
        req = urllib2.Request(_common.UCSC_RENDERTRACKS + "?" + query_string)
        req.add_header("cookie", "hguid=" + self.sessid)
    
        response = urllib2.urlopen(req)
    
        return response.read()
        
