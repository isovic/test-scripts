#!/usr/bin/env python

"""Convert vcf entry, i.e. assume ref reads were mapped against alt"""

import gzip
import sys
import os
for p in ['/mnt/pnsg10_home/wilma/local/src/compbio-utils-git/', '/home/wilma/local/src/compbio-utils-git/']:
    if os.path.exists(p):
        sys.path.insert(0, p)
import vcf_read_support

def simple_vcf_write(var):
    info = ";".join(["%s:%s" % (k,v) for (k,v) in var.info.items()])
    return "%s\t%s" % ("\t".join([var.chrom, "%d" % (var.pos+1), var.id, var.ref, var.alt, "%s" % (var.qual), var.filter]), info)

fh = gzip.open(sys.argv[1])
fho = sys.stdout

offset = 0
for var in vcf_read_support.simple_vcf_reader(fh):
    var = var._replace(ref=var.alt, alt=var.ref, pos=var.pos+offset)
    # for mutatrix swap TYPE:ins and TYPE:del
    fho.write(simple_vcf_write(var) + "\n")
    offset += (len(var.ref)-1)
    offset -= (len(var.alt)-1)
if fho != sys.stdout:
    fho.close()

sys.stderr.write("NOTE: doesn't reverse indel type if listed in INFO field\n")