#!/usr/bin/env python

import pysam
import json
import argparse, sys, os.path
import math, time, re
from collections import Counter
from argparse import RawTextHelpFormatter
from functools import wraps

import svtyper.version

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
svtyper\n\
author: " + svtyper.version.__author__ + "\n\
version: " + svtyper.version.__version__ + "\n\
description: Compute genotype of structural variants based on breakpoint depth")
    parser.add_argument('-i', '--input_vcf', metavar='FILE', type=argparse.FileType('r'), default=None, help='VCF input (default: stdin)')
    parser.add_argument('-o', '--output_vcf', metavar='FILE',  type=argparse.FileType('w'), default=sys.stdout, help='output VCF to write (default: stdout)')
    parser.add_argument('-B', '--bam', metavar='FILE', type=str, required=True, help='BAM or CRAM file(s), comma-separated if genotyping multiple samples')
    parser.add_argument('-T', '--ref_fasta', metavar='FILE', type=str, required=False, default=None, help='Indexed reference FASTA file (recommended for reading CRAM files)')
    parser.add_argument('-S', '--split_bam', type=str, required=False, help=argparse.SUPPRESS)
    parser.add_argument('-l', '--lib_info', metavar='FILE', dest='lib_info_path', type=str, required=False, default=None, help='create/read JSON file of library information')
    parser.add_argument('-m', '--min_aligned', metavar='INT', type=int, required=False, default=20, help='minimum number of aligned bases to consider read as evidence [20]')
    parser.add_argument('-n', dest='num_samp', metavar='INT', type=int, required=False, default=1000000, help='number of reads to sample from BAM file for building insert size distribution [1000000]')
    parser.add_argument('-q', '--sum_quals', action='store_true', required=False, help='add genotyping quality to existing QUAL (default: overwrite QUAL field)')
    parser.add_argument('--max_reads', metavar='INT', type=int, default=None, required=False, help='maximum number of reads to assess at any variant (reduces processing time in high-depth regions, default: unlimited)')
    parser.add_argument('--split_weight', metavar='FLOAT', type=float, required=False, default=1, help='weight for split reads [1]')
    parser.add_argument('--disc_weight', metavar='FLOAT', type=float, required=False, default=1, help='weight for discordant paired-end reads [1]')
    parser.add_argument('-w', '--write_alignment', metavar='FILE', dest='alignment_outpath', type=str, required=False, default=None, help='write relevant reads to BAM file')
    parser.add_argument('--debug', action='store_true', help=argparse.SUPPRESS)

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input_vcf == None:
        if not sys.stdin.isatty():
            args.input_vcf = sys.stdin

    # send back the user input
    return args

# ==================================================
# Memoization Helpers
# ==================================================
def memoize_bam_search(func):
    cache = {}

    @wraps(func)
    def wrap(*args):
        # skip the fragment_dict input to gather_reads
        (sample, chrom, lpos, rpos) = args
        inputs = (sample.name, chrom, lpos, rpos)
        if inputs not in cache:
            sys.stderr.write('inovking bam search fn -- {} -- '.format(inputs))
            cache[inputs] = func(*args)
            sys.stderr.write('length: {} \n'.format(len(cache[inputs])))
        else:
            sys.stderr.write('using cached results -- {}\n'.format(inputs))
        return cache[inputs]

    return wrap

# ==================================================
# VCF parsing tools
# ==================================================

class Vcf(object):
    def __init__(self):
        self.file_format = 'VCFv4.2'
        # self.fasta = fasta
        self.reference = ''
        self.sample_list = []
        self.info_list = []
        self.format_list = []
        self.alt_list = []
        self.add_format('GT', 1, 'String', 'Genotype')
        self.header_misc = []

    def add_header(self, header):
        for line in header:
            if line.split('=')[0] == '##fileformat':
                self.file_format = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##reference':
                self.reference = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##INFO':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_info(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##ALT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_alt(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##FORMAT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_format(*[b.split('=')[1] for b in r.findall(a)])
            elif line[0] == '#' and line[1] != '#':
                self.sample_list = line.rstrip().split('\t')[9:]
            elif line.startswith('##fileDate='):
                pass
            else:
                self.header_misc.append(line.rstrip())

    # return the VCF header
    def get_header(self):
        header = '\n'.join(['##fileformat=' + self.file_format,
                            '##fileDate=' + time.strftime('%Y%m%d'),
                            '##reference=' + self.reference] + \
                           [i.hstring for i in self.info_list] + \
                           [a.hstring for a in self.alt_list] + \
                           [f.hstring for f in self.format_list] + \
                           self.header_misc + \
                           ['\t'.join([
                               '#CHROM',
                               'POS',
                               'ID',
                               'REF',
                               'ALT',
                               'QUAL',
                               'FILTER',
                               'INFO',
                               'FORMAT'] + \
                                      self.sample_list
                                  )])
        return header

    def add_info(self, id, number, type, desc):
        if id not in [i.id for i in self.info_list]:
            inf = self.Info(id, number, type, desc)
            self.info_list.append(inf)

    def add_alt(self, id, desc):
        if id not in [a.id for a in self.alt_list]:
            alt = self.Alt(id, desc)
            self.alt_list.append(alt)

    def add_format(self, id, number, type, desc):
        if id not in [f.id for f in self.format_list]:
            fmt = self.Format(id, number, type, desc)
            self.format_list.append(fmt)

    def add_sample(self, name):
        self.sample_list.append(name)

    # get the VCF column index of a sample
    # NOTE: this is zero-based, like python arrays
    def sample_to_col(self, sample):
        return self.sample_list.index(sample) + 9

    class Info(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##INFO=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

    class Alt(object):
        def __init__(self, id, desc):
            self.id = str(id)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##ALT=<ID=' + self.id + ',Description=\"' + self.desc + '\">'

    class Format(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##FORMAT=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

class Variant(object):
    def __init__(self, var_list, vcf):
        self.chrom = var_list[0]
        self.pos = int(var_list[1])
        self.var_id = var_list[2]
        self.ref = var_list[3]
        self.alt = var_list[4]
        if var_list[5] == '.':
            self.qual = 0
        else:
            self.qual = float(var_list[5])
        self.filter = var_list[6]
        self.sample_list = vcf.sample_list
        self.info_list = vcf.info_list
        self.info = dict()
        self.format_list = vcf.format_list
        self.active_formats = list()
        self.gts = dict()

        # fill in empty sample genotypes
        if len(var_list) < 8:
            sys.stderr.write('Error: VCF file must have at least 8 columns\n')
            exit(1)
        if len(var_list) < 9:
            var_list.append("GT")

        # make a genotype for each sample at variant
        for s in self.sample_list:
            try:
                s_gt = var_list[vcf.sample_to_col(s)].split(':')[0]
                self.gts[s] = Genotype(self, s, s_gt)
                # import the existing fmt fields
                for j in zip(var_list[8].split(':'), var_list[vcf.sample_to_col(s)].split(':')):
                    self.gts[s].set_format(j[0], j[1])
            except IndexError:
                self.gts[s] = Genotype(self, s, './.')

        self.info = dict()
        i_split = [a.split('=') for a in var_list[7].split(';')] # temp list of split info column
        for i in i_split:
            if len(i) == 1:
                i.append(True)
            self.info[i[0]] = i[1]

    def set_info(self, field, value):
        if field in [i.id for i in self.info_list]:
            self.info[field] = value
        else:
            sys.stderr.write('Error: invalid INFO field, \"' + field + '\"\n')
            exit(1)

    def get_info(self, field):
        return self.info[field]

    def get_info_string(self):
        i_list = list()
        for info_field in self.info_list:
            if info_field.id in self.info.keys():
                if info_field.type == 'Flag':
                    i_list.append(info_field.id)
                else:
                    i_list.append('%s=%s' % (info_field.id, self.info[info_field.id]))
        return ';'.join(i_list)

    def get_format_string(self):
        f_list = list()
        for f in self.format_list:
            if f.id in self.active_formats:
                f_list.append(f.id)
        return ':'.join(f_list)

    def genotype(self, sample_name):
        if sample_name in self.sample_list:
            return self.gts[sample_name]
        else:
            sys.stderr.write('Error: invalid sample name, \"' + sample_name + '\"\n')

    def get_var_string(self):
        s = '\t'.join(map(str,[
            self.chrom,
            self.pos,
            self.var_id,
            self.ref,
            self.alt,
            '%0.2f' % self.qual,
            self.filter,
            self.get_info_string(),
            self.get_format_string(),
            '\t'.join(self.genotype(s).get_gt_string() for s in self.sample_list)
        ]))
        return s

class Genotype(object):
    def __init__(self, variant, sample_name, gt):
        self.format = dict()
        self.variant = variant
        self.set_format('GT', gt)

    def set_format(self, field, value):
        if field in [i.id for i in self.variant.format_list]:
            self.format[field] = value
            if field not in self.variant.active_formats:
                self.variant.active_formats.append(field)
                # sort it to be in the same order as the format_list in header
                self.variant.active_formats.sort(key=lambda x: [f.id for f in self.variant.format_list].index(x))
        else:
            sys.stderr.write('Error: invalid FORMAT field, \"' + field + '\"\n')
            exit(1)

    def get_format(self, field):
        return self.format[field]

    def get_gt_string(self):
        g_list = list()
        for f in self.variant.active_formats:
            if f in self.format:
                if type(self.format[f]) == float:
                    g_list.append('%0.2f' % self.format[f])
                else:
                    g_list.append(self.format[f])
            else:
                g_list.append('.')
        return ':'.join(map(str,g_list))

# ==================================================
# BAM and JSON output
# ==================================================

# write read to BAM file, checking whether read is already written
def write_alignment(read, bam, written_reads, is_alt=None):
    read.query_sequence = None
    read_hash = (read.query_name, read.flag)

    if bam is None or read_hash in written_reads:
        return written_reads
    else:
        bam.write(read)
        written_reads.add(read_hash)
        return written_reads

# dump the sample and library info to a file
def write_sample_json(sample_list, lib_info_file):
    lib_info = {}
    for sample in sample_list:
        s = {}
        s['sample_name'] = sample.name
        s['bam'] = sample.bam.filename
        s['libraryArray'] = []
        s['mapped'] = sample.bam.mapped
        s['unmapped'] = sample.bam.unmapped

        for lib in sample.lib_dict.values():
            l = {}
            l['library_name'] = lib.name
            l['readgroups'] = lib.readgroups
            l['read_length'] = lib.read_length
            l['mean'] = lib.mean
            l['sd'] = lib.sd
            l['prevalence'] = lib.prevalence
            l['histogram'] = lib.hist

            s['libraryArray'].append(l)

        lib_info[sample.name] = s

    # write the json file
    json.dump(lib_info, lib_info_file, indent=4)
    lib_info_file.close()


# ==================================================
# Statistical tools
# ==================================================

# efficient combinatorial function to handle extremely large numbers
def log_choose(n, k):
    r = 0.0
    # swap for efficiency if k is more than half of n
    if k * 2 > n:
        k = n - k

    for  d in xrange(1,k+1):
        r += math.log(n, 10)
        r -= math.log(d, 10)
        n -= 1

    return r

# return the genotype and log10 p-value
def bayes_gt(ref, alt, is_dup):
    # probability of seeing an alt read with true genotype of of hom_ref, het, hom_alt respectively
    if is_dup: # specialized logic to handle non-destructive events such as duplications
        p_alt = [1e-2, 1/3.0, 0.5]
    else:
        p_alt = [1e-3, 0.5, 0.9]

    total = ref + alt
    log_combo = log_choose(total, alt)

    lp_homref = log_combo + alt * math.log(p_alt[0], 10) + ref * math.log(1 - p_alt[0], 10)
    lp_het = log_combo + alt * math.log(p_alt[1], 10) + ref * math.log(1 - p_alt[1], 10)
    lp_homalt = log_combo + alt * math.log(p_alt[2], 10) + ref * math.log(1 - p_alt[2], 10)

    return (lp_homref, lp_het, lp_homalt)

# get the number of entries in the set
def countRecords(myCounter):
    numRecords = sum(myCounter.values())
    return numRecords

# median is approx 50th percentile, except when it is between
# two values in which case it's the mean of them.
def median(myCounter):
    #length is the number of bases we're looking at
    numEntries = countRecords(myCounter)

    # the ordinal value of the middle element
    # if 2 middle elements, then non-integer
    limit = 0.5 * numEntries

    # a list of the values, sorted smallest to largest
    # note that this list contains unique elements only
    valueList = list(myCounter)
    valueList.sort()

    # number of entries we've gone through
    runEntries = 0
    # index of the current value in valueList
    i = 0
    # initiate v, in case list only has one element
    v = valueList[i]

    # move through the value list, iterating by number of
    # entries for each value
    while runEntries < limit:
        v = valueList[i]
        runEntries += myCounter[v]
        i += 1
    if runEntries == limit:
        return (v + valueList[i]) / 2.0
    else:
        return v

# calculate upper median absolute deviation
def upper_mad(myCounter, myMedian):
    residCounter = Counter()
    for x in myCounter:
        if x > myMedian:
            residCounter[abs(x - myMedian)] += myCounter[x]
    return median(residCounter)

# sum of the entries
def sumRecords(myCounter):
    mySum = 0.0
    for c in myCounter:
        mySum += c * float(myCounter[c])
    return mySum

# calculate the arithmetic mean, given a counter and the
# length of the feature (chromosome or genome)
# for x percentile, x% of the elements in the set are
# <= the output value
def mean(myCounter):
    # the number of total entries in the set is the
    # sum of the occurrences for each value
    numRecords = countRecords(myCounter)

    # u holds the mean
    u = float()

    u = sumRecords(myCounter) / numRecords
    return u

def stdev(myCounter):
    # the number of total entries in the set is the
    # sum of the occurrences for each value
    numRecords = countRecords(myCounter)

    # u holds the mean
    u = mean(myCounter)
    sumVar = 0.0

    # stdev is sqrt(sum((x-u)^2)/#elements)
    for c in myCounter:
        sumVar += myCounter[c] * (c - u)**2
    myVariance = float(sumVar) / numRecords
    stdev = myVariance**(0.5)
    return stdev

# ==================================================
# Library parsing
# ==================================================

# holds a library's insert size and read length information
class Library(object):
    def __init__(self,
                 name,
                 bam,
                 readgroups,
                 read_length,
                 hist,
                 dens,
                 mean,
                 sd,
                 prevalence,
                 num_samp):

        # parse arguments
        self.name = name
        self.bam = bam
        self.num_samp = num_samp
        self.readgroups = readgroups
        self.read_length = read_length
        self.hist = hist
        self.dens = dens
        self.mean = mean
        self.sd = sd
        self.prevalence = prevalence

        # if information is missing, compute it
        if self.read_length is None:
            self.calc_read_length()
        if self.hist is None:
            self.calc_insert_hist()
        if self.dens is None:
            self.calc_insert_density()
        if self.prevalence is None:
            self.calc_lib_prevalence()

    @classmethod
    def from_lib_info(cls,
                      sample_name,
                      lib_index,
                      bam,
                      lib_info):

        lib = lib_info[sample_name]['libraryArray'][lib_index]

        # convert the histogram keys to integers (from strings in JSON)
        lib_hist = {int(k):int(v) for k,v in lib['histogram'].items()}

        return cls(lib['library_name'],
                   bam,
                   lib['readgroups'],
                   int(lib['read_length']),
                   lib_hist,
                   None,
                   float(lib['mean']),
                   float(lib['sd']),
                   float(lib['prevalence']),
                   0)

    @classmethod
    def from_bam(cls,
                 lib_name,
                 bam,
                 num_samp):

        # get readgroups that comprise the library
        readgroups = []
        for r in bam.header['RG']:
            try:
                in_lib = r['LB'] == lib_name
            except KeyError, e:
                in_lib = lib_name == ''

            if in_lib:
                readgroups.append(r['ID'])

        return cls(lib_name,
                   bam,
                   readgroups,
                   None,
                   None,
                   None,
                   None,
                   None,
                   None,
                   num_samp)

    # calculate the library's prevalence in the BAM file
    def calc_lib_prevalence(self):
        max_count = 100000
        lib_counter = 0
        read_counter = 0

        for read in self.bam.fetch():
            if read_counter == max_count:
                break
            if read.get_tag('RG') in self.readgroups:
                lib_counter += 1
            read_counter += 1

        self.prevalence = float(lib_counter) / read_counter

    # get read length
    def calc_read_length(self):
        max_rl = 0
        counter = 0
        num_samp = 10000
        for read in self.bam.fetch():
            if read.get_tag('RG') not in self.readgroups:
                continue
            if read.infer_query_length() > max_rl:
                max_rl = read.infer_query_length()
            if counter == num_samp:
                break
            counter += 1
        self.read_length = max_rl

    # generate empirical histogram of the sample's insert size distribution
    def calc_insert_hist(self):
        counter = 0
        skip = 0
        skip_counter = 0
        mads = 10
        ins_list = []

        # Each entry in valueCounts is a value, and its count is
        # the number of instances of that value observed in the dataset.
        # So valueCount[5] is the number of times 5 has been seen in the data.
        valueCounts = Counter()
        for read in self.bam.fetch():
            if skip_counter < skip:
                skip_counter += 1
                continue
            if (read.is_reverse
                or not read.mate_is_reverse
                or read.is_unmapped
                or read.mate_is_unmapped
                or not is_primary(read)
                or read.template_length <= 0
                or read.get_tag('RG') not in self.readgroups):
                continue
            else:
                valueCounts[read.template_length] += 1
                counter += 1
            if counter == self.num_samp:
                break

        if len(valueCounts) == 0:
            sys.stderr.write('Error: failed to build insert size histogram for paired-end reads.\n\
Please ensure BAM file (%s) has inward facing, paired-end reads.\n' % self.bam.filename)
            exit(1)

        # remove outliers
        med = median(valueCounts)
        u_mad = upper_mad(valueCounts, med)
        for x in [x for x in list(valueCounts) if x > med + mads * u_mad]:
            del valueCounts[x]

        self.hist = valueCounts
        self.mean = mean(self.hist)
        self.sd = stdev(self.hist)

    # calculate the density curve for and insert size histogram
    def calc_insert_density(self):
        dens = Counter()
        for i in list(self.hist):
            dens[i] = float(self.hist[i])/countRecords(self.hist)
        self.dens = dens

# ==================================================
# Sample parsing
# ==================================================

# holds each sample's BAM and library information
class Sample(object):
    # general constructor
    def __init__(self,
                 name,
                 bam,
                 num_samp,
                 lib_dict,
                 rg_to_lib,
                 min_lib_prevalence,
                 bam_mapped,
                 bam_unmapped):

        self.name = name
        self.bam = bam
        self.lib_dict = lib_dict
        self.rg_to_lib = rg_to_lib
        self.bam_mapped = bam_mapped
        self.bam_unmapped = bam_unmapped

        # get active libraries
        self.active_libs = []
        for lib in lib_dict.values():
            if lib.prevalence >= min_lib_prevalence:
                self.active_libs.append(lib)

    # constructor from supplied JSON descriptor
    @classmethod
    def from_lib_info(cls,
                      bam,
                      lib_info,
                      min_lib_prevalence):

        name = bam.header['RG'][0]['SM']
        num_samp = 0
        rg_to_lib = {}
        lib_dict = {}

        try:
            for i in xrange(len(lib_info[name]['libraryArray'])):
                lib = lib_info[name]['libraryArray'][i]
                lib_name = lib['library_name']

                # make library object
                lib_dict[lib_name] = Library.from_lib_info(name,
                                                           i,
                                                           bam,
                                                           lib_info)

                # make a map from readgroup IDs to library objects
                for rg in lib['readgroups']:
                    rg_to_lib[rg] = lib_dict[lib_name]
        except KeyError:
            sys.stderr.write('Error: sample %s not found in JSON library file.\n' % name)
            exit(1)

        return cls(name,
                   bam,
                   num_samp,
                   lib_dict,
                   rg_to_lib,
                   min_lib_prevalence,
                   lib_info[name]['mapped'],
                   lib_info[name]['unmapped'])

    # constructor for empirical distributions
    @classmethod
    def from_bam(cls,
                 bam,
                 num_samp,
                 min_lib_prevalence):

        name = bam.header['RG'][0]['SM']
        rg_to_lib = {}
        lib_dict = {}

        for r in bam.header['RG']:
            try:
                lib_name=r['LB']
            except KeyError, e:
                lib_name=''

            # add the new library
            if lib_name not in lib_dict:
                new_lib = Library.from_bam(lib_name, bam, num_samp)
                lib_dict[lib_name] = new_lib
            rg_to_lib[r['ID']] = lib_dict[lib_name]

        return cls(name,
                   bam,
                   num_samp,
                   lib_dict,
                   rg_to_lib,
                   min_lib_prevalence,
                   bam.mapped,
                   bam.unmapped)

    # get the maximum fetch flank for reading the BAM file
    def get_fetch_flank(self, z):
        return max([lib.mean + (lib.sd * z) for lib in self.lib_dict.values()])

    # return the library object for a specified read group
    def get_lib(self, readgroup):
        return self.rg_to_lib[readgroup]

    # return the expected spanning coverage at any given base
    def set_exp_spanning_depth(self, min_aligned):
        genome_size = float(sum(self.bam.lengths))
        weighted_mean_span = sum([(lib.mean - 2 * lib.read_length + 2 * min_aligned) * lib.prevalence for lib in self.lib_dict.values()])
        exp_spanning_depth = (weighted_mean_span * self.bam_mapped) / genome_size

        self.exp_spanning_depth = exp_spanning_depth

        return

    # return the expected sequence coverage at any given base
    def set_exp_seq_depth(self, min_aligned):
        genome_size = float(sum(self.bam.lengths))
        weighted_mean_read_length = sum([(lib.read_length - 2 * min_aligned) * lib.prevalence for lib in self.lib_dict.values()])
        exp_seq_depth = (weighted_mean_read_length * self.bam_mapped) / genome_size

        self.exp_seq_depth = exp_seq_depth

        return

# ==================================================
# Class for SAM fragment, containing all alignments
# from a single molecule
# ==================================================

class SamFragment(object):
    def __init__(self, read, lib):
        self.lib = lib
        self.primary_reads = []
        self.split_reads = []
        self.read_set = set()
        self.num_primary = 0
        self.query_name = read.query_name

        self.readA = None
        self.readB = None
        self.ispan = None
        self.ospan = None

        self.add_read(read)

    def add_read(self, read):
        # ensure we don't add the same read twice
        read_hash = read.__hash__()
        if read_hash in self.read_set:
            return
        else:
            self.read_set.add(read_hash)

        if is_primary(read):
            # add the primary reads
            self.primary_reads.append(read)
            self.num_primary += 1

            # make split candidate and check whether it's valid
            split_candidate = SplitRead(read, self.lib)
            if split_candidate.is_valid():
                self.split_reads.append(split_candidate)

            # complete set of primaries
            if self.num_primary == 2:
                self.readA, self.readB = self.primary_reads

    # tag the read with R (ref), A (alt), or U (unknown) XV tag
    def tag_span(self, p_alt=None):
        if p_alt is None:
            value = 'U'
        elif p_alt > 0:
            value = 'A'
        else:
            value = 'R'
        for read in self.primary_reads:
            if not read.has_tag('XV'):
                read.set_tag('XV', value)

        return

    # get inner span
    def get_ispan(self, min_aligned):
        ispan1 = self.readA.reference_start + min_aligned
        ispan2 = self.readB.reference_end - min_aligned - 1

        return (ispan1, ispan2)

    # get outer span
    def get_ospan(self):
        ospan1 = self.readA.reference_start
        ospan2 = self.readB.reference_end

        return (ospan1, ospan2)

    # returns boolean of whether a single read crosses
    # a genomic reference point with appropriate flanking
    # sequence
    def is_ref_seq(self,
                   read,
                   chrom, pos, ci,
                   min_aligned):
        # check chromosome matching
        # Note: this step is kind of slow
        if read.reference_name != chrom:
            return False

        # ensure there is min_aligned on both sides of position
        if read.get_overlap(max(0, pos - min_aligned), pos + min_aligned) < 2 * min_aligned:
            return False

        return True


    # returns boolean of whether the primary pair of a
    # fragment straddles a genomic point
    def is_pair_straddle(self,
                         chromA, posA, ciA,
                         chromB, posB, ciB,
                         o1_is_reverse, o2_is_reverse,
                         min_aligned,
                         lib):
        if self.num_primary != 2:
            return False

        # check orientation
        if self.readA.is_reverse != o1_is_reverse:
            return False
        if self.readB.is_reverse != o2_is_reverse:
            return False

        # check chromosome matching
        # Note: this step is kind of slow
        if self.readA.reference_name != chromA:
            return False
        if self.readB.reference_name != chromB:
            return False

        # get the inner span
        ispan = self.get_ispan(min_aligned)

        # check orientations and positions
        flank = lib.mean + lib.sd * 3
        if not o1_is_reverse and (ispan[0] > posA + ciA[1] or ispan[0] < posA + ciA[0] - flank):
            return False
        if o1_is_reverse and (ispan[0] < posA + ciA[0] or ispan[0] > posA + ciA[1] + flank):
            return False
        if not o2_is_reverse and (ispan[1] > posB + ciB[1] or ispan[1] < posB + ciB[0] - flank):
            return False
        if o2_is_reverse and (ispan[1] < posB + ciB[0] or ispan[1] > posB + ciB[1] + flank):
            return False

        return True

    # calculate the probability that a read pair is concordant at a breakpoint,
    # given the putative variant size and insert distribution of the library.
    def p_concordant(self, var_length=None):
        # a priori probability that a read-pair is concordant
        disc_prior = 0.05
        conc_prior = 1 - disc_prior

        ospan = self.get_ospan()

        # outer span length
        ospan_length = abs(ospan[1] - ospan[0])

        # if no variant length (such as in the case of a non-deletion variant)
        # default to mean plus 3 stdev
        z = 3
        if var_length is None:
            var_length = self.lib.mean + self.lib.sd * z

        try:
            p = float(self.lib.dens[ospan_length]) * conc_prior / (conc_prior * self.lib.dens[ospan_length] + disc_prior * (self.lib.dens[ospan_length - var_length]))
        except ZeroDivisionError:
            p = None

        return p > 0.5

# ==================================================
# Class for a split-read, containing all alignments
# from a single chimeric read
# ==================================================

# each SplitRead object has a left and a right SplitPiece
# (reads with more than 2 split alignments are discarded)
class SplitRead(object):
    def __init__(self, read, lib):
        self.query_name = read.query_name
        self.read = read
        self.lib = lib
        self.sa = None
        self.q1 = None
        self.q2 = None
        self.is_soft_clip = False

    # the piece of each split alignemnt
    class SplitPiece(object):
        def __init__(self,
                     chrom,
                     reference_start,
                     is_reverse,
                     cigar,
                     mapq):
            self.chrom = chrom
            self.reference_start = reference_start
            self.reference_end = None
            self.is_reverse = is_reverse
            self.cigar = cigar
            self.mapq = mapq
            self.left_query = None
            self.right_query = None

            # get query positions
            self.query_pos = get_query_pos_from_cigar(self.cigar, self.is_reverse)

        def set_reference_end(self, reference_end):
            self.reference_end = reference_end

    # check if passes QC, and populate with necessary info
    def is_valid(self,
                 min_non_overlap = 20,
                 min_indel = 50,
                 max_unmapped_bases = 50):
        # check for SA tag
        if not self.read.has_tag('SA'):
            # Include soft-clipped reads that didn't generate split-read alignments
            if is_clip_op(self.read.cigar[0][0]) or is_clip_op(self.read.cigar[-1][0]):
                clip_length = max(self.read.cigar[0][1] * is_clip_op(self.read.cigar[0][0]), self.read.cigar[-1][1] * is_clip_op(self.read.cigar[-1][0]))
                # Only count if the longest clipping event is greater than the cutoff and we have mapped a reasonable number of bases
                if clip_length > 0 and (self.read.query_length - self.read.query_alignment_length) <= max_unmapped_bases:
                    a = self.SplitPiece(self.read.reference_name,
                                        self.read.reference_start,
                                        self.read.is_reverse,
                                        self.read.cigar,
                                        self.read.mapping_quality)
                    a.set_reference_end(self.read.reference_end)
                    b = self.SplitPiece(None,
                                        1,
                                        self.read.is_reverse,
                                        self.read.cigar,
                                        0)
                    b.set_reference_end(1)
                    self.set_order_by_clip(a, b)
                    self.is_soft_clip = True
                    return True
                else:
                    return False

            return False

        # parse SA tag
        sa_list = self.read.get_tag('SA').rstrip(';').split(';')
        if len(sa_list) > 1:
            return False
        else:
            self.sa = sa_list[0].split(',')
            mate_chrom = self.sa[0]
            mate_pos = int(self.sa[1]) - 1 # SA tag is one-based, while SAM is zero-based
            mate_is_reverse = self.sa[2] == '-'
            mate_cigar = cigarstring_to_tuple(self.sa[3])
            mate_mapq = int(self.sa[4])

        # make SplitPiece objects
        a = self.SplitPiece(self.read.reference_name,
                            self.read.reference_start,
                            self.read.is_reverse,
                            self.read.cigar,
                            self.read.mapping_quality)
        a.set_reference_end(self.read.reference_end)

        b = self.SplitPiece(mate_chrom,
                            mate_pos,
                            mate_is_reverse,
                            mate_cigar,
                            mate_mapq)
        b.set_reference_end(get_reference_end_from_cigar(b.reference_start, b.cigar))

        # set query_left and query_right splitter by alignment position on the reference
        # this is used for non-overlap and off-diagonal filtering
        # (query_left and query_right are random when on different chromosomes
        if self.read.reference_name == mate_chrom:
            if self.read.pos > mate_pos:
                self.query_left = b
                self.query_right = a
            else:
                self.query_left = a
                self.query_right = b
        else:
            self.set_order_by_clip(a, b)

        # check non-overlap
        if self.non_overlap() < min_non_overlap:
            return False

        # check off-diagonal distance and desert
        # only relevant when split pieces are on the same chromosome and strand
        if (self.query_left.chrom == self.query_right.chrom
            and self.query_left.is_reverse == self.query_right.is_reverse):
            # use end diagonal on left and start diagonal on right since
            # the start and end diags might be different if there is an
            # indel in the alignments
            if self.query_left.is_reverse:
                left_diag = get_start_diagonal(self.query_left)
                right_diag = get_end_diagonal(self.query_right)
                ins_size = right_diag - left_diag
            else:
                left_diag = get_end_diagonal(self.query_left)
                right_diag = get_start_diagonal(self.query_right)
                ins_size = left_diag - right_diag
            if abs(ins_size) < min_indel:
                return False

            # check for desert gap of indels
            desert = self.query_right.query_pos.query_start - self.query_left.query_pos.query_end - 1
            if desert > 0 and desert - max(0, ins_size) > max_unmapped_bases:
                return False

        # passed all checks. valid split-read
        return True

    def non_overlap(self):
        # get overlap of aligned query positions
        overlap = get_query_overlap(self.query_left.query_pos.query_start,
                                    self.query_left.query_pos.query_end,
                                    self.query_right.query_pos.query_start,
                                    self.query_right.query_pos.query_end)

        # get minimum non-overlap
        left_non_overlap = 1 + self.query_left.query_pos.query_end - self.query_left.query_pos.query_start - overlap
        right_non_overlap = 1 + self.query_right.query_pos.query_end - self.query_right.query_pos.query_start - overlap
        non_overlap = min(left_non_overlap, right_non_overlap)

        return non_overlap

    @staticmethod
    def check_split_support(split, chrom, pos, is_reverse, split_slop):
        if split.chrom != chrom:
            return False

        if is_reverse:
            coord = split.reference_start
        else:
            coord = split.reference_end

        if (coord > pos + split_slop
                or coord < pos - split_slop):
            return False
        return True

    def is_split_straddle(self,
                          chromA, posA, ciA,
                          chromB, posB, ciB,
                          o1_is_reverse, o2_is_reverse,
                          svtype, split_slop):

        # arrange the SV breakends from left to right
        if (chromA != chromB
            or (chromA == chromB and posA > posB)):
            chrom_left = chromB
            pos_left = posB
            ci_left = ciB
            is_reverse_left = o2_is_reverse
            chrom_right = chromA
            pos_right = posA
            ci_right = ciA
            is_reverse_right = o1_is_reverse
        else:
            chrom_left = chromA
            pos_left = posA
            ci_left = ciA
            is_reverse_left = o1_is_reverse
            chrom_right = chromB
            pos_right = posB
            ci_right = ciB
            is_reverse_right = o2_is_reverse


        # check split chromosomes against variant
        left_split = False
        right_split = False

        if (not self.is_soft_clip) or svtype == 'DEL' or svtype == 'INS':
            left_split = self.check_split_support(self.query_left,
                    chrom_left,
                    pos_left,
                    is_reverse_left,
                    split_slop)
            right_split = self.check_split_support(self.query_right,
                    chrom_right,
                    pos_right,
                    is_reverse_right,
                    split_slop)
        elif svtype == 'DUP':
            left_split = self.check_split_support(self.query_left,
                    chrom_right,
                    pos_right,
                    is_reverse_right,
                    split_slop)
            right_split = self.check_split_support(self.query_right,
                    chrom_left,
                    pos_left,
                    is_reverse_left,
                    split_slop)
        elif svtype == 'INV':
                # check all possible sides
            left_split_left = self.check_split_support(self.query_left,
                    chrom_left,
                    pos_left,
                    is_reverse_left,
                    split_slop)
            left_split_right = self.check_split_support(self.query_left,
                    chrom_right,
                    pos_right,
                    is_reverse_right,
                    split_slop)
            left_split = left_split_left or left_split_right
            right_split_left = self.check_split_support(self.query_right,
                    chrom_left,
                    pos_left,
                    is_reverse_left,
                    split_slop)
            right_split_right = self.check_split_support(self.query_right,
                    chrom_right,
                    pos_right,
                    is_reverse_right,
                    split_slop)
            right_split = right_split_left or right_split_right

        return (left_split, right_split)

    # tag the read with R (ref), A (alt), or U (unknown) XV tag
    def tag_split(self, p_alt=None):
        if p_alt is None:
            value = 'U'
        elif p_alt > 0:
            value = 'A'
        else:
            value = 'R'

        self.read.set_tag('XV', value)

        return

    def set_order_by_clip(self, a, b):
        '''
        Determine which SplitPiece is the leftmost based
        on the side of the longest clipping operation
        '''
        if is_left_clip(a.cigar):
            self.query_left = b
            self.query_right = a
        else:
            self.query_left = a
            self.query_right = b


# ==================================================
# Miscellaneous methods for manipulating SAM alignments
# ==================================================

# adapted from Matt Shirley (http://coderscrowd.com/app/public/codes/view/171)
def cigarstring_to_tuple(cigarstring):
    cigar_dict = {'M':0, 'I':1,'D':2,'N':3, 'S':4, 'H':5, 'P':6, '=':7, 'X':8}
    pattern = re.compile('([MIDNSHPX=])')
    values = pattern.split(cigarstring)[:-1] ## turn cigar into tuple of values
    paired = (values[n:n+2] for n in xrange(0, len(values), 2)) ## pair values by twos
    return [(cigar_dict[pair[1]], int(pair[0])) for pair in paired]

def is_left_clip(cigar):
    '''
    whether the left side of the read (w/ respect to reference) is clipped.
    Clipping side is determined as the side with the longest clip.
    Adjacent clipping operations are not considered
    '''
    left_tuple = cigar[0]
    right_tuple = cigar[-1]
    left_clipped = is_clip_op(left_tuple[0])
    right_clipped = is_clip_op(right_tuple[0])

    return (left_clipped and not right_clipped) or (left_clipped and right_clipped and left_tuple[1] > right_tuple[1])

def is_clip_op(op):
    '''
    whether the CIGAR OP code represents a clipping event
    '''
    return op == 4 or op == 5

# reference position where the alignment would have started
# if the entire query sequence would have aligned
def get_start_diagonal(split_piece):
    sclip = split_piece.query_pos.query_start
    if split_piece.is_reverse:
        sclip = split_piece.query_pos.query_length - split_piece.query_pos.query_end
    return split_piece.reference_start - sclip

# reference position where the alignment would have ended
# if the entire query sequence would have aligned
def get_end_diagonal(split_piece):
    query_aligned = split_piece.query_pos.query_end
    if split_piece.is_reverse:
        query_aligned = split_piece.query_pos.query_length - split_piece.query_pos.query_start
    return split_piece.reference_end - query_aligned

def get_reference_end_from_cigar(reference_start, cigar):
    '''
    This returns the coordinate just past the last aligned base.
    This matches the behavior of pysam's reference_end method
    '''
    reference_end = reference_start

    # iterate through cigartuple
    for i in xrange(len(cigar)):
        k, n = cigar[i]
        if k in (0,2,3,7,8): # M, D, N, =, X
            reference_end += n
    return reference_end

# get the positions of the query that are aligned
def get_query_pos_from_cigar(cigar, is_reverse):
    query_start = 0
    query_end = 0
    query_length = 0

    # flip if negative strand
    if is_reverse:
        cigar = cigar[::-1]

    # iterate through cigartuple
    for i in xrange(len(cigar)):
        k, n = cigar[i]
        if k in (4,5): # H, S
            if i == 0:
                query_start += n
                query_end += n
                query_length += n
            else:
                query_length += n
        elif k in (0,1,7,8): # M, I, =, X
            query_end += n
            query_length +=n

    d = QueryPos(query_start, query_end, query_length);
    return d

# structure to hold query position information
class QueryPos (object):
    """
    struct to store the start and end positions of query CIGAR operations
    """
    def __init__(self, query_start, query_end, query_length):
        self.query_start = int(query_start)
        self.query_end = int(query_end)
        self.query_length  = int(query_length)

def get_query_overlap(s1, e1, s2, e2):
    o = 1 + min(e1, e2) - max(s1, s2)
    return max(0, o)

# read is neither secondary nor supplementary
def is_primary(read):
    return (not read.is_supplementary and not read.is_secondary)

# get the non-phred-scaled mapq of a read
def prob_mapq(read):
    return 1 - 10 ** (-read.mapq / 10.0)

# methods to grab reads from region of interest in BAM file
def gather_all_reads(sample, chromA, posA, ciA, chromB, posB, ciB, z, max_reads):
    # grab batch of reads from both sides of breakpoint
    read_batch = {}
    read_batch, many = gather_reads(sample, chromA, posA, ciA, z, read_batch, max_reads)
    if many:
        return {}, True

    read_batch, many = gather_reads(sample, chromB, posB, ciB, z, read_batch, max_reads)
    if many:
        return {}, True

    return read_batch, many

def gather_reads(sample,
                 chrom, pos, ci,
                 z,
                 fragment_dict,
                 max_reads):

    # the distance to the left and right of the breakpoint to scan
    # (max of mean + z standard devs over all of a sample's libraries)
    fetch_flank = sample.get_fetch_flank(z)
    chrom_length = sample.bam.lengths[sample.bam.gettid(chrom)]

    many = False

    reads = fetch_reads_from_bam(
        sample,
        chrom,
        max(pos + ci[0] - fetch_flank, 0),
        min(pos + ci[1] + fetch_flank, chrom_length)
    )

    for i, read in enumerate(reads):
        if read.is_unmapped or read.is_duplicate:
            continue

        lib = sample.get_lib(read.get_tag('RG'))
        if lib not in sample.active_libs:
            continue

        # read.query_sequence = "*"
        # read.query_qualities = "*"
        if max_reads is not None and i > max_reads:
            many = True
            break

        if read.query_name in fragment_dict:
            fragment_dict[read.query_name].add_read(read)
        else:
            fragment_dict[read.query_name] = SamFragment(read, lib)

    return fragment_dict, many

@memoize_bam_search
def fetch_reads_from_bam(sample, chrom, left_pos, right_pos):
    return list(sample.bam.fetch(chrom, left_pos, right_pos))

# ==================================================
# Genotyping function
# ==================================================

def sv_genotype(bam_string,
                vcf_in,
                vcf_out,
                min_aligned,
                split_weight,
                disc_weight,
                num_samp,
                lib_info_path,
                debug,
                alignment_outpath,
                ref_fasta,
                sum_quals,
                max_reads):

    # parse the comma separated inputs
    bam_list = []
    for b in bam_string.split(','):
        if b.endswith('.bam'):
            bam_list.append(pysam.AlignmentFile(b, mode='rb'))
        elif b.endswith('.cram'):
            bam_list.append(pysam.AlignmentFile(b, mode='rc', reference_filename=ref_fasta))
        else:
            sys.stderr.write('Error: %s is not a valid alignment file (*.bam or *.cram)\n' % b)
            exit(1)
            
    min_lib_prevalence = 1e-3 # only consider libraries that constitute at least this fraction of the BAM

    # parse lib_info_path JSON
    lib_info = None
    if lib_info_path is not None and os.path.isfile(lib_info_path):
        lib_info_file = open(lib_info_path, 'r')
        lib_info = json.load(lib_info_file)

    if vcf_in is None:
        sys.stderr.write('Warning: VCF not found.\n')

    # build the sample libraries, either from the lib_info JSON or empirically from the BAMs
    sample_list = list()
    for i in xrange(len(bam_list)):
        if lib_info is None:
            sys.stderr.write('Calculating library metrics from %s...' % bam_list[i].filename)
            sample = Sample.from_bam(bam_list[i], num_samp, min_lib_prevalence)
        else:
            sys.stderr.write('Reading library metrics from %s...' % lib_info_path)
            sample = Sample.from_lib_info(bam_list[i], lib_info, min_lib_prevalence)

        sample.set_exp_seq_depth(min_aligned)
        sample.set_exp_spanning_depth(min_aligned)
        sample_list.append(sample)
    sys.stderr.write(' done\n')

    # diagnostic dump of relevant BAM reads
    if alignment_outpath is not None:
        # create a BAM file of the reads used for genotyping
        out_bam_written_reads = set()
        template_bam = pysam.AlignmentFile(bam_string.split(',')[0], 'rb')
        out_bam = pysam.AlignmentFile(alignment_outpath, 'wb', template_bam)
        template_bam.close()

    # write the JSON for each sample's libraries
    if lib_info_path is not None and not os.path.isfile(lib_info_path):
        sys.stderr.write('Writing library metrics to %s...' % lib_info_path)
        lib_info_file = open(lib_info_path, 'w')
        write_sample_json(sample_list, lib_info_file)
        lib_info_file.close()
        sys.stderr.write(' done\n')

    # quit early if VCF absent
    if vcf_in is None:
        if alignment_outpath is not None:
            out_bam.close()
        return

    # set variables for genotyping
    z = 3
    split_slop = 3 # amount of slop around breakpoint to count splitters
    in_header = True
    header = []
    breakend_dict = {} # cache to hold unmatched generic breakends for genotyping
    vcf = Vcf()

    # read input VCF
    for line in vcf_in:
        if in_header:
            if line[0] == '#':
                header.append(line)
                if line[1] != '#':
                    vcf_samples = line.rstrip().split('\t')[9:]
                continue
            else:
                in_header = False
                vcf.add_header(header)
                # if detailed:
                vcf.add_info('SVTYPE', 1, 'String', 'Type of structural variant')
                vcf.add_format('GQ', 1, 'Integer', 'Genotype quality')
                vcf.add_format('SQ', 1, 'Float', 'Phred-scaled probability that this site is variant (non-reference in this sample')
                vcf.add_format('GL', 'G', 'Float', 'Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy')
                vcf.add_format('DP', 1, 'Integer', 'Read depth')
                vcf.add_format('RO', 1, 'Integer', 'Reference allele observation count, with partial observations recorded fractionally')
                vcf.add_format('AO', 'A', 'Integer', 'Alternate allele observations, with partial observations recorded fractionally')
                vcf.add_format('QR', 1, 'Integer', 'Sum of quality of reference observations')
                vcf.add_format('QA', 'A', 'Integer', 'Sum of quality of alternate observations')
                vcf.add_format('RS', 1, 'Integer', 'Reference allele split-read observation count, with partial observations recorded fractionally')
                vcf.add_format('AS', 'A', 'Integer', 'Alternate allele split-read observation count, with partial observations recorded fractionally')
                vcf.add_format('ASC', 'A', 'Integer', 'Alternate allele clipped-read observation count, with partial observations recorded fractionally')
                vcf.add_format('RP', 1, 'Integer', 'Reference allele paired-end observation count, with partial observations recorded fractionally')
                vcf.add_format('AP', 'A', 'Integer', 'Alternate allele paired-end observation count, with partial observations recorded fractionally')
                vcf.add_format('AB', 'A', 'Float', 'Allele balance, fraction of observations from alternate allele, QA/(QR+QA)')


                # add the samples in the BAM files to the VCF output
                for sample in sample_list:
                    if sample.name not in vcf.sample_list:
                        vcf.add_sample(sample.name)

                # write the output header
                vcf_out.write(vcf.get_header() + '\n')


        v = line.rstrip().split('\t')
        var = Variant(v, vcf)
        var_length = None # var_length should be None except for deletions
        if not sum_quals:
            var.qual = 0

        # genotype generic breakends
        try:
            svtype = var.get_info('SVTYPE')
        except KeyError:
            sys.stderr.write('Warning: SVTYPE missing at variant %s. Skipping.\n' % (var.var_id))
            vcf_out.write(var.get_var_string() + '\n')
            continue
            
        # print original line if unsupported svtype
        if svtype not in ('BND', 'DEL', 'DUP', 'INV'):
            sys.stderr.write('Warning: Unsupported SVTYPE at variant %s (%s). Skipping.\n' % (var.var_id, svtype))
            vcf_out.write(var.get_var_string() + '\n')
            continue

        if svtype == 'BND':
            if var.info['MATEID'] in breakend_dict:
                var2 = var
                var = breakend_dict[var.info['MATEID']]
                chromA = var.chrom
                chromB = var2.chrom
                posA = var.pos
                posB = var2.pos
                # confidence intervals
                ciA = map(int, var.info['CIPOS'].split(','))
                ciB = map(int, var2.info['CIPOS'].split(','))

                # infer the strands from the alt allele
                if var.alt[-1] == '[' or var.alt[-1] == ']':
                    o1_is_reverse = False
                else: o1_is_reverse = True
                if var2.alt[-1] == '[' or var2.alt[-1] == ']':
                    o2_is_reverse = False
                else: o2_is_reverse = True

                # remove the BND from the breakend_dict
                # to free up memory
                del breakend_dict[var.var_id]
            else:
                breakend_dict[var.var_id] = var
                continue
        else:
            chromA = var.chrom
            chromB = var.chrom
            posA = var.pos
            posB = int(var.get_info('END'))
            # confidence intervals
            ciA = map(int, var.info['CIPOS'].split(','))
            ciB = map(int, var.info['CIEND'].split(','))
            if svtype == 'DEL':
                var_length = posB - posA
                o1_is_reverse, o2_is_reverse =  False, True
            elif svtype == 'DUP':
                o1_is_reverse, o2_is_reverse =  True, False
            elif svtype == 'INV':
                o1_is_reverse, o2_is_reverse = False, False

        # increment the negative strand values (note position in VCF should be the base immediately left of the breakpoint junction)
        if o1_is_reverse: posA += 1
        if o2_is_reverse: posB += 1

        for sample in sample_list:
            # grab reads from both sides of breakpoint
            read_batch, many = gather_all_reads(sample, chromA, posA, ciA, chromB, posB, ciB, z, max_reads)
            if many:
                var.genotype(sample.name).set_format('GT', './.')
                continue

            # initialize counts to zero
            ref_span, alt_span = 0, 0
            ref_seq, alt_seq = 0, 0
            alt_clip = 0

            # ref_ciA = ciA
            # ref_ciB = ciB
            ref_ciA = [0,0]
            ref_ciB = [0,0]

            for query_name in read_batch:
                fragment = read_batch[query_name]
                # boolean on whether to write the fragment
                write_fragment = False

                # -------------------------------------
                # Check for split-read evidence
                # -------------------------------------

                # get reference sequences
                for read in fragment.primary_reads:
                    if (fragment.is_ref_seq(read, chromA, posA, ciA, min_aligned)
                        or fragment.is_ref_seq(read, chromB, posB, ciB, min_aligned)):

                        p_reference = prob_mapq(read)
                        ref_seq += p_reference

                        read.set_tag('XV', 'R')
                        write_fragment = True

                # get non-reference split-read support
                for split in fragment.split_reads:

                    split_lr = split.is_split_straddle(chromA, posA, ciA,
                                                       chromB, posB, ciB,
                                                       o1_is_reverse, o2_is_reverse,
                                                       svtype, split_slop)
                    # p_alt = prob_mapq(split.query_left) * prob_mapq(split.query_right)
                    p_alt = (prob_mapq(split.query_left) * split_lr[0] + prob_mapq(split.query_right) * split_lr[1]) / 2.0
                    if split.is_soft_clip:
                        alt_clip += p_alt
                    else:
                        alt_seq += p_alt

                    if p_alt > 0:
                        split.tag_split(p_alt)
                        write_fragment = True

                # -------------------------------------
                # Check for paired-end evidence
                # -------------------------------------

                # tally spanning alternate pairs
                if svtype == 'DEL' and posB - posA < 2 * fragment.lib.sd:
                    alt_straddle = False
                else:
                    alt_straddle = fragment.is_pair_straddle(chromA, posA, ciA,
                                                             chromB, posB, ciB,
                                                             o1_is_reverse, o2_is_reverse,
                                                             min_aligned,
                                                             fragment.lib)

                # check both sides if inversion (perhaps should do this for BND as well?)
                if svtype in ('INV'):
                    alt_straddle_reciprocal = fragment.is_pair_straddle(chromA, posA, ciA,
                                                                        chromB, posB, ciB,
                                                                        not o1_is_reverse,
                                                                        not o2_is_reverse,
                                                                        min_aligned,
                                                                        fragment.lib)
                else:
                    alt_straddle_reciprocal = False

                if alt_straddle or alt_straddle_reciprocal:
                    if svtype == 'DEL':
                        p_conc = fragment.p_concordant(var_length)
                        if p_conc is not None:
                            p_alt = (1 - p_conc) * prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                            alt_span += p_alt

                            # # since an alt straddler is by definition also a reference straddler,
                            # # we can bail out early here to save some time
                            # p_reference = p_conc * prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                            # ref_span += p_reference
                            # continue

                            fragment.tag_span(p_alt)
                            write_fragment = True

                    else:
                        p_alt = prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                        alt_span += p_alt

                        fragment.tag_span(p_alt)
                        write_fragment = True

                # # tally spanning reference pairs
                if svtype == 'DEL' and posB - posA < 2 * fragment.lib.sd:
                    ref_straddle_A = False
                    ref_straddle_B = False
                else:
                    ref_straddle_A = fragment.is_pair_straddle(chromA, posA, ref_ciA,
                                                               chromA, posA, ref_ciA,
                                                               False, True,
                                                               min_aligned,
                                                               fragment.lib)
                    ref_straddle_B = fragment.is_pair_straddle(chromB, posB, ref_ciB,
                                                               chromB, posB, ref_ciB,
                                                               False, True,
                                                               min_aligned,
                                                               fragment.lib)

                if ref_straddle_A or ref_straddle_B:
                    # don't allow the pair to jump the entire variant, except for
                    # length-changing SVs like deletions
                    if not (ref_straddle_A and ref_straddle_B) or svtype == 'DEL':
                        p_conc = fragment.p_concordant(var_length)
                        if p_conc is not None:
                            p_reference = p_conc * prob_mapq(fragment.readA) * prob_mapq(fragment.readB)
                            ref_span += (ref_straddle_A + ref_straddle_B) * p_reference / 2

                            fragment.tag_span(1 - p_conc)
                            write_fragment = True

                # write to BAM if requested
                if alignment_outpath is not None and  write_fragment:
                    for read in fragment.primary_reads + [split.read for split in fragment.split_reads]:
                        out_bam_written_reads = write_alignment(read, out_bam, out_bam_written_reads)

            if debug:
                print '--------------------------'
                print 'ref_span:', ref_span
                print 'alt_span:', alt_span
                print 'ref_seq:', ref_seq
                print 'alt_seq:', alt_seq
                print 'alt_clip:', alt_clip

            # in the absence of evidence for a particular type, ignore the reference
            # support for that type as well
            if (alt_seq + alt_clip) < 0.5 and alt_span >= 1:
                alt_seq = 0
                alt_clip = 0
                ref_seq = 0
            if alt_span < 0.5 and (alt_seq + alt_clip) >= 1:
                alt_span = 0
                ref_span = 0

            if ref_seq + alt_seq + ref_span + alt_span + alt_clip > 0:
                # get bayesian classifier
                if var.info['SVTYPE'] == "DUP": is_dup = True
                else: is_dup = False
                alt_splitters = alt_seq + alt_clip
                QR = int(split_weight * ref_seq) + int(disc_weight * ref_span)
                QA = int(split_weight * alt_splitters) + int(disc_weight * alt_span)
                gt_lplist = bayes_gt(QR, QA, is_dup)
                gt_idx = gt_lplist.index(max(gt_lplist))

                # print log probabilities of homref, het, homalt
                if debug:
                    print gt_lplist

                # set the overall variant QUAL score and sample specific fields
                var.genotype(sample.name).set_format('GL', ','.join(['%.0f' % x for x in gt_lplist]))
                var.genotype(sample.name).set_format('DP', int(ref_seq + alt_seq + alt_clip + ref_span + alt_span))
                var.genotype(sample.name).set_format('RO', int(ref_seq + ref_span))
                var.genotype(sample.name).set_format('AO', int(alt_seq + alt_clip + alt_span))
                var.genotype(sample.name).set_format('QR', QR)
                var.genotype(sample.name).set_format('QA', QA)
                # if detailed:
                var.genotype(sample.name).set_format('RS', int(ref_seq))
                var.genotype(sample.name).set_format('AS', int(alt_seq))
                var.genotype(sample.name).set_format('ASC', int(alt_clip))
                var.genotype(sample.name).set_format('RP', int(ref_span))
                var.genotype(sample.name).set_format('AP', int(alt_span))
                try:
                    var.genotype(sample.name).set_format('AB', '%.2g' % (QA / float(QR + QA)))
                except ZeroDivisionError:
                    var.genotype(sample.name).set_format('AB', '.')


                # assign genotypes
                gt_sum = 0
                for gt in gt_lplist:
                    try:
                        gt_sum += 10**gt
                    except OverflowError:
                        gt_sum += 0
                if gt_sum > 0:
                    gt_sum_log = math.log(gt_sum, 10)
                    sample_qual = abs(-10 * (gt_lplist[0] - gt_sum_log)) # phred-scaled probability site is non-reference in this sample
                    if 1 - (10**gt_lplist[gt_idx] / 10**gt_sum_log) == 0:
                        phred_gq = 200
                    else:
                        phred_gq = abs(-10 * math.log(1 - (10**gt_lplist[gt_idx] / 10**gt_sum_log), 10))
                    var.genotype(sample.name).set_format('GQ', int(phred_gq))
                    var.genotype(sample.name).set_format('SQ', sample_qual)
                    var.qual += sample_qual
                    if gt_idx == 1:
                        var.genotype(sample.name).set_format('GT', '0/1')
                    elif gt_idx == 2:
                        var.genotype(sample.name).set_format('GT', '1/1')
                    elif gt_idx == 0:
                        var.genotype(sample.name).set_format('GT', '0/0')
                else:
                    var.genotype(sample.name).set_format('GQ', '.')
                    var.genotype(sample.name).set_format('SQ', '.')
                    var.genotype(sample.name).set_format('GT', './.')
            else:
                var.genotype(sample.name).set_format('GT', './.')
                var.qual = 0
                var.genotype(sample.name).set_format('GQ', '.')
                var.genotype(sample.name).set_format('SQ', '.')
                var.genotype(sample.name).set_format('GL', '.')
                var.genotype(sample.name).set_format('DP', 0)
                var.genotype(sample.name).set_format('AO', 0)
                var.genotype(sample.name).set_format('RO', 0)
                # if detailed:
                var.genotype(sample.name).set_format('AS', 0)
                var.genotype(sample.name).set_format('ASC', 0)
                var.genotype(sample.name).set_format('RS', 0)
                var.genotype(sample.name).set_format('AP', 0)
                var.genotype(sample.name).set_format('RP', 0)
                var.genotype(sample.name).set_format('QR', 0)
                var.genotype(sample.name).set_format('QA', 0)
                var.genotype(sample.name).set_format('AB', '.')

        # after all samples have been processed, write
        vcf_out.write(var.get_var_string() + '\n')
        if var.info['SVTYPE'] == 'BND':
            var2.qual = var.qual
            var2.active_formats = var.active_formats
            var2.genotype = var.genotype
            vcf_out.write(var2.get_var_string() + '\n')

    # close the files
    vcf_in.close()
    vcf_out.close()
    if alignment_outpath is not None:
        out_bam.close()

    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    if args.split_bam is not None:
        sys.stderr.write('Warning: --split_bam (-S) is deprecated. Ignoring %s.\n' % args.split_bam)

    # call primary function
    sv_genotype(args.bam,
                args.input_vcf,
                args.output_vcf,
                args.min_aligned,
                args.split_weight,
                args.disc_weight,
                args.num_samp,
                args.lib_info_path,
                args.debug,
                args.alignment_outpath,
                args.ref_fasta,
                args.sum_quals,
                args.max_reads)

# --------------------------------------
# command-line/console entrypoint

def cli():
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

# initialize the script
if __name__ == '__main__':
    cli()


# NOTES:
#http://code.activestate.com/recipes/578231-probably-the-fastest-memoization-decorator-in-the-/
#http://kitchingroup.cheme.cmu.edu/blog/2013/06/20/Memoizing-expensive-functions-in-python-and-saving-results/
#https://stackoverflow.com/questions/4566769/can-i-memoize-a-python-generator
#https://docs.python.org/2/library/itertools.html
#https://stackoverflow.com/questions/5234090/how-to-take-the-first-n-items-from-a-generator-or-list-in-python
#https://stackoverflow.com/questions/5105517/deep-copy-of-a-dict-in-python
#https://stackoverflow.com/questions/41029158/creating-a-tuple-out-of-args
# http://pysam.readthedocs.io/en/latest/api.html
#     import pysam
#     import itertools
#     s = pysam.AlignmentFile("./tests/data/NA12878.target_loci.sorted.bam", 'rb')
#     foo = list(itertools.islice(s.fetch('2', 2798415, 2798700), 15))
#     len(foo)
#     foo500 = list(itertools.islice(s.fetch('2', 2798415, 2798700), 500))
#     foo15 = list(itertools.islice(s.fetch('2', 2798415, 2798700), 15))
#     foo500[0].__hash__()
#     foo15[0].__hash__()
#     foo15[0].is_unmapped
#     foo15[0].is_duplicate
