#!/bin/env python

# Parse the txt files produced by DumpPostFitHistos, and make a yield table
#
# Example usage:
# $ python/yield_table_from_txt.py -f batch/ll_only/2015-12-13-a/ -s  batch/ll_only/2015-12-13-b/
#
# TODO: detailed description
# TODO: pdg rounding
#
# expected inputs:
#   batch/ll_only/
#   |-- 2015-12-11-a      # free fit
#   |   |-- taue
#   |   |   |-- fitemu_l1pt0_condMu0.txt
#   |   |   |-- fitemu_l1pt0_condMu1.txt
#   |   |   `-- ...
#   |   `-- taumu
#   |       `-- ...
#   `-- 2015-12-11-b      # stat only fit
#       |-- taue
#       |   `-- ...
#       `-- taumu
#           `-- ...
# 
#
# davide.gerbaudo@gmail.com
# Dec 2015



import os, sys, time
from optparse import OptionParser
from math import sqrt

def main():
    parser = OptionParser()
    parser.add_option('-f', '--input-float', help='dir with files from free fit')
    parser.add_option('-s', '--input-stat', help='dir with files from --stat-only fit')
    
    (options, args) = parser.parse_args()
    
    if not options.input_float: parser.error('input-float not given')
    if not options.input_stat: parser.error('input-stat not given')
    
    class Yield(object):
        def __init__(self, sample='', region='', l1ptbin='', channel='',
                     value_toterr=None,
                     value_staterr=None):
            self.sample = sample # signal, fake, bkg
            self.region = region  # 'nojets' or 'jets'
            self.l1ptbin = l1ptbin # l1pt[012]
            self.channel = channel # emu or mue
            if value_toterr:
                self.value, self.tot_err = value_toterr
                self.stat_err = 0.0
            elif value_staterr:
                self.value, self.stat_err = value_staterr
                self.tot_err = 0.0
            else:
                raise StandardError(self.name+
                                    " must specify either value and tot_err or value and stat_err"
                                    ", for example as (0.0, 0.0)"
                                    "(got value_toterr=%s, value_staterr=%s"
                                    %(str(value_toterr), str(value_staterr)))
        @property
        def syst_err(self):
            return sqrt(self.tot_err*self.tot_err - self.stat_err*self.stat_err)
        @property
        def name(self):            
            return '_'.join([self.sample, self.region, self.l1ptbin, self.channel])        
        @property
        def value_as_str(self):
            # return "%s +/- %s (tot) +/- %s (stat)" % (self.value, self.tot_err, self.stat_err)
            return "%s +/- %s (tot) +/- %s (syst) +/- %s (stat)" % (self.value, self.tot_err,
                                                                    self.syst_err, self.stat_err)

        def combine(self, other):
            "combine tot err from one with stat err from other"
            miss_tot = self.tot_err==0.0
            miss_stat = self.stat_err==0.0
            same_name = self.name==other.name
            one_tot = (self.tot_err!=0.0) is not (other.tot_err!=0.0)
            one_stat = (self.stat_err!=0.0) is not (other.stat_err!=0.0)
            if not (one_tot and one_stat and same_name):
                reason = ('not one_tot' if not one_tot else
                          'not one_stat' if not one_stat else
                          'different name' if not same_name else
                          'unknown')
                raise StandardError("cannot combine %s and %s: %s ; %s (reason: %s)"
                                    % (self.name, other.name,
                                       self.value_as_str, other.value_as_str, reason))
            # perhaps also check consistent yields?
            print self.name,' miss_stat ',miss_stat,' miss_tot ',miss_tot
            self.tot_err = other.tot_err if miss_tot else self.tot_err
            self.stat_err = other.stat_err if miss_stat else self.stat_err
            return self

    data_token = 'Data = '
    fake_token = 'Fakes__ = '
    bkg_token = 'Bkd = '
    sig_token = 'Signal = '
    def get_yield(filename, token='Something = '):
        with open(filename) as input_file:
            for line in input_file.readlines():
                if token in line:
                    return int(line.split('=')[1])
    def get_yield_err(filename, token='Something = '):
        "parse yield and error, separated by +/-"
        with open(filename) as input_file:
            for line in input_file.readlines():
                if token in line:
                    val_err = line.split('=')[1]
                    val = float(val_err.split('+-')[0])
                    err = float(val_err.split('+-')[1])
                    return val, err
            raise StandardError("cannot find token '%s' in '%s'"%(token, filename))

    floatall_dir = options.input_float
    floatstat_dir = options.input_stat

    regions = ['nojets', 'jets']
    l1ptbins = ['l1pt0', 'l1pt1', 'l1pt2']
    channels = ['emu', 'mue']

    # base_dir = floatall_dir+'/'+'taue'
    def filenameMuhat(base, channel, l1pt, region):
        return base+'/'+"fit%s_%s%s.txt" % (channel, l1pt, '' if region is 'nojets' else '_jets')
    def filenameMu0(base, channel, l1pt, region):
        return base+'/'+"fit%s_%s%s_condMu0.txt" % (channel, l1pt, '' if region is 'nojets' else '_jets')
    def filenameMu1(base, channel, l1pt, region):
        return base+'/'+"fit%s_%s%s_condMu1.txt" % (channel, l1pt, '' if region is 'nojets' else '_jets')

    # fitemu_l1pt0_condMu0.txt
    # batch/ll_only/2015-12-11-b/taue/fitemu_l1pt1_jets.root
    # batch/ll_only/2015-12-11-b/taue/fitemu_l1pt1.root
    data_yields = [Yield(sample='data', region=r, l1ptbin=l, channel=c,
                         value_toterr=(get_yield(filenameMuhat(floatall_dir+'/'+'taue', c, l, r),
                                                 data_token), 0.0))
                   for r in regions
                   for l in l1ptbins
                   for c in channels]

    # base_dir = floatall_dir+'/'+'taue'

    # bkg_yields_t = [Yield(sample='bkg', region=r, l1ptbin=l, channel=c,
    #                       value_toterr=get_yield_err(filenameMu0(floatall_dir+'/'+'taue', c, l, r),
    #                                                  bkg_token))
    #                for r in regions
    #                for l in l1ptbins
    #                for c in channels]
    # bkg_yields_s = [Yield(sample='bkg', region=r, l1ptbin=l, channel=c,
    #                       value_staterr=get_yield_err(filenameMu0(floatstat_dir+'/'+'taue', c, l, r),
    #                                                   bkg_token))
    #                for r in regions
    #                for l in l1ptbins
    #                for c in channels]
    # fake_yields_t = [Yield(sample='fake', region=r, l1ptbin=l, channel=c,
    #                        value_toterr=get_yield_err(filenameMu0(floatall_dir+'/'+'taue', c, l, r),
    #                                                   fake_token))
    #                  for r in regions
    #                  for l in l1ptbins
    #                  for c in channels]
    # fake_yields_s = [Yield(sample='fake', region=r, l1ptbin=l, channel=c,
    #                        value_staterr=get_yield_err(filenameMu0(floatstat_dir+'/'+'taue',
    #                                                                c, l, r),
    #                                                    fake_token))
    #                  for r in regions
    #                  for l in l1ptbins
    #                  for c in channels]
    sig_yields_t = [Yield(sample='sig', region=r, l1ptbin=l, channel=c,
                          value_toterr=get_yield_err(filenameMu1(floatall_dir+'/'
                                                                 +('taue' if c is 'emu' else 'taumu'),
                                                                 c, l, r),
                                                     sig_token))
                   for r in regions
                   for l in l1ptbins
                   for c in channels]
    sig_yields_s = [Yield(sample='sig', region=r, l1ptbin=l, channel=c,
                          value_staterr=get_yield_err(filenameMu1(floatstat_dir+'/'
                                                                  +('taue' if c is 'emu' else 'taumu'),
                                                                  c, l, r),
                                                      sig_token))
                   for r in regions
                   for l in l1ptbins
                   for c in channels]
    # bkg_yields = [t.combine(s) for t,s in zip(bkg_yields_t, bkg_yields_s)]
    # fake_yields = [t.combine(s) for t,s in zip(fake_yields_t, fake_yields_s)]
    sig_yields = [t.combine(s) for t,s in zip(sig_yields_t, sig_yields_s)]
    # for dy in data_yields:
    #     print dy.name,' : ',dy.value_as_str
    # for by in bkg_yields:
    #     print by.name,' : ',by.value_as_str
    # for fy in fake_yields:
    #     print fy.name,' : ',fy.value_as_str
    for sy in sig_yields:
        print sy.name,' : ',sy.value_as_str

    
if __name__=='__main__':
    main()
