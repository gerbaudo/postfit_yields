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
from pdgRounding import pdgRound

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
        def raw_value_as_str(self):
            "values stored internally"
            return "%s +/- %s (tot) +/- %s (stat)" % (self.value, self.tot_err, self.stat_err)
        @property
        def value_as_str(self):
            return "%s +/- %s (syst) +/- %s (stat)" % (self.value, self.syst_err, self.stat_err)
        @property
        def value_as_float(self):
            return "%.2f +/- %.2f (stat) +/- %.2f (syst)" % (float(self.value), float(self.stat_err), float(self.syst_err))
        @property
        def has_tot_err(self):
            return self.tot_err!=0.0
        @property
        def has_stat_err(self):
            return self.stat_err!=0.0

        def combine(self, other):
            "combine tot err from one with stat err from other"
            same_name = self.name==other.name
            one_tot = (self.has_tot_err) is not (other.has_tot_err)
            one_stat = (self.has_stat_err) is not (other.has_stat_err)
            if not (one_tot and one_stat and same_name):
                reason = ('not one_tot' if not one_tot else
                          'not one_stat' if not one_stat else
                          'different name' if not same_name else
                          'unknown')
                raise StandardError("cannot combine %s and %s: %s ; %s (reason: %s)"
                                    % (self.name, other.name,
                                       self.raw_value_as_str, other.raw_value_as_str, reason))
            # perhaps also check consistent yields?
            self.tot_err = self.tot_err if self.has_tot_err else other.tot_err
            self.stat_err = self.stat_err if self.has_stat_err else other.stat_err
            return self

        @property
        def rounded(self):
            "report rounded val +/- stat +/- syst"
            is_data = not self.has_tot_err and not self.has_stat_err
            if is_data:
                return "%d"%self.value
            else:
                max_err = max([self.syst_err, self.stat_err])
                min_err = min([v for v in [self.syst_err, self.stat_err] if v!=0.0])
                err = min_err # not sure whether we want to round on max or min
                value, _    = pdgRound(self.value,    err)
                stat_err, _ = pdgRound(self.stat_err, err)
                syst_err, _ = pdgRound(self.syst_err, err)
                return "%s +/- %s +/- %s"%(value, stat_err, syst_err)


    data_token = 'Data = '
    fake_token = 'Fakes__ = '
    bkg_token = 'Bkd = '
    sig_token = 'Signal = '
    def get_yield(filename, token='Something = '):
        "parse yield, formatted as 'token' value"
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

    # main idea:
    # data (doesn't depend on fit), bkg (from mu=0 fit), fake (from mu=0 fit), signal (from mu=1 fit)
    # fill once `value_toterr` and once `value_staterr`
    def filenameMuhat(base, channel, l1pt, region):
        return base+'/'+"fit%s_%s%s.txt" % (channel, l1pt, '' if region is 'nojets' else '_jets')
    def filenameMu0(base, channel, l1pt, region):
        return base+'/'+"fit%s_%s%s_condMu0.txt" % (channel, l1pt, '' if region is 'nojets' else '_jets')
    def filenameMu1(base, channel, l1pt, region):
        return base+'/'+"fit%s_%s%s_condMu1.txt" % (channel, l1pt, '' if region is 'nojets' else '_jets')

    data_yields = [Yield(sample='data', region=r, l1ptbin=l, channel=c,
                         value_toterr=(get_yield(filenameMuhat(floatall_dir+'/'+'taue', c, l, r),
                                                 data_token), 0.0))
                   for r in regions for l in l1ptbins for c in channels]

    bkg_tot = [Yield(sample='bkg', region=r, l1ptbin=l, channel=c,
                     value_toterr=get_yield_err(filenameMu0(floatall_dir+'/'+'taue', c, l, r),
                                                bkg_token))
               for r in regions for l in l1ptbins for c in channels]
    bkg_stat = [Yield(sample='bkg', region=r, l1ptbin=l, channel=c,
                      value_staterr=get_yield_err(filenameMu0(floatstat_dir+'/'+'taue', c, l, r),
                                                  bkg_token))
                for r in regions for l in l1ptbins for c in channels]
    fake_tot = [Yield(sample='fake', region=r, l1ptbin=l, channel=c,
                      value_toterr=get_yield_err(filenameMu0(floatall_dir+'/'+'taue', c, l, r),
                                                 fake_token))
                for r in regions for l in l1ptbins for c in channels]
    fake_stat = [Yield(sample='fake', region=r, l1ptbin=l, channel=c,
                       value_staterr=get_yield_err(filenameMu0(floatstat_dir+'/'+'taue', c, l, r),
                                                   fake_token))
                 for r in regions for l in l1ptbins for c in channels]
    sig_tot = [Yield(sample='sig', region=r, l1ptbin=l, channel=c,
                     value_toterr=get_yield_err(filenameMu1(floatall_dir+'/'
                                                            +('taue' if c is 'emu' else 'taumu'),
                                                            c, l, r),
                                                sig_token))
               for r in regions for l in l1ptbins for c in channels]
    sig_stat = [Yield(sample='sig', region=r, l1ptbin=l, channel=c,
                      value_staterr=get_yield_err(filenameMu1(floatstat_dir+'/'
                                                              +('taue' if c is 'emu' else 'taumu'),
                                                              c, l, r),
                                                  sig_token))
                for r in regions for l in l1ptbins for c in channels]

    bkg_yields = [t.combine(s) for t,s in zip(bkg_tot, bkg_stat)]
    fake_yields = [t.combine(s) for t,s in zip(fake_tot, fake_stat)]
    # sig_yields = [t.combine(s) for t,s in zip(sig_tot, sig_stat)]
    # ---> issue to be fixed w/ sig breakdown...report tot
    sig_yields = [t for t,s in zip(sig_tot, sig_stat)]

    for dy in data_yields:
        print dy.name,' : ',dy.raw_value_as_str,' : ',dy.rounded
    for by in bkg_yields:
        print by.name,' : ',by.raw_value_as_str,' : ',by.rounded
    for fy in fake_yields:
        print fy.name,' : ',fy.raw_value_as_str,' : ',fy.rounded
    for sy in sig_yields:
        print sy.name,' : ',sy.raw_value_as_str,' : ',sy.rounded

    print '-'*10
    print 'before rounding'
    print '\t'.join(['region', 'data', 'bkg', 'fake', 'sig'])
    for dy, by, fy, sy in zip(data_yields, bkg_yields, fake_yields, sig_yields):
        print dy.name,'  ',dy.raw_value_as_str,'  ',by.raw_value_as_str,'  ',fy.raw_value_as_str,'  ',sy.raw_value_as_str
    print '-'*10
    print 'after rounding'
    print '\t'.join(['region', 'data', 'bkg', 'fake', 'sig'])
    for dy, by, fy, sy in zip(data_yields, bkg_yields, fake_yields, sig_yields):
        print dy.name,'  ',dy.rounded,'  ',by.rounded,'  ',fy.rounded,'  ',sy.rounded

if __name__=='__main__':
    main()
