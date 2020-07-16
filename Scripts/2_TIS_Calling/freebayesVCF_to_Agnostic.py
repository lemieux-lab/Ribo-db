# -*- coding: utf-8 -*-
import sys, os, traceback
from optparse import OptionParser


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg


def main():
    try:
        usage = "usage: %prog [options] <freebayes_vcf_file>"  # review with extra args
        parser = OptionParser(usage)

        # Define optparse arguments, see https://docs.python.org/2/library/optparse.html
        # option processing -- Please revise accordingly !
        parser.add_option("-s", "--somatic", action="store_true", dest="somatic", default=False,
                          help="Only export somatic SNPs")

        (options, args) = parser.parse_args()

        i = 1
        with open(args[0].replace('.vcf', '.pga'), 'w') as outfile:
            print 'Infile ',args[0]
            print 'Outfile ',args[0].replace('.vcf', '.pga')

            outfile.write('\t'.join(['chromosomeNumber', 'uniqueId', 'start', 'end',
                                     'ref', 'alleles', 'quality', 'caller']) + '\n')
            with open(args[0], 'r') as infile:
                for line in infile:
                    if '#' == line[0]:
                        continue

                    sl = [s for s in line.strip().split('\t')]
                    # print sl
                    alleles = sl[9].split(':')[0]
                    # print alleles
                    if '1' not in alleles:  # 0/0, sample in homozygous for REF allele
                        continue

                    chrom, start, ID, ref, alt, qual = sl[0], int(sl[1]) - 1, sl[2], sl[3], sl[4], sl[5]
                    end = start + len(alt.split(',')[0])
                    if '0' in alleles:
                        alt = '{0},{1}'.format(ref, alt)

                    if 1 < len(ref) or (',' not in alt and 1 < len(alt)):  # only retain SNPs
                        continue

                    if options.somatic and '.' != ID:  # Non annotated SNP, must be somatic !
                        continue

                    if 0.001 > float(qual):
                        qual = '0.001'

                    outfile.write('\t'.join([chrom.strip('chr'), str(i), str(start), str(end),
                                             ref, alt, qual, 'freebayes']) + '\n')

                    i += 1

    except Usage, err:
        print >> sys.stderr, ''
        print >> sys.stderr, sys.argv[0].split('/')[-1] + ': ' + str(err.msg)
        print >> sys.stderr, '\t for help use --help'
        return 2

    except Error, err:
        print >> sys.stderr, ''
        print >> sys.stderr, sys.argv[0].split('/')[-1] + ': ' + str(err.msg)
        print >> sys.stderr, '\t for help use --help'
        return 1

    except Exception:
        print >> sys.stderr, ''
        print >> sys.stderr, traceback.format_exc(2)
        print >> sys.stderr, __doc__
        print >> sys.stderr, 'Last line:\n' + '\t'.join(sl)
        return 1


if __name__ == '__main__':
    sys.exit(main())
