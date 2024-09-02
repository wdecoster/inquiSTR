# example
# python ~/repositories/inquiSTR/scripts/multiple_hits.py --min-hits {params.min_hits} --samples {params.typeB_list} {input} | bgzip > {output} 2> {log}"


from argparse import ArgumentParser
import gzip
import os
import sys

def main():
    args = get_args()
    samples = [l.rstrip() for l in open(args.samples) if l]
    sys.stderr.write(f'Parsed file, {len(samples)} samples to search for\n')
    if os.path.splitext(args.outliers)[1] == '.gz':
        handle = gzip.open(args.outliers, 'rt')
    else:
        handle = open(args.outliers, 'r')
    for line in handle:
        carriers = [s for s in line.split('\t')[3].split(',') if s in samples]
        sys.stderr.write(f'Found {len(carriers)} samples in line\n')
        if len(carriers) >= args.min_hits:
            if args.max_others is None or len(line.split('\t')[3].split(',')) - len(carriers) <= args.max_others:
                print(line.rstrip() + f'\t{len(carriers)}', end='\n')


def get_args():
    parser = ArgumentParser(description='Multiple hits')
    parser.add_argument('outliers', help='inquistr outlier file')
    parser.add_argument(
        "-m", "--min-hits", help="Minimal number of samples from list", default=2, type=int)
    parser.add_argument("--max-others", help="Maximal number of other samples", default=None, type=int)
    parser.add_argument(
        "-s", "--samples", help="File with list of samples to find", required=True)
    return parser.parse_args()

if __name__ == '__main__':
    main()