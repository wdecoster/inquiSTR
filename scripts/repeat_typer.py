from pysam import AlignmentFile
from argparse import ArgumentParser
import mappy as mp
from pyfaidx import Fasta
import tempfile
import os
import sys


class Repeat(object):
    def __init__(self, read, position, length, type, aligner=None, softclip_side=None):
        self.name = read.query_name
        self.length = -length if type == "contraction" else length
        self.type = type
        self.softclip_side = softclip_side
        self.phase = read.get_tag("HP") if read.has_tag("HP") else "Unp"
        self.phaseset = read.get_tag("PS") if read.has_tag("PS") else "Unp"
        self.seq = self.get_seq(read, position, aligner)

    def get_seq(self, read, read_position, aligner):
        if self.type == "contraction":
            return ""
        else:
            seq = read.query_sequence[read_position : read_position + self.length]
        if self.type == "softclip":
            seq, type = self.realign_softclip(seq, read, aligner)
            if seq:
                self.length = len(seq)
                self.type = type
        return seq

    def realign_softclip(self, seq, read, aligner):
        """
        Check if the end of the softclip doesn't align in the region
        Using a small fasta of ~20kb with the simple repeat masked
        If it does align, it is an insertion/spanning alignment or duplex
        """
        try:
            hit = next(aligner.map(seq))
        except StopIteration:  # If no alignment is found the seq remains as is
            return seq, "softclip"
        orig_strand = "-1" if read.is_reverse else "1"
        # indicates softclip is the duplex of the original alignment
        if orig_strand == "1" and hit.strand != 1:
            return None, None
        # softclip is the duplex of the original alignment, but seq is revcomp for rev reads
        if orig_strand == "-1" and hit.strand == -1:
            return None, None
        not_aligned = [seq[: hit.q_st], seq[hit.q_en :]]
        longest_seq = max(not_aligned, key=len)
        return longest_seq, "aligned_softclip"


def main():
    args = get_args()
    chrom, start, end = process_region(args.region, wobble=args.wobble)
    realign_fasta = make_realignment_fasta(args.ref, contig=chrom, start=start, stop=end)
    a = mp.Aligner(realign_fasta, preset="map-ont")
    for read in AlignmentFile(args.bam).fetch(contig=chrom, start=start, stop=end):
        if (
            start < read.reference_start < end
            or start < read.reference_end < end
            or read.mapping_quality <= 10
        ):
            continue
        softclip_side = "left"
        calls = []
        read_position = 0
        reference_position = read.reference_start + 1
        for operation, length in read.cigartuples:
            if operation in [0, 7, 8]:
                read_position += length
                reference_position += length
            elif operation == 3:
                reference_position += length
            elif operation == 2:
                if length >= args.minlen and start < reference_position < end:
                    calls.append(Repeat(read, read_position, length, "contraction"))
                reference_position += length
            elif operation == 4:
                if length >= args.minlen and start < reference_position < end:
                    calls.append(
                        Repeat(
                            read,
                            read_position,
                            length,
                            "softclip",
                            aligner=a,
                            softclip_side=softclip_side,
                        )
                    )
                read_position += length
            elif operation == 1:
                if length >= args.minlen and start < reference_position < end:
                    calls.append(Repeat(read, read_position, length, "insert"))
                read_position += length
            softclip_side = "right"
        calls = [
            c for c in calls if c.seq is not None
        ]  # That would be a duplex read [I don't remember this comment]
        if calls:
            output_calls(
                calls,
                locus=args.locus,
                ignore_contractions=args.ignore_contractions,
            )


def process_region(region, wobble):
    """parse a region string and extend the start and begin by a wobble space
    defined by a wobble fraction relative to the interval length"""
    chrom = region.split(":")[0]
    start = int(region.split(":")[1].split("-")[0])
    end = int(region.rstrip().split(":")[1].split("-")[1])
    wobble_length = int((end - start) * wobble / 2)
    return chrom, start - wobble_length, end + wobble_length


def output_calls(calls, locus, ignore_contractions=False):
    metadata = {
        "locus": locus,
        "reverse": calls[0].read.is_reverse,
        "length": sum([c.length for c in calls]),
        "phase": calls[0].phase,
        "phaseset": calls[0].phaseset,
        "softclip_side": get_softclip_side(calls),
    }

    metadata["type"], metadata["flag"] = resolve_type(calls, metadata["length"])

    if ignore_contractions and metadata["type"] == "contraction":
        return None

    seq = "".join([c.seq for c in calls])
    metadata_string = " ".join([f"{k}={v}" for k, v in metadata.items() if v is not None])
    print(f">{calls[0].name} {metadata_string}\n{seq}")


def get_softclip_side(calls):
    sides = [c.softclip_side for c in calls if c.softclip_side]
    return sides[0] if sides else None


def resolve_type(calls, length):
    types = [c.type for c in calls]
    if "aligned_softclip" in types:
        flag = "aligned_softclip"
        types = ["insert" if t == "aligned_softclip" else t for t in types]
    else:
        flag = None
    if len(set(types)) == 1:
        return types[0], flag
    elif "softclip" in types:
        return "softclip", flag
    elif "insert" in types and "contraction" in types:
        return ("insert", flag) if length > 0 else ("contraction", flag)
    else:
        sys.stderr.write(f"Unexpected types for resolve_type(): {types}\n")
        return "undefined", flag


def make_realignment_fasta(ref, contig, start, stop):
    fa = Fasta(ref)
    upstream = fa[contig][start - 10000 : start].seq
    downstream = fa[contig][stop : stop + 10000].seq
    masked = "N" * (stop - start)
    handle, name = tempfile.mkstemp(suffix=".vcf")
    out = open(name, "w")
    out.write(f">realignment_sequence\n{upstream}{masked}{downstream}\n")
    os.close(handle)
    return name


def get_args():
    parser = ArgumentParser(description="sum insertions")
    parser.add_argument("--bam", help="bam file to get inserted sequence from", required=True)
    parser.add_argument("--region", help="region string to call expansion in", required=True)
    parser.add_argument("--locus", help="name of the locus in --region", required=True)
    parser.add_argument("--ref", help="reference genome", required=True)
    parser.add_argument(
        "--minlen",
        help="minimal length of insertion/deletion operation",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--wobble",
        help="fraction to extend the region intervals",
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "--ignore_contractions",
        help="only look at expansions",
        action="store_true",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
