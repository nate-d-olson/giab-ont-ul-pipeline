import numpy
import pysam
import sys
import tqdm

try:
    from biorpy import r
except ImportError:
    print("biorpy isn't installed; figures won't be generated "
          "(install from https://github.com/nspies/biorpy)")
    r = None

import qc
    




def do_qc(path, pdf_path=None):
    inf = pysam.AlignmentFile(path)
    
    read_lengths = []
    unmapped_count = 0
    unmapped_bases = 0

    for count, read in tqdm.tqdm(enumerate(inf)):
        if read.is_secondary or read.is_supplementary:
            continue
        if read.is_unmapped:
            unmapped_bases += read.query_length
            unmapped_count += 1
            continue

        read_lengths.append(read.query_alignment_length)
        # if count > 2000:
        #     print("ENDING EARLY:", str(read)[:1000])
        #     break
        
    read_lengths = numpy.array(read_lengths)

    print(f"Number of mapped reads: {len(read_lengths):,} (excludes supplementary and secondary alignments)")
    print(f"Number of unmapped reads: {unmapped_count:,}")
    print(f"Number of unmapped bases: {unmapped_bases:,}")

    print(f"N50: {qc.n50(read_lengths):,}")
    cutoffs = [0, 10e3, 25e3, 50e3, 100e3, 250e3, 500e3]
    covs = qc.coverages(read_lengths, cutoffs)
    for cutoff, cov in zip(cutoffs, covs):
        print(f" coverage by reads >= {int(cutoff):>10,}: {cov:.3f}x ({cov/covs[0]:6.1%})")
    
    print("Top read lengths:")
    read_lengths.sort()
    
    for length in read_lengths[:-11:-1]:
        print(f" {length:,}")
    
    if pdf_path and r:
        r.pdf(pdf_path)
        qc.plot_coverages(read_lengths)
        r.devoff()


if __name__ == "__main__":
    do_qc(*sys.argv[1:])
