import sys
import pandas
import pysam

def get_bam_stats(bam_path, stats_path):
    
    aln_lengths = []
    ref_lengths = []

    bases = []
    read_ids = []

    try:
        bam = pysam.AlignmentFile(bam_path)
    except ValueError:
        print(f"Error reading from {bam_path}...")
        raise

    for read in bam:
        if read.is_unmapped or read.is_secondary: continue
        
        read_ids.append(read.query_name)        
        aln_lengths.append(read.query_alignment_length)
        ref_lengths.append(read.reference_length)

        if not read.is_supplementary:
            bases.append(read.query_length)
        else:
            bases.append(None)
            
    result = pandas.DataFrame({"read_id":read_ids,
                               "aln_length":aln_lengths, 
                               "ref_length":ref_lengths, 
                               "bases":bases})
    result = result.groupby("read_id").agg({"aln_length":[sum, max, "count"], 
                                            "ref_length":[sum, max, "count"], 
                                            "bases":sum})
            
    result.columns = pandas.Series(result.columns.tolist()).apply(pandas.Series).sum(axis=1)
    result = result.rename(columns={"aln_lengthcount":"aln_count", 
                                    "basessum":"bases"})
    
    result.to_csv(stats_path, sep="\t", compression="gzip")


def main():
    if len(sys.argv) == 3:
        bam_path = sys.argv[1]
        out_path = sys.argv[2]
    else:
        print("[ERROR] bam and output path must be provided")


    print("Getting read length stats ...")
    get_bam_stats(bam_path, out_path)

if __name__ == '__main__':
    main()