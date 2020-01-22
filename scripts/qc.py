import datetime
import glob
import gzip
import numpy
import os
import pandas
import pysam


#from biorpy import r, iimage, plotting
#from ont_fast5_api.fast5_file import Fast5File


try:
    from biorpy import r
except ImportError:
    r = None

#BASE_PATH = "/scratch/groups/msalit/nanopore"
BASE_PATH = ""

def plot_coverages(lengths, main=""):
    x = numpy.arange(0,0.5e6+1,1e4)
    y = coverages(lengths, x)

    r.plot(x/1e3, y, type="l", xlim=[0,575], lwd=2, cex=1.5,
           xlab="Read length (kb)", ylab="Coverage by reads > length",
           main=main)

    highlight_x = [0, 50e3, 100e3, 250e3, 500e3]
    highlight_y = coverages(lengths, highlight_x)
    r.points(numpy.array(highlight_x)/1e3, highlight_y, pch=20)
    r.text(numpy.array(highlight_x)/1e3, highlight_y, 
           ["{:.2f} ({:.1%})".format(i, i/highlight_y[0]) for i in highlight_y], pos=4)
    r.mtext(f"{lengths.sum()/1e6:,.1f}mb total")


def n50(values, fraction=0.5):
    if len(values) < 5:
        return numpy.nan
    values = values.copy()
    values.sort()
    
    cumsum = numpy.cumsum(values)
    sum_ = cumsum[-1]
    
    i = numpy.where(cumsum>=fraction*sum_)[0][0]
    
    return values[i]
    

def coverages(lengths, cutoffs):
    covs = []
    for cutoff in cutoffs:
        cur_lengths = lengths[lengths >= cutoff]
        covs.append(cur_lengths.sum() / 3.1e9)
    return covs


def analysis_needs_rerun(raw_path, analysis_path):
    if not os.path.exists(analysis_path):
        return True
    if os.path.getmtime(raw_path) > os.path.getmtime(analysis_path):
        return True
    return False


def fastq_paths(flowcell_name, run_name):
    fastq_dir = f"{BASE_PATH}/raw/{run_name}/fastq"
    
    return glob.glob(f"{fastq_dir}/*.fastq.gz")
    
def fastq_iter(inpath):
    try:
        fastq_file = gzip.open(inpath, "rt")
        fastq_file.read(100)
        fastq_file.seek(0)
    except OSError:
        fastq_file = open(inpath)

    while True:
        lines = [fastq_file.readline().strip() for i in range(4)]
        if sum(map(len, lines)) == 0:
            break
        yield lines
def get_date(header_line):
    for field in header_line.strip().split():
        if field.startswith("start_time"):
            value = field.split("=")[1]
            return datetime.datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ")


def parse_fastq_header_line(header_line):
    fields = header_line[1:].strip().split()
    read_id = fields[0]
    
    seq_date = run_id = None
    for field in fields:
        if field.startswith("start_time"):
            value = field.split("=")[1]
            seq_date = datetime.datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ")
        elif field.startswith("runid"):
            run_id = field.split("=")[1]
    
    return read_id, run_id, seq_date
                
        
def get_fastq_stats(path):
    stats_path = path + ".stats.tsv.gz"
    if not analysis_needs_rerun(path, stats_path):
        result = pandas.read_table(stats_path, sep="\t", index_col=0)
        return result
    
    results = []
    for lines in fastq_iter(path):
        read_id, run_id, seq_date = parse_fastq_header_line(lines[0])
        read_length = len(lines[1])
        results.append([read_id, run_id, seq_date, read_length])
        
    result = pandas.DataFrame(results, columns=["read_id", "run_id", "date","length"])
    result.to_csv(stats_path, sep="\t", compression="gzip")
    return result


def bam_path(run_name, ref_name):
    path = f"{BASE_PATH}/raw/{run_name}/aln_{ref_name}/{run_name}.combined.sorted.bam"
    return path
        
def get_bam_stats(bam_path):
    stats_path = bam_path + ".stats.tsv.gz"
    if not analysis_needs_rerun(bam_path, stats_path):
        result = pandas.read_table(stats_path, sep="\t", index_col=0)
        return result
    
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
            
    result = pandas.DataFrame({"read_id":read_ids, "aln_length":aln_lengths, "ref_length":ref_lengths, "bases":bases})
    result = result.groupby("read_id").agg({"aln_length":[sum, max, "count"], "ref_length":[sum, max, "count"], "bases":sum})
            
    result.columns = pandas.Series(result.columns.tolist()).apply(pandas.Series).sum(axis=1)
    result = result.rename(columns={"aln_lengthcount":"aln_count", "basessum":"bases"})
    
    result.to_csv(stats_path, sep="\t", compression="gzip")

    return result
