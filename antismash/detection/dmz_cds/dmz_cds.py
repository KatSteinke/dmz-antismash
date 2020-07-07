"""Finds a known cluster boundary in an organism by identifying the corresponding genes.
"""
import os

from typing import Dict, List, Union

from antismash.common import fasta, path, secmet, subprocessing


# CLUSTERBLAST DUPLICATES
def get_cds_lengths(record: secmet.Record) -> Dict[str, int]:
    """ Calculates the lengths of each CDS feature in a Record.
        Arguments:
            record: the Record to gather CDS features from
        Returns:
            a dictionary mapping CDS accession to length of the CDS
    """
    lengths = {}
    for cds in record.get_cds_features():
        lengths[cds.get_name()] = len(cds.translation)
    return lengths


def make_blastdb(inputfile: str, db_prefix: str) -> subprocessing.RunResult:
    """ Runs makeblastdb on the inputs to create a blast protein database
        makeblastdb will create 3 files with the given prefix and the extensions:
            .pin, .phr, .psq
        Arguments:
            inputfile: the input filename
            db_prefix: the prefix to use for the created database
        Returns:
            a subprocessing.RunResult instance
    """
    command = ["makeblastdb", "-in", inputfile, "-out", db_prefix, "-dbtype", "prot"]
    result = subprocessing.execute(command)
    if not result.successful():
        raise RuntimeError("makeblastdb failed to run: %s -> %s" % (command, result.stderr[-100:]))
    return result


def run_blast(query: str, database: str) -> str:
    """ Runs blastp, comparing the given query with the given database
        An output file will be created, using the name of the query but with the
        extension changed to .out
        Arguments:
            query: the path of query sequence file
            database: the path of the database to compare to
        Returns:
            the name of the created output file
    """
    out_file = query.rsplit(".", 1)[0] + ".out"
    command = ["blastp", "-db", database, "-query", query, "-outfmt", "6",
               "-max_target_seqs", "10000", "-evalue", "1e-05",
               "-out", out_file]
    res = subprocessing.execute(command)
    if not res.successful():
        raise RuntimeError("blastp run failed: %s..." % res.stderr[-200:])
    return out_file


def make_db_from_record(record: secmet.Record, database_name: str):
    """ Creates a blast database from a record.

    Arguments:
        record: the antiSMASH Record to turn into a database
        database_name: the name of the database to save the record as

    """
    # extract protein sequences from record and write to file
    output_dir = path.get_full_path(database_name)
    record_fasta = os.path.join(output_dir, "record_fasta_for_splits.fa")
    with open(record_fasta, "w") as outfile:
        outfile.write(fasta.get_fasta_from_record(record))
    make_blastdb(record_fasta, database_name)


def blast_splitter_gene(gene_name: str, database_name: str) -> str:
    """Runs blast on a database generated from a record.
    Arguments:
        gene_name: the filename of the splitter gene
        database_name: the filename of the database to blast against
    Returns:
        a string containing all blastp output
    """
    blast_hits = run_blast(gene_name, database_name)
    blast_out = ""
    with open(blast_hits, "r") as infile:
        for line in infile:     # infile.read() produces a faulty line for mystery reasons
            blast_out += line
    return blast_out


def get_top_splitter_name(blast_string: str, record: secmet.Record, min_seq_coverage: float = 90.0,
                          min_perc_identity: float = 90.0) -> Union[str, bool]:
    """Get the name of the CDS most similar to the splitter gene from blast results
     if it passes coverage and identity thresholds.
    Arguments:
        blast_string: the text of the blast search result in tab-separated format
        record: the record the gene was blasted against, for sequence lengths
        min_seq_coverage: the exclusive lower bound of sequence coverage for a match
        min_perc_identity: the exclusive lower bound of identity similarity for a match
    Returns:
        The name of the most similar CDS if one passes the filter, else False.
    """
    filtered_lines = []
    # catch the case of filter not being passed
    top_splitter_name = False
    seqlengths = get_cds_lengths(record)
    blastlines = [line.split("\t") for line in blast_string.rstrip().splitlines()]
    for hit_line in blastlines:
        if len(hit_line) < 12:
            raise ValueError("Malformed blast hit: {}".format("\t".join(hit_line)))
        hit_name = hit_line[1]
        perc_identity = float(hit_line[2])
        align_length = float(hit_line[3])
        if hit_name in seqlengths:
            perc_coverage = (align_length / seqlengths[hit_name]) * 100
        else:
            try:
                seqlength = len(record.get_cds_by_name(hit_name).translation)
            except KeyError:
                raise KeyError("Blast hit {} not found in genome.".format(hit_name))
            perc_coverage = (align_length / seqlength) * 100
        if perc_identity >= min_perc_identity and perc_coverage >= min_seq_coverage:
            filtered_lines.append(hit_line)
    # check that at least one passes
    if filtered_lines:
        # get highest bit score
        filtered_lines.sort(key=lambda x: x[11], reverse=True)
        top_scoring = filtered_lines[0]
        # return the name of the gene
        top_splitter_name = top_scoring[1]
    return top_splitter_name


def find_all_splitters(record: secmet.Record, splitter_genes: List[str], out_dir: str, min_seq_coverage: float = 90.0,
                       min_perc_identity: float = 90.0) -> List[Union[secmet.features.CDSFeature, bool]]:
    """Finds CDS features matching genes known to split two regions in a record.
    Arguments:
        record: the record to search
        splitter_genes: a list of filenames of know splitter genes
        out_dir: directory for outputting dmz_cds files
        min_seq_coverage: the exclusive lower bound of sequence coverage for a match
        min_perc_identity: the exclusive lower bound of identity similarity for a match
    Returns:
        for each splitter gene, the most similar CDS if there is one, else False
    """
    # generate db name
    db_name = os.path.join(out_dir, record.name + "_db_for_split")
    make_db_from_record(record, db_name)
    splitter_cds = []
    for splitter in splitter_genes:
        splitter_hits = blast_splitter_gene(splitter, db_name)
        splitter_name = get_top_splitter_name(splitter_hits, record, min_seq_coverage, min_perc_identity)
        top_splitter = record.get_cds_by_name(splitter_name) if splitter_name else False
        splitter_cds.append(top_splitter)
    return splitter_cds



def get_split_locations(splitter_cdses: List[Union[secmet.features.CDSFeature, bool]]) -> List[int]:
    """Finds the location at which the region should be split, either as the
    start of a single gene or the middle of a "block" of splitter genes.
    Arguments:
        splitter_cdses: A list of splitting CDSs identified by find_all_splitters
    Returns:
        The location(s) where the region should be split
    """
    splitter_locations = [splitter.location.start for splitter in splitter_cdses if splitter]
    return splitter_locations


def find_splitter_cds(record: secmet.Record, splitter_genes: List[str], out_dir: str, min_seq_coverage: float = 90.0,
                      min_perc_identity: float = 90.0) -> List[int]:
    """From a record and a list of genes known to split a region in it,
    finds the location at which to split the region.
    Arguments:
        record: the record to search
        splitter_genes: a list of filenames of known splitter genes
        out_dir: directory for outputting dmz_cds files
        min_seq_coverage: the exclusive lower bound of sequence coverage for a match
        min_perc_identity: the exclusive lower bound of identity similarity for a match
    Returns:
            The location(s) where the region should be split
    """
    all_splitters = find_all_splitters(record, splitter_genes, out_dir, min_seq_coverage, min_perc_identity)
    split_locations = get_split_locations(all_splitters)
    return split_locations


def find_default_splitters(record: secmet.Record, base_out_dir: str) -> List[int]:
    """Convenience function for running search on files in data directory rather than as extra option.
    Arguments:
        record: the Record to run on
        base_out_dir: the antiSMASH output directory
    Returns:
        The location(s) where the region should be split as a list of integers
    """
    data_dir = path.get_full_path(__file__, "data")
    dmz_cds_out_dir = os.path.join(base_out_dir, "dmz_cds")
    if not os.path.exists(dmz_cds_out_dir):
        os.mkdir(dmz_cds_out_dir)
    default_splitters = [splitter_file.path for splitter_file in os.scandir(data_dir)
                         if splitter_file.is_file and not splitter_file.name.endswith(".out")]
    default_splitter_location = find_splitter_cds(record, default_splitters, dmz_cds_out_dir)
    return default_splitter_location
