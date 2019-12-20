
# ccfind - Circular Complete sequence FINDer

## description
* `ccfind` is a general tool to detect circular complete sequences with clues of terminal redundancy.
* `ccfind` was originally intended for identification of complete virus genomes from metagenome assembly (see `citation`).
* `ccfind` can be used for any contig/genome, but cutoff values should be carefully considered.

## requirements
* ssearch (in [FASTA v36 package](http://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml))
* blastn (in BLAST+ package)
* prodigal
* Ruby (ver >=2.0)
* GNU parallel (if you need parallel computing)

## algorithm
1. extract the first and last `N` bp (`N` is specified by `--terminal-fragment-size`) for each input sequence. Note that if sequence length is less than `N` x 2, the sequence will be ignored.
2. detect sequence homology between the two fragments by `blastn` (as primary screening).
3. For pairs that passed the screening, compute alignment between them by `ssearch` (Smith-Waterman algorithm).
4. If the alignment does not reach the start of the first part, extend alignments to the start, with regarding the exended region as mismatch.
5. If the alignment does not reach the end of the last part, extend alignments to the end, with regarding the exended region as mismatch.
6. Terminal redundancy (and thus circular completeness) is detected when the alignment meets criteria of `--min-percent-identity` and `--min-aligned-length`.

## usage 
```
### ccfind ver 1.4.0 (2019-12-20) ###

[description]
ccfind - Circular Complete sequence FINDer.

ccfind detects completion of nucleotide sequences with clues of terminal redundancy.
ccfind was originally intended for identification of complete virus genomes from metagenome assembly.
ccfind can be used for any contig/genome, but cutoff values should be carefully considered (see README).

[dependencies]
- ssearch        -- Smith-Waterman alignment, included in FASTA program
                    (http://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml)
- blastn         -- nucleotide blast, included in the BLAST+ program
                    (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- prodigal       -- a gene prediction tool, used to detect a well-arranged start position of a circular complete sequence (do not split ORF).
                    (https://github.com/hyattpd/Prodigal)
- ruby (ver >=2.0)

[usage]
$ ccfind <input fasta> <output dir> [options]

- <output dir> should not exist.

[options]
  (general)
    -h, --help
    -v, --version

  (cutoff values)
    -L, --terminal-fragment-size  [int] (default: 500) -- should be longer than estimated terminal redundancy length. If sequence length is less than 2 x (this value), the sequence will be ignored and listed in "result/intermediate/01.skipped.list".
    -i, --min-percent-identity    [int] (default: 94)
    -l, --min-aligned-length      [int] (default: 50)

  (use GNU parallel)
    --ncpus        [int]        -- number of jobs in parallel

[output files]
result/circ.detected.list    -- list of circular complete sequences.
result/circ.fasta            -- circular sequences in fasta format. Terminal redundancy is preserved (without modification).
result/circ.noTR.fasta       -- circular sequences in fasta format. Terminal redundancy is trimmed (redundant region at the end is removed).
result/circ.noTR.cPerm.fasta -- circular sequences in fasta format. Terminal redundancy is trimmed and circular permutation is perfomred, trying not to split any ORF into two fragments (the last intergenic region detected by prodigal is selected as termini).
```

## citation
If you use results genereted by ccfind in your research, please cite:
```
Environmental Viral Genomes Shed New Light on Virus-Host Interactions in the Ocean.
mSphere 2(2):e00359-16 (2017), doi:10.1128/mSphere.00359-16
Yosuke Nishimura, Hiroyasu Watai, Takashi Honda, Tomoko Mihara, Kimiho Omae, Simon Roux, Romain Blanc-Mathieu, Keigo Yamamoto, Pascal Hingamp, Yoshihiko Sako, Matthew B. Sullivan, Susumu Goto, Hiroyuki Ogata, Takashi Yoshida
```
