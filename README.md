
# ccfind - circular complete sequence finder by detecting terminal redundancy

## description
ccfind is originally designed to find circular complete contigs/scaffolds generated by the SPAdes assembler. 
ccfind will be applicable for sequences generated by other assemblers, but not tested.

## requirements
* ssearch36 (included in [FASTA version 36](http://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml))
* Ruby (ver >=2.0)

## usage 
```
### ccfind ver 1.1.1 (2018-02-08) ###

[description]
ccfind - tool to find circular complete sequence by detecting terminal redundancy

ccfind is originally designed to find circular complete contigs/scaffolds generated by the SPAdes assembler.
ccfind may be applicable for sequences generated by other assemblers, but not tested.

[dependencies]
- ssearch36        -- included in FASTA program version 36
                      (http://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml)
- ruby (ver >=2.0)

[usage]
$ ccfind <input fasta> <output dir> [-i <min %identity> -l <min aligned len> -L <terminal fragment size>]

[options]
    -h, --help
    -v, --version
    -i, --min-percent-identity    [int] (default: 94)
    -l, --min-aligned-length      [int] (default: 50)
    -L, --terminal-fragment-size  [int] (default: 500)

[output files]
result/05.circ.detected.out -- list of sequences with TR
result/05.circ.noTR.fasta   -- modified input sequences by removing TR in the end of sequence.

Please note that very short sequences that are less than 2 x terminal fragment size (1,000 nt in case of default value) will not be tested.
These sequences will be logged in the file: result/01.skipped.list
```

## citation
If you use results genereted by ccfind in your research, please cite:
```
Environmental Viral Genomes Shed New Light on Virus-Host Interactions in the Ocean.
mSphere 2(2):e00359-16 (2017), doi:10.1128/mSphere.00359-16
Yosuke Nishimura, Hiroyasu Watai, Takashi Honda, Tomoko Mihara, Kimiho Omae, Simon Roux, Romain Blanc-Mathieu, Keigo Yamamoto, Pascal Hingamp, Yoshihiko Sako, Matthew B. Sullivan, Susumu Goto, Hiroyuki Ogata, Takashi Yoshida
```
