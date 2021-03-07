ChIP-R ("chipper")
==================

ChIP-R uses an adaptation of the rank product statistic to assess the reproducibility of ChIP-seq peaks by incorporating information from multiple ChIP-seq replicates and "fragmenting" peak locations to better combine the information present across the replicates.

Install
-------

- [Python3.x](https://www.python.org/getit/) with the following packages:
- Numpy
- Scipy

Installation for ChIP-R has been made easy with a range of different install methods available. Installing via any of 
these options will handle all dependencies for you as well.

#### Option 1: via Conda

    conda install chip-r

#### Option 2: via PyPi
    
    pip install ChIP-R

#### Option 3: from source

If you want to install from source:

    git clone https://github.com/rhysnewell/ChIP-R.git
    cd ChIP-R
    python3 setup.py install



Usage
-----

ChIP-R requires only a single input type: A set of any number of BED file regions. Typically the output of peak calling from 
ChIP-seq peak calling on transcription factor or histone mark samples. Alternatively, ChIP-R can also be used on 
ATAC-seq peaks to retrieve reproducible peaks across ATAC-seq experiments.


#### Input

The input BED files must follow ENCODE narrowPeak or broadPeak format specifications. Typically, this format is the default
for peak callers such as MACS2. 

#### Peak calling

ChIP-R is compatible with the output peaks for any peak caller as long as the output is in the correct narrowPeak or broadPeak
format. Additionally, there is no need to call peaks with relaxed thresholds when using your chosen peak caller as is the suggested
by IDR.

#### Parameters

ChIP-R is fairly light on parameters that need to be chosen by the user. A couple of options that users may want to play with is
`minentries` and `size`. 

`minentries` determines the number of peak overlaps required to start calling a peak "reproducible". 
The default of 2 typically provides the best results in our benchmarks but there may be a case where a user requires 
ChIP-R to call peaks within a much stricter window.

`size` determines the minimum peak size during peak output. Transcription factors generally want more punctate peaks, and 
so the default value of 20 may be sufficient. However, histone marks may require a much larger value be set for this depending
on how broad you expect the histone mark to be. Generally, if you find ChIP-R produces too many small noisy peaks then this 
value can be increased to filter them out.

Example
------
    $ chipr -i sample1.bed sample2.bed sample3.bed sample4.bed -m 2 -o output_prefix   

In the command line, type in **'chipr -h '** for detailed usage.

    $ chipr -h
    
    usage: chipr [-h] -i INPUT [INPUT ...] [-o OUTPUT] [-m MINENTRIES]
             [--rankmethod RANKMETHOD] [--duphandling DUPHANDLING]
             [--seed RANDOM_SEED] [-a ALPHA]

    Combine multiple ChIP-seq files and return a union of all peak locations and a
    set confident, reproducible peaks as determined by rank product analysis

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                            ChIP-seq input files. These files must be in either
                            narrowPeak, broadPeak, or regionPeak format. Multiple
                            inputs are separeted by a single space
      -o OUTPUT, --output OUTPUT
                            ChIP-seq output filename prefix
      -m MINENTRIES, --minentries MINENTRIES
                            The minimum peaks between replicates required to form
                            an intersection of the peaks Default: 1
      --rankmethod RANKMETHOD
                            The ranking method used to rank peaks within
                            replicates. Options: 'signalvalue', 'pvalue',
                            'qvalue'. Default: pvalue
      --duphandling DUPHANDLING
                            Specifies how to handle entries that are ranked
                            equally within a replicate Can either take the
                            'average' ranks or a 'random' rearrangement of the
                            ordinal ranks Options: 'average', 'random' Default:
                            'average'
      --seed RANDOM_SEED    Specify a seed to be used in conjunction with the
                            'random' option for -duphandling Must be between 0 and
                            1 Default: 0.5
      -a ALPHA, --alpha ALPHA
                            Alpha specifies the user cut-off value for set of
                            reproducible peaks The analysis will still produce
                            results including peaks within the threshold
                            calculated using the binomial method Default: 0.05
      -s SIZE, --size SIZE  Sets the default minimum peak size when peaks are
                            reconnected after fragmentation. Usually the minimum
                            peak size is determined by the size of surrounding
                            peaks, but in the case that there are no surrounding
                            peaks this value will be used Default: 20


#### Entrypoints

I get that the naming convention for ChIP-R kind of sucks, and have to remember the capital letters and the hyphen can
be frustrating as such you can call ChIP-R from three different entrypoints:
    
    # Easiest
    chipr -h
    
    # Lowercase
    chip-r -h
    
    # Hard mode
    ChIP-R -h

I learned my lesson when naming programs and will never do this again. My apologies :P

Output
------

Important result files:

- **prefixname_ALL.bed**: All intersected peaks, ordered from most significant to least (10 columns)
- **prefixname_T2.bed**: The tier 2 intersected peaks, the peaks that fall within the binomial threshold (10 columns)
- **prefixname_T1.bed**: The tier 1 intersected peaks, the peaks that fall within the user defined threshold (10 columns)
- **prefixname_log.txt**: A log containing the number of peaks appearing in each tier.


prefixname.bed file has 10 columns. The output follows the standard peak format for bed files, with the addition of a 10th column that specifies the ranks of the peaks that produced this possible peak. See the toy example below.

|chr |start|end  |name |score |strand  |signalValue |p-value |q-value|
|----|-----|-----|----|------|-----|------|------|------|
|chr1|9118 |10409|T3_peak_87823|	491|	.	|15.000000	| 0.113938|0.712353	|


Tutorial
--------

The following is a short tutorial on how to generate a set of peaks that can be used as input for ChIP-R.
Usually, an experiment begins with the generation of HTS reads which will require mapping against a reference genome.
Here we will be using dummy values for the reference genome and reads. This example also makes use of BWA for the read mapping
and MACS2 for the peak calling. Alternatives to these programs is perfectly acceptable when using ChIP-R. 

#### Read Mapping

##### Replicate 1
```
bwa mem reference_genome.fasta sample_1.1.fastq sample_1.2.fastq | samtools view -h -b -S -F4 | samtools sort > aln_pe_rep1.bam
```

##### Control 1
```
bwa mem reference_genome.fasta control_1.1.fastq control_1.2.fastq | samtools view -h -b -S -F4 | samtools sort > aln_pe_input1.bam
```

##### Replicate 2
```
bwa mem reference_genome.fasta sample_2.1.fastq sample_2.2.fastq | samtools view -h -b -S -F4 | samtools sort > aln_pe_rep2.bam
```

##### Control 2
```
bwa mem reference_genome.fasta control_2.1.fastq control_2.2.fastq | samtools view -h -b -S -F4 | samtools sort > aln_pe_input2.bam
```

This will produce four BAM files, two experimental replicate read mapping files `aln_pe_rep1.bam` & `aln_pe_rep2.bam`
and two control read mapping files used for the differential peak calling `aln_pe_input1.bam` & `aln_pe_input2.bam`

#### Peak Calling

We keep experiment 1 with control 1, and experiment 2 with control 2.

##### Replicate 1

```
macs callpeak -t aln_pe_rep1.bam -c aln_pe_input1.bam -n rep1
```

##### Replicate 2

```
macs callpeak -t aln_pe_rep2.bam -c aln_pe_input2.bam -n rep2
```

This will result in multiple different output files from MACS2 but the files of interest will be
`rep1_macs2_peaks.bed` & `rep2_macs2_peaks.bed`. We can pass these peaks directly to ChIP-R.

#### Reproducibility analysis

You may want to change the `m` values for your specific experiment. Lower `m` values will produce finer peak boundaries
but there will be many fragmented peaks. However, `m` must always be less than or equal to the number of input files.
```commandline
chipr -i rep1_macs_peaks.bed rep2_macs_peak.bed -m 2 -o output
```

ChIP-R will produce three output files including a log file and two BED files. The BED files consist of 
the set of ALL peak fragments and then the set of optimal peak fragments.

Citation
--------

**Preprint available on bioarxiv**
https://www.biorxiv.org/content/10.1101/2020.11.24.396960v1



Contact
-------

Authors: Rhys Newell, Michael Piper, Mikael Boden, Alexandra Essebier

Contact:  rhys.newell(AT)hdr.qut.edu.au
