# DToL insects
Collection of scripts and data useful for analyses of insect genomes produced by the [Darwin Tree of Life project](https://www.darwintreeoflife.org/).

## Download DToL genomes

`get_genomes.py` is a script to pull down DToL genomes from ncbi using tsv file from the [ena project page](https://www.ebi.ac.uk/ena/browser/text-search?query=darwin%20tree%20of%20life)

### Requirements
`-e, --ena` flag requires a tsv file from ENA with all species genome information (see data folder for example).

`-g, --group` flag requires name of insect groups you wish to download data eg. `Insecta` or `Hemiptera, Hymenoptera, Diptera` or `Lepidoptera` (NOTE if Insecta is used, there is no need to specify individual insect orders).

### Output
The output is the genomes from the desired insect group in fasta format.

If the `-i, --info` flag is used, a summary tsv file called `genome_summary.tsv` is also produced which conatins info on insect Order, species name, shortened DToL species name, genome size and chromosome count (number of scaffolds assigned to chromosomes).

### Full usage

```
usage: get_genomes.py [-h] -e ENA -g GROUP [GROUP ...] [-i]

optional arguments:
  -h, --help            show this help message and exit
  -e ENA, --ena ENA     ena tsv file with genome info
  -g GROUP [GROUP ...], --group GROUP [GROUP ...]
                        list of insect orders to download - eg. Insect,
			Odonata, Ephemeroptera, Coleoptera, Hymenoptera,
			Lepidoptera...
  -i, --info            flag to output information on genome data - eg. genome
      			size, chromosome count
```
