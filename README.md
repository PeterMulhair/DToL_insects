# DToL insects
Collection of scripts and data useful for analyses of insect genomes produced by the [Darwin Tree of Life project](https://www.darwintreeoflife.org/).

## Download DToL genomes

`get_genomes.py` is a script to pull down DToL genomes for any desired insect order from ncbi using xml file from the [ncbi project page](https://www.ncbi.nlm.nih.gov/assembly?LinkName=bioproject_assembly_all&from_uid=667893) using the `Send to: File` drop down in the top right of the page and selecting XML from the format list.

### Requirements
`-i, --input` flag requires an xml file from ncbi DToL page with all species genome information (see `data/` folder for example).

`-g, --group` flag requires name of insect groups you wish to download data eg. `Insecta` or `Hemiptera, Hymenoptera, Diptera` or `Lepidoptera` (NOTE if Insecta is used, there is no need to specify individual insect orders). This can also be a name of a file with a list of species of interest (see `data/sp_list.txt` as an example). 

* If the `-a, --annotation` flag is used only genomes for species that have annotation data available will be downloaded - requires csv file from [ensembl rapid release site](https://rapid.ensembl.org/info/about/species.html) (see `data/` folder for example)

### Output
The output is the genomes from the desired insect group in fasta format.

If the `-i, --info` flag is used, a summary tsv file called `genome_summary.tsv` is also produced which conatins info on insect Order, species name, shortened DToL species name, genome size and chromosome count (number of scaffolds assigned to chromosomes).

### Full usage

```
usage: get_genomes.py [-h] -i INPUT -g GROUP [GROUP ...] [-a ANNOTATION] [-d]
                      [-t THREADS]
required arguments:
  -i INPUT, --input INPUT
  -g GROUP [GROUP ...], --group GROUP [GROUP ...]
                        list of insect orders to download - eg. Insect,
			Odonata, Ephemeroptera, Coleoptera, Hymenoptera,
	 		Lepidoptera... Or name of text file with species
			names eg. see data/sp_list.txt

optional arguments:
  -h, --help            show this help message and exit
  -a ANNOTATION, --annotation ANNOTATION
                        flag to only download genomes that have annotation
			data available, requires ensembl csv data
  -d, --info            flag to output information on genome data - eg. genome
                        size, chromosome count
  -t THREADS, --threads THREADS
                        number of threads to use i.e. number of genomes to
			download at once, default is 1
```

## Download DToL gene annotation data

`get_genes.py` is a script to pull down DToL annotation data (proteins, cds, gff, gtf) for all available species. Requires xml file (`assembly_result.xml`) from the [ncbi project page](https://www.ncbi.nlm.nih.gov/bioproject/667893) and csv file (`Species.csv`) from [ensembl rapid release site](https://rapid.ensembl.org/info/about/species.html) (which is obtained using the `Download whole table` tab in the top right corner of the page).

### Requirements
`--genomes` flag requires an xml file from ncbi with all species genome information (see `data/` folder for example).

`--genes` flag requires a csv file from ensembl with all species genome information (see `data/` folder for example).

Optional: `--group` flag takes a text file of species names you wish to download data for (see `data/sp_list.txt` folder for example).

### Output
The output is four directories containing peptides, cds, gff and gtf files for each species.

### Full usage

```
usage: get_genes.py [-h] --genomes GENOMES --genes GENES [--group sp_list.txt]

required arguments:
  --genomes GENOMES  ncbi xml file with genome info
  --genes GENES      ensembl csv file with gene annotation info

optional arguments:
  -h, --help         show this help message and exit
  --group            Name of text file with species names
                     to download eg. see data/sp_list.txt
```
