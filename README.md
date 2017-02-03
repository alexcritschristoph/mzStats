# mzStats
Statistical tools for comparing MS2 spectra of complex biological samples

Requirements:
*Python: pyteomics, matplotlib, numpy, scipy, lxml, sklearn, argparse.*

*R: ggplot2, vegan, optparse.*

##External resources

[Dereplicator](http://cab.spbu.ru/software/dereplicator/) - for comparing MS2 to predicted spectra of all known NRPS. To use it on GNPS, [click here](http://gnps.ucsd.edu/ProteoSAFe/static/gnps-theoretical.jsp).

[norine](http://bioinfo.lifl.fr/norine/result.jsp?ID=NOR01048) - a database of known NRPs.

[PRISM](https://github.com/magarveylab/prism-releases) - for predicting compound structures from genomes.

Skinnider MA, Dejong CA, Rees PN, Johnston CW, Li H, Webster ALH, Wyatt MA, Magarvey NA (2015). Genomes to natural products PRediction Informatics for Secondary Metabolomes (PRISM). Nucleic Acids Research, 16, 9645-62. doi: 10.1093/nar/gkv1012

[CFM-ID](http://cfmid.wishartlab.com/predict) - for predicting the fragmentation spectrum for a compound. Super cool! For an example [see this](http://cfmid.wishartlab.com/queries/259302b9cf3cde8a2443f0cc3806db2582ce01d5).

Allen, Felicity, et al. "CFM-ID: a web server for annotation, spectrum prediction and metabolite identification from tandem mass spectra." Nucleic acids research 42.W1 (2014): W94-W99.

[ISDB](http://oolonek.github.io/ISDB/) - a database of in silico predicted fragmentation spectra for Natural Products.

Allard, Pierre-Marie, et al. "Integration of molecular networking and in-silico MS/MS fragmentation for natural products dereplication." Analytical chemistry 88.6 (2016): 3317-3323.

[The Secondary Metabolite Bioinformatics Portal](http://www.secondarymetabolites.org/)

[pep2path](http://pep2path.sourceforge.net/)

## Tutorial
 
#### Viewing your data in MZmine

[MZmine](http://mzmine.github.io/) is a program for viewing and exploring tandem MS-MS data. It has been installed on the lab Mac computer.
To explore your data, open up the program. then click on the menu "Raw data methods" and click "Raw data import". Import your mzXML files here. It will take a few minutes to process the samples. In the left hand portion of your screen, you should see your files listed. Clicking the arrows next to them will show you a list of scans sorted by *scan number* in the file, like this:
![image1](http://i.imgur.com/bgdfnBc.png)

You can double click on each scan to see a view of the scan, like this:
![image2](http://i.imgur.com/vxXN4RU.png)

You can also see a very useful 3d view of mz vs retention time by clicking on an entire file, and then clicking on Visualization, and then selecting 3D visualizer. Hit auto-range and select ok:
![image2](http://i.imgur.com/JedavdE.png)


MZmine is really useful for exploring your scans and looking at an overview of your data.

#### Make your mapping file
See `mapping_example.txt` for an example. The mapping file should have at least 3 columns: file, sample, and the first grouping column. You can add additional columns. I like to edit my mapping files in Excel, then copy and paste them into a text editor and save the file. Make sure that tabs and not spaces separate columns in the mapping file.

#### Preprocess data to generate the compounds table.

```
python preprocess_data.py -i mapping_example.txt -o compounds_table.txt
```

Help: `python preprocess_data.py --help`.

This command should generate the file `compounds_table.txt` in the current working directory. Open it in your text editor or Excel for browsing. The value for each compound in each sample is the *largest MS1 precursor peak* among all scans for that sample which were identified as this compound. 

#### Filter your compounds table. The filter command is: 

```
usage: filter_data.py [-h] -i INPUT -m MAPPING [-c CONTROL] [-z MIN_MZ]
                      [--filter_singletons]

Filters an MZstats compound table by removing compounds found in controls,
small compounds, and singletons.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to input compound table
  -m MAPPING, --mapping MAPPING
                        Path to input mapping file
  -c CONTROL, --control CONTROL
                        Name of control samples in mapping file's third column
                        - subtracts any compound found in these controls
  -z MIN_MZ, --min_mz MIN_MZ
                        Minimum mz to retain in the compound table
  --filter_singletons
```

For example, to remove any compounds found in Media samples (a grouping in column 3 of our mapping file from the example), all singletons, and all compounds smaller than 400 mz:

```
python filter_data.py -i compounds_table.txt -m mapping_example.txt -c Media -z 400 --filter_singletons
```

That command will create the file: `compound_table_filtered.txt` in the local directory. This is the new compound table with just the compounds that pass the filter. If it already exists, it will overwrite it. 

#### Find compounds unique to specific conditions.
To find compounds which are present in one set of samples but not in others (say, to identify "species-specific" compounds, or compounds producced only in interactions), use the command:

```
python unique_compounds.py -i compounds_table_filtered.txt -c VR4 -m mapping_example.txt
```

Where "VR4" is a group in the mapping file's 3rd column. This command will identify compounds found in ALL VR4 samples and *no other samples*. If perhaps you want to identify compounds found in at least 2 VR4 samples, you could run:

```
python unique_compounds.py -i compounds_table_filtered.txt -c VR4 -m mapping_example.txt -n 2
```

Note that you can also run this command on an *unfiltered* compound table, like so:


```
python unique_compounds.py -i compounds_table.txt -c VR4 -m mapping_example.txt -n 2
```

#### Create PCoA's and run statistical tests

This is the only portion that requires R (and the R packages vegan, ggplot2, and optparse). It creates a PCoA of the chemical composition data using two distance metrics: bray-curtis (which treats your values as continuous) or jaccard (which treats them as either a 0 or a 1). Jaccard is recommended. This script also runs a statistical test called [ANOSIM](https://en.wikipedia.org/wiki/Analysis_of_similarities) which tests the null hypothesis that "the similarity between groups is greater or equal to the similarity within the groups". In other words, it tests whether your groups are more similar to themselves (within group) than to the other groups in terms of chemical composition. 

To create Principle Coordinate Analysis plots of your data, run:

```
Rscript statistical_tests.-f compounds_table.txt -m mapping_example.txt -c grouping -g grouping
```

Where `-c` refers to the "category" **column name** that the ANOSIM test will be run on, and `-g` refers to the "grouping" **column name** that the PCoA will be colored by. Note that in this case it is "grouping" - what we've named the third column in the mapping file. Down below is an example where we run the same script except on the "bacillus_vs_all" column. 

You can also do your statistics on a filtered compounds table:

```
Rscript statistical_tests.-f compounds_table_filtered.txt -m mapping_example.txt -c bacillus_vs_all -g bacillus_vs_all
```

If you'd like to change the look or style of the PCoA's feel free to edit this R script! It is very simple and easy to edit if you've used R or ggplot2 before.


#### Example: Finding compounds produced during interactions.

We can combine the above in a practical example: identifying compounds produced in an interaction that are not produced in either of the two individual species. In our dataset we have 3 replicates from a species VR22 (called "Streptomyces_sp_C38" in the mapping file) and 3 replicates from a species VR11 (called "Lysinibacillus_sp_CH19" in the mapping file), along with 2 replicates of their interacting colonies called "VR11-VR22". To identify compounds found in the interaction replicates but not the others, run:

```
python unique_compounds.py -i compounds_table.txt -c VR4 -m mapping_example.txt -n 2
```

And you will find that there are two compounds found in both interaction replicates but not in other samples in the table.

#### How do I see the scan for a compound?

`Preprocess_data.py` created your compounds table and it also created another file called *your_compounds_file_name.txt_spectra.txt* at the same time. If we open this file, the lines will look something like:

```
compound_7|159.013472116+-0.0|../MS/011316_Alex_41.mzXML$485|../MS/011316_Alex_25.mzXML$302,15,208,113,394,671,486,767|../MS/011316_Alex_15.mzXML$44,252,553,353,151,453,797|../MS/011316_Alex_26.mzXML$216,788,691,588|../MS/011316_Alex_21.mzXML$495,403,309,587,680|../MS/011316_Alex_16.mzXML$311,120,216,22,497,405,801,589,684|../MS/011316_Alex_23.mzXML$401,581,765|../MS/011316_Alex_19.mzXML$19,121,497,686
```

This means that compound_7 with that mass was found in the listed files; and the numbers after the "$" symbol for each file are the MS2 scan numbers in that file that correspond to this compound. You should be able to open up these scan numbers in a program like mzMINE.
