Summarise taxonmic distribution in the EMP dataset
================
true
23 July, 2021

-   [Summarise the taxonmic distribution across different sample types
    in the Earth Microbiome
    Project](#summarise-the-taxonmic-distribution-across-different-sample-types-in-the-earth-microbiome-project)
    -   [Set main parameters and download data tables from the EMP FTP
        site](#set-main-parameters-and-download-data-tables-from-the-emp-ftp-site)
    -   [Generate a phyloseq object from the biom and (modified) map
        files](#generate-a-phyloseq-object-from-the-biom-and-modified-map-files)
    -   [Generate a summary table of taxonomic distributions across
        sample
        types](#generate-a-summary-table-of-taxonomic-distributions-across-sample-types)
    -   [Plot heatmaps](#plot-heatmaps)
    -   [Plot a bar chart](#plot-a-bar-chart)
    -   [Plot a violin and dot chart](#plot-a-violin-and-dot-chart)
-   [Generate heat-tree plots using Metacoder (Foster, Sharpton, and
    Grünwald
    2017)](#generate-heat-tree-plots-using-metacoder3-foster_metacoder_2017)
-   [Colophon](#colophon)
-   [References](#references)

## Summarise the taxonmic distribution across different sample types in the Earth Microbiome Project

The Earth Microbiome Project ([EMP](https://earthmicrobiome.org/)) is
collaborative effort to curate and characterise microbial taxonomy and
function on Earth ([Thompson et al.
2017](#ref-thompson_communal_2017)).  
The aim of this script is to summarise the taxonomic distribution across
different sample types, as defined by the EMP. The script summarises
phyla by default but can also summarise the sequences at different
taxonomy levels. The script uses the 3<sup>rd</sup> classification level
of the [EMP
Ontology](https://earthmicrobiome.org/protocols-and-standards/empo/)
(`empo_3`), but this can also be altered.

### Set main parameters and download data tables from the EMP FTP site

``` r
ftp <- "ftp://ftp.microbio.me/emp/release1/"
biom_file <- "emp_cr_silva_16S_123.subset_10k.biom"
map_file <- "emp_qiime_mapping_subset_10k.tsv"
new_map_file <- "emp_qiime_mapping_subset_10k_noBLNK.csv"
Rank <- "Phylum"
Grouping_var <- "empo_3"
EMPO_levels2keep <- c("Soil (non-saline)", "Plant rhizosphere", "Animal corpus", "Water (saline)", "Sediment (saline)", "Water (non-saline)", "Hypersaline (saline)", "Sediment (non-saline)", "Animal distal gut", "Animal proximal gut", "Aerosol (non-saline)", "Animal secretion", "Animal surface", "Plant corpus", "Surface (non-saline)", "Plant surface", "Surface (saline)") # remove unwanted environments
s_size = 1000 # genomes to sample per community
```

``` r
biom_ftp <- paste0(ftp, "otu_tables/closed_ref_silva/", biom_file)
map_ftp <- paste0(ftp, "mapping_files/", map_file)
download.file(biom_ftp, biom_file)
download.file(map_ftp, map_file)

# fill-in missing values with NA, convert to CSV
system("awk 'BEGIN { FS = OFS = \"\t\" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = \"NA\" }; 1' emp_qiime_mapping_subset_10k.tsv  |  awk 'BEGIN { FS=\",\"; OFS=\";\" } {$1=$1; print}' |  awk 'BEGIN { FS=\"\t\"; OFS=\",\" } {$1=$1; print}' > emp_qiime_mapping_subset_10k_noBLNK.csv")
```

### Generate a phyloseq object from the biom and (modified) map files

``` r
(EMP_otu_tax <- import_biom(biom_file, parseFunction = parse_taxonomy_silva_128, 
                            parallel = TRUE))
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:          [ 126730 taxa and 10000 samples ]:
    ## tax_table()   Taxonomy Table:     [ 126730 taxa by 7 taxonomic ranks ]:
    ## taxa are rows

``` r
read_csv(new_map_file,
                    trim_ws = TRUE) %>% 
  column_to_rownames("#SampleID") %>% 
  sample_data() ->
  EMP_map

EMP <- merge_phyloseq(EMP_otu_tax, EMP_map)

# for debugging and testing
# EMP %>%
#   subset_samples(empo_3 %in% EMPO_levels2keep) %>%
#   prune_samples(sample_names(.)[sample(1:nsamples(.), 100)], .) %>%
#   phyloseq::filter_taxa(., function(x) sum(x) > 0, TRUE) %>%
#   prune_taxa(taxa_names(.)[sample(1:ntaxa(.), 1000)], .) ->
#   EMP_subset

EMP %>% 
  PS_merge_samples(grouping_name = Grouping_var, fun = "mean") ->
  EMP_merge

save(EMP, file = paste0("EMP.RDS"))
save(EMP_merge, file = paste0("EMP_merge.RDS"))
# load("./EMP.RDS")
# load("./EMP_merge.RDS")
```

### Generate a summary table of taxonomic distributions across sample types

``` r
# summarise abundance by environment and phylum
EMP %>%
  tax_glom(taxrank = Rank) %>%                     # agglomerate at 'Rank' level
  psmelt() %>%                                         # Melt to long format
  arrange(Rank) %>%                                  # arrange by 'Rank'
  group_by_at(c(Grouping_var, Rank)) %>% 
  # summarise and create a column with the relative abundance
  summarise(`Abundance (mean)` = mean(Abundance),
            n = n(),  # calculates the sample size per group
            SD = sd(Abundance, na.rm = T)) %>%  # calculates the standard error of each group
  mutate(Freq = (100 * `Abundance (mean)` / sum(`Abundance (mean)`)),
         `Freq SD` = (100 * SD / sum(`Abundance (mean)`))) -> # abundance per empo_3 cat
  taxa_sum_df

# determine the phylum order by decreasing abundance (for plotting)
EMP_merge %>%
  tax_glom(taxrank = Rank) %>%                     # agglomerate at 'Rank' level
    transform_sample_counts(., function(x) x / sum(x) * 100) %>% 
  psmelt() %>%                                        # Melt to long format
  arrange(Rank)  %>%                                  # arrange by 'Rank'
  group_by_at(c("OTU", Rank)) %>% 
  summarise(`Abundance (mean)` = mean(Abundance)) %>%
  arrange(`Abundance (mean)`) ->
  Rank_order

# print summary table
taxa_sum_df %>% 
  filter(`Abundance (mean)` != 0) %>%
  kable(., digits = c( 0, 1)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)
```

<table class="table table-striped table-hover table-condensed table-responsive" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
empo\_3
</th>
<th style="text-align:left;">
Phylum
</th>
<th style="text-align:right;">
Abundance (mean)
</th>
<th style="text-align:right;">
n
</th>
<th style="text-align:right;">
SD
</th>
<th style="text-align:right;">
Freq
</th>
<th style="text-align:right;">
Freq SD
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
811
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
1437
</td>
<td style="text-align:right;">
0.8
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
8865
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
7179
</td>
<td style="text-align:right;">
8.7
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
142
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
87
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
171
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
10186
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
9304
</td>
<td style="text-align:right;">
10.0
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
41
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
110
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
278
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
344
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
17601
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
46093
</td>
<td style="text-align:right;">
17.3
</td>
<td style="text-align:right;">
45
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
139
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
114
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
29439
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
27758
</td>
<td style="text-align:right;">
29.0
</td>
<td style="text-align:right;">
27
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
877
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
1248
</td>
<td style="text-align:right;">
0.9
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
215
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
318
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
138
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
221
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
325
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
31337
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
32672
</td>
<td style="text-align:right;">
30.8
</td>
<td style="text-align:right;">
32
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
184
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
435
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
97
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
179
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
269
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
76
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
192
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
73
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
272
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
78
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
162
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
730
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
1062
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
86
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Aerosol (non-saline)
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
122
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
1136
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
3285
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
9205
</td>
<td style="text-align:right;">
11.3
</td>
<td style="text-align:right;">
32
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
102
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
928
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
6033
</td>
<td style="text-align:right;">
3.2
</td>
<td style="text-align:right;">
21
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
962
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
3178
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
4863
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
6151
</td>
<td style="text-align:right;">
16.7
</td>
<td style="text-align:right;">
21
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
277
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
84
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
18728
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
15935
</td>
<td style="text-align:right;">
64.5
</td>
<td style="text-align:right;">
55
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
79
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
291
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal corpus
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
322
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
405
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
2515
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
9155
</td>
<td style="text-align:right;">
3.8
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
34
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
13581
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
19188
</td>
<td style="text-align:right;">
20.6
</td>
<td style="text-align:right;">
29
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
89
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
113
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
604
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
1015
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
6320
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
102
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
218
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
611
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
138
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
27539
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
31837
</td>
<td style="text-align:right;">
41.8
</td>
<td style="text-align:right;">
48
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
1273
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
5624
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
122
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
368
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
56
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Parvarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
103
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
546
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
17547
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
39480
</td>
<td style="text-align:right;">
26.7
</td>
<td style="text-align:right;">
60
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
RsaHF231
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
323
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
994
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
83
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
718
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
2305
</td>
<td style="text-align:right;">
1.1
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
161
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1906
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
694
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
2031
</td>
<td style="text-align:right;">
1.1
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal distal gut
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
2280
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
7379
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
11
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
19428
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
23899
</td>
<td style="text-align:right;">
28.5
</td>
<td style="text-align:right;">
35
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
951
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
2916
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
266
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
135
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
429
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
32369
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
30522
</td>
<td style="text-align:right;">
47.5
</td>
<td style="text-align:right;">
45
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
271
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
1866
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
187
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
559
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
367
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
1367
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
7614
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
12998
</td>
<td style="text-align:right;">
11.2
</td>
<td style="text-align:right;">
19
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
1673
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
3375
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
41
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
106
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
275
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
668
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
2322
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
5153
</td>
<td style="text-align:right;">
3.4
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal proximal gut
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
354
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
124
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
522
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
5146
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
7378
</td>
<td style="text-align:right;">
10.5
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
134
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
5469
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
7229
</td>
<td style="text-align:right;">
11.2
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
243
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
132
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Crenarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
640
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
2632
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
131
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Dictyoglomi
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
15784
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
16637
</td>
<td style="text-align:right;">
32.3
</td>
<td style="text-align:right;">
34
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
1970
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
3010
</td>
<td style="text-align:right;">
4.0
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
225
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
OC31
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
113
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
16921
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
21567
</td>
<td style="text-align:right;">
34.6
</td>
<td style="text-align:right;">
44
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
RsaHF231
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
53
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
SM1K20
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
453
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
2426
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
5875
</td>
<td style="text-align:right;">
5.0
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
220
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal secretion
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
917
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
136
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
324
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
9877
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
15253
</td>
<td style="text-align:right;">
14.0
</td>
<td style="text-align:right;">
22
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Aenigmarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Ancient Archaeal Group(AAG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
86
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
5947
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
10204
</td>
<td style="text-align:right;">
8.5
</td>
<td style="text-align:right;">
15
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
216
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Crenarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
1479
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
5813
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
109
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
566
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
41
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
15695
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
21419
</td>
<td style="text-align:right;">
22.3
</td>
<td style="text-align:right;">
30
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
527
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
1585
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
194
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Marine Hydrothermal Vent Group(MHVG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
OC31
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
66
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
167
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
34831
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
42565
</td>
<td style="text-align:right;">
49.5
</td>
<td style="text-align:right;">
61
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
RsaHF231
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
SM1K20
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
1176
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
4374
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
34
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Thermodesulfobacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1974
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
203
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
573
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Animal surface
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
987
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
239
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
410
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
1489
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
3558
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
3660
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
4586
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Aenigmarchaeota
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
98
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
319
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
83
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
178
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
8816
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
5918
</td>
<td style="text-align:right;">
12.4
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
1450
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
2570
</td>
<td style="text-align:right;">
2.0
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Crenarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
627
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
1253
</td>
<td style="text-align:right;">
0.9
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
338
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
413
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
23421
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
23500
</td>
<td style="text-align:right;">
32.9
</td>
<td style="text-align:right;">
33
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
62
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
510
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
786
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
374
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
725
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
100
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
104
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
192
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
LD1-PA38
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
243
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
404
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
124
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
244
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
41
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Nanohaloarchaeota
</td>
<td style="text-align:right;">
135
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
233
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
127
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
316
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
975
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
1604
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
25655
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
27096
</td>
<td style="text-align:right;">
36.1
</td>
<td style="text-align:right;">
38
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
455
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
739
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
81
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
189
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
1038
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
1989
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
435
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
628
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
208
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
366
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Hypersaline (saline)
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
128
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
95
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
601
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
3284
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
9843
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
43
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
1306
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
3687
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
168
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
100494
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
61858
</td>
<td style="text-align:right;">
77.5
</td>
<td style="text-align:right;">
48
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
77
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
1442
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
4930
</td>
<td style="text-align:right;">
1.1
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
104
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
393
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
89
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
108
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
22764
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
20513
</td>
<td style="text-align:right;">
17.6
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
122
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
202
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
104
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant corpus
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
15054
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
13169
</td>
<td style="text-align:right;">
14.1
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
4532
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
5551
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Aenigmarchaeota
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
456
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
249
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
10383
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
7540
</td>
<td style="text-align:right;">
9.7
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
69
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Candidate division WS6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
64
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
1307
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
1250
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
4736
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
3491
</td>
<td style="text-align:right;">
4.4
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
4226
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
4549
</td>
<td style="text-align:right;">
4.0
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
43
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Diapherotrites
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Dictyoglomi
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
251
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
201
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
906
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
994
</td>
<td style="text-align:right;">
0.8
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
176
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
147
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
2173
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
3244
</td>
<td style="text-align:right;">
2.0
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
1556
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
2367
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
40
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
36
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
377
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
290
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
43
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Marine Hydrothermal Vent Group(MHVG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
166
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
205
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Miscellaneous Euryarchaeotic Group(MEG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
1628
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
1534
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
OC31
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
56
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Parvarchaeota
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
2599
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
1945
</td>
<td style="text-align:right;">
2.4
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
48525
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
40482
</td>
<td style="text-align:right;">
45.4
</td>
<td style="text-align:right;">
38
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
SM1K20
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
566
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
649
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
543
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
1206
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
64
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
6071
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
3839
</td>
<td style="text-align:right;">
5.7
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant rhizosphere
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
85
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
300
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1582
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
3281
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
6332
</td>
<td style="text-align:right;">
3.9
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
13694
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
15776
</td>
<td style="text-align:right;">
16.4
</td>
<td style="text-align:right;">
19
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
96
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
538
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
6554
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
11644
</td>
<td style="text-align:right;">
7.9
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
67
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
233
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
223
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1028
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
129
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
336
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
52
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
6198
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
6476
</td>
<td style="text-align:right;">
7.4
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
48527
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
35238
</td>
<td style="text-align:right;">
58.2
</td>
<td style="text-align:right;">
42
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
RsaHF231
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
95
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1906
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
4280
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
7052
</td>
<td style="text-align:right;">
5.1
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Plant surface
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
5513
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
4427
</td>
<td style="text-align:right;">
8.9
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
2577
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
4477
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Aenigmarchaeota
</td>
<td style="text-align:right;">
284
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
353
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
1134
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
1224
</td>
<td style="text-align:right;">
1.8
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Ancient Archaeal Group(AAG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
146
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
125
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
209
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
621
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
3768
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
3871
</td>
<td style="text-align:right;">
6.1
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
575
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
528
</td>
<td style="text-align:right;">
0.9
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
943
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
766
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
4979
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
3587
</td>
<td style="text-align:right;">
8.0
</td>
<td style="text-align:right;">
6
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Crenarchaeota
</td>
<td style="text-align:right;">
120
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
174
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
1056
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
2384
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
67
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
363
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
93
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
99
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Diapherotrites
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Dictyoglomi
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
121
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
92
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
2604
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
2594
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
84
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
85
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
590
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
2563
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
933
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
1378
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
34
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
371
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
253
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
LD1-PA38
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
130
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
123
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Marine Hydrothermal Vent Group(MHVG)
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
3493
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
5555
</td>
<td style="text-align:right;">
5.6
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Miscellaneous Euryarchaeotic Group(MEG)
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
1896
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
1700
</td>
<td style="text-align:right;">
3.1
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
OC31
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
69
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Parvarchaeota
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
1524
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
918
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
24543
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
15795
</td>
<td style="text-align:right;">
39.6
</td>
<td style="text-align:right;">
25
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
SM1K20
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
316
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
282
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
530
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
552
</td>
<td style="text-align:right;">
0.9
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
591
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
819
</td>
<td style="text-align:right;">
1.0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1088
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
2295
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
2823
</td>
<td style="text-align:right;">
3.7
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (non-saline)
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
145
</td>
<td style="text-align:right;">
544
</td>
<td style="text-align:right;">
153
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
2849
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
6898
</td>
<td style="text-align:right;">
2.9
</td>
<td style="text-align:right;">
7
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
2110
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
2923
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Aenigmarchaeota
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
86
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
1039
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Aigarchaeota
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
190
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
318
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Ancient Archaeal Group(AAG)
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
61
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
1938
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
28882
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:right;">
29
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
9904
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
12878
</td>
<td style="text-align:right;">
9.9
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
145
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
212
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
41
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Candidate division WS6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
40
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
368
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
748
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
2271
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
2829
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Chrysiogenetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
70
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Crenarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
1174
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
3824
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
251
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
289
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
97
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Diapherotrites
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Dictyoglomi
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
68
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
579
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
2577
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
169
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
1696
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
4052
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
92
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
894
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
649
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
684
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
65
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
115
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
70
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
140
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
107
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
447
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
458
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
LD1-PA38
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
387
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
856
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Marine Hydrothermal Vent Group(MHVG)
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
88
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
131
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
331
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
79
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
189
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Miscellaneous Euryarchaeotic Group(MEG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Nanohaloarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
501
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
638
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
OC31
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
89
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
58
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
113
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Parvarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
77
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
2122
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
2384
</td>
<td style="text-align:right;">
2.1
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
66006
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
49322
</td>
<td style="text-align:right;">
66.3
</td>
<td style="text-align:right;">
50
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
RsaHF231
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
SBYG-2791
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
29
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
SM1K20
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
374
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
771
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
124
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
178
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
90
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
254
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
1904
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
3852
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
43
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
84
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
1082
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
uncultured bacterium
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
2251
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
5036
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Sediment (saline)
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
106
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
474
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
21343
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
22301
</td>
<td style="text-align:right;">
20.6
</td>
<td style="text-align:right;">
22
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
9428
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
12099
</td>
<td style="text-align:right;">
9.1
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Aenigmarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Aigarchaeota
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
88
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
98
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
162
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
323
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
391
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
680
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
11108
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
14660
</td>
<td style="text-align:right;">
10.7
</td>
<td style="text-align:right;">
14
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
64
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
467
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
42
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
152
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Candidate division WS6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
67
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
153
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
365
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
1038
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
2393
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
3026
</td>
<td style="text-align:right;">
2.3
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Crenarchaeota
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
360
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
630
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
2090
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
103
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Diapherotrites
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Dictyoglomi
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
215
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
281
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
193
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
656
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
78
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
299
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
4339
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
12302
</td>
<td style="text-align:right;">
4.2
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
28
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
288
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
2564
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
3689
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
116
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
234
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
474
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
89
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Marine Hydrothermal Vent Group(MHVG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
110
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
1010
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Miscellaneous Euryarchaeotic Group(MEG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
984
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
1565
</td>
<td style="text-align:right;">
0.9
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
OC31
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
61
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
98
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Parvarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
1987
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
2464
</td>
<td style="text-align:right;">
1.9
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
34608
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
28725
</td>
<td style="text-align:right;">
33.4
</td>
<td style="text-align:right;">
28
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
RsaHF231
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
SM1K20
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
64
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
61
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
265
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
34
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
615
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
816
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
2454
</td>
<td style="text-align:right;">
0.8
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
51
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
99
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1908
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
10922
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
12722
</td>
<td style="text-align:right;">
10.5
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
84
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
240
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
1304
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Soil (non-saline)
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
1638
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
4822
</td>
<td style="text-align:right;">
2.8
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
7052
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
9451
</td>
<td style="text-align:right;">
12.2
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Aenigmarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Aigarchaeota
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
286
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
988
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
451
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
739
</td>
<td style="text-align:right;">
0.8
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
3868
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
5950
</td>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
10
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Calescamantes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
128
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Candidate division WS6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
27
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
219
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
332
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
635
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
875
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1584
</td>
<td style="text-align:right;">
1.5
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Chrysiogenetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Crenarchaeota
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
4064
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
5189
</td>
<td style="text-align:right;">
7.0
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
364
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
757
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Diapherotrites
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Dictyoglomi
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
67
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
11291
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
16134
</td>
<td style="text-align:right;">
19.5
</td>
<td style="text-align:right;">
28
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
230
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
706
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
248
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
736
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
62
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
298
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Miscellaneous Euryarchaeotic Group(MEG)
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Nanoarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
222
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1966
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
OC31
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Parvarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
445
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1198
</td>
<td style="text-align:right;">
0.8
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
21389
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
25991
</td>
<td style="text-align:right;">
37.0
</td>
<td style="text-align:right;">
45
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
RsaHF231
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
S2R-29
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
SM1K20
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
56
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
3885
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
7795
</td>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
13
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
62
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
543
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Thermodesulfobacteria
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
140
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
310
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1906
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
749
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
1918
</td>
<td style="text-align:right;">
1.3
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (non-saline)
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
953
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
378
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
1374
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
3427
</td>
<td style="text-align:right;">
0.9
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
5011
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
7693
</td>
<td style="text-align:right;">
3.3
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Aenigmarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Ancient Archaeal Group(AAG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
108
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
963
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
124
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
398
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
169
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
20892
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
17969
</td>
<td style="text-align:right;">
13.7
</td>
<td style="text-align:right;">
12
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
55
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
304
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
46
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
181
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
384
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
820
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
1807
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
22891
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
43067
</td>
<td style="text-align:right;">
15.0
</td>
<td style="text-align:right;">
28
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
105
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
438
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
174
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Diapherotrites
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Dictyoglomi
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
110
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
271
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
71
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
10261
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
23703
</td>
<td style="text-align:right;">
6.7
</td>
<td style="text-align:right;">
16
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
190
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
1620
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
267
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
642
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
104
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
358
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
112
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
287
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
LD1-PA38
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
255
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
724
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Marine Hydrothermal Vent Group(MHVG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
43
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
138
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Miscellaneous Euryarchaeotic Group(MEG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
1321
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
5461
</td>
<td style="text-align:right;">
0.9
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
OC31
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
277
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
829
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
60
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
130
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Parvarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
61
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
3328
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
6570
</td>
<td style="text-align:right;">
2.2
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
80232
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
46373
</td>
<td style="text-align:right;">
52.5
</td>
<td style="text-align:right;">
30
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
RsaHF231
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
40
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
SM1K20
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
101
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
272
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
16
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
24
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
154
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
292
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
847
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
53
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
234
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
3967
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
6640
</td>
<td style="text-align:right;">
2.6
</td>
<td style="text-align:right;">
4
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
47
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Surface (saline)
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
53
</td>
<td style="text-align:right;">
117
</td>
<td style="text-align:right;">
134
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
1590
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
4700
</td>
<td style="text-align:right;">
1.2
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
24041
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
51850
</td>
<td style="text-align:right;">
17.9
</td>
<td style="text-align:right;">
39
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Aenigmarchaeota
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
41
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
279
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
229
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
1351
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
26067
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
50661
</td>
<td style="text-align:right;">
19.4
</td>
<td style="text-align:right;">
38
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
161
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
423
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Candidate division WS6
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
104
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
1488
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
785
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
2952
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
1471
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
3662
</td>
<td style="text-align:right;">
1.1
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
21
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Crenarchaeota
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
14013
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
33077
</td>
<td style="text-align:right;">
10.4
</td>
<td style="text-align:right;">
25
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
38
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
19
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Diapherotrites
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Dictyoglomi
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
128
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
465
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
99
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
496
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
39
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
120
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
1851
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
10266
</td>
<td style="text-align:right;">
1.4
</td>
<td style="text-align:right;">
8
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
188
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
332
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
1244
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
63
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
152
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
89
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
517
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Marine Hydrothermal Vent Group(MHVG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
88
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
241
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Miscellaneous Euryarchaeotic Group(MEG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
195
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
563
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
OC31
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
77
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
220
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
242
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
1059
</td>
<td style="text-align:right;">
0.2
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Parvarchaeota
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
2245
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
3986
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
54265
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
86306
</td>
<td style="text-align:right;">
40.3
</td>
<td style="text-align:right;">
64
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
RsaHF231
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
132
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
SM1K20
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
45
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
161
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
207
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
31
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
308
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
49
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
180
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1908
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
5904
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
12303
</td>
<td style="text-align:right;">
4.4
</td>
<td style="text-align:right;">
9
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
44
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
124
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (non-saline)
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
954
</td>
<td style="text-align:right;">
80
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Acetothermia
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Acidobacteria
</td>
<td style="text-align:right;">
708
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
9258
</td>
<td style="text-align:right;">
0.4
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Actinobacteria
</td>
<td style="text-align:right;">
4280
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
8209
</td>
<td style="text-align:right;">
2.5
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Aenigmarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Aerophobetes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
34
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Aigarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Aminicenantes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
30
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Ancient Archaeal Group(AAG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
aquifer1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Aquificae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Armatimonadetes
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Atribacteria
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
20
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Bacteroidetes
</td>
<td style="text-align:right;">
49447
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
67127
</td>
<td style="text-align:right;">
28.8
</td>
<td style="text-align:right;">
39
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Caldiserica
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Candidate division OP3
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
85
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Candidate division SR1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Chlamydiae
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Chlorobi
</td>
<td style="text-align:right;">
32
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
352
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Chloroflexi
</td>
<td style="text-align:right;">
496
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1814
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
CKC4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Cloacimonetes
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
158
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Crenarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Cyanobacteria
</td>
<td style="text-align:right;">
21919
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
50275
</td>
<td style="text-align:right;">
12.8
</td>
<td style="text-align:right;">
29
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Deferribacteres
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
85
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Deinococcus-Thermus
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Diapherotrites
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Dictyoglomi
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Elusimicrobia
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
11
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Euryarchaeota
</td>
<td style="text-align:right;">
840
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
2879
</td>
<td style="text-align:right;">
0.5
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Fibrobacteres
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
54
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Firmicutes
</td>
<td style="text-align:right;">
210
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
986
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Fusobacteria
</td>
<td style="text-align:right;">
62
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1043
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
GAL08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Gemmatimonadetes
</td>
<td style="text-align:right;">
100
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
379
</td>
<td style="text-align:right;">
0.1
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
GOUTA4
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Gracilibacteria
</td>
<td style="text-align:right;">
26
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
48
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Hyd24-12
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
62
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Hydrogenedentes
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
JL-ETNP-Z39
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Kazan-3B-09
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
13
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Latescibacteria
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
74
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
LCP-89
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
9
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
LD1-PA38
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Lentisphaerae
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
102
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Marine Hydrothermal Vent Group(MHVG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Marinimicrobia (SAR406 clade)
</td>
<td style="text-align:right;">
3425
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
9294
</td>
<td style="text-align:right;">
2.0
</td>
<td style="text-align:right;">
5
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Microgenomates
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Miscellaneous Crenarchaeotic Group
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Miscellaneous Euryarchaeotic Group(MEG)
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Nitrospirae
</td>
<td style="text-align:right;">
57
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
270
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
OC31
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Omnitrophica
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
12
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Parcubacteria
</td>
<td style="text-align:right;">
7
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Parvarchaeota
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
PAUC34f
</td>
<td style="text-align:right;">
52
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
221
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Planctomycetes
</td>
<td style="text-align:right;">
949
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
3828
</td>
<td style="text-align:right;">
0.6
</td>
<td style="text-align:right;">
2
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Proteobacteria
</td>
<td style="text-align:right;">
86484
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
90270
</td>
<td style="text-align:right;">
50.4
</td>
<td style="text-align:right;">
53
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
RsaHF231
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Saccharibacteria
</td>
<td style="text-align:right;">
18
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
191
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
SBYG-2791
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
SHA-109
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
SM1K20
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
SM2F11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Spirochaetae
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
412
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Synergistetes
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
TA06
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
61
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Tenericutes
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
33
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Thaumarchaeota
</td>
<td style="text-align:right;">
1250
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
5228
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Thermotogae
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
TM6
</td>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
22
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
uncultured archaeon
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1364
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Verrucomicrobia
</td>
<td style="text-align:right;">
1230
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
5687
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
3
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
WCHB1-60
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
WD272
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
Water (saline)
</td>
<td style="text-align:left;">
Woesearchaeota (DHVEG-6)
</td>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
682
</td>
<td style="text-align:right;">
25
</td>
<td style="text-align:right;">
0.0
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

``` r
# save summary table
taxa_sum_df %>% 
  filter(`Abundance (mean)` != 0) %>% 
  write.csv(file = paste0("EMP_", Grouping_var,"_", Rank,".csv"))
```

### Plot heatmaps

``` r
# Agglomerate taxa
EMP %>% 
  tax_glom(Rank, NArm = TRUE) %>% 
  transform_sample_counts(function(x) x / sum(x) * 100) ->
  Ps_obj_filt_glom

plot_heatmap(
  Ps_obj_filt_glom,
  # method = "NMDS",
  # distance = "bray",
  trans = log_trans(4),
  sample.order = Grouping_var,
  sample.label = Grouping_var,
  taxa.label = Rank,
  taxa.order = Rank_order$OTU
) + 
  theme_bw(base_size = 18) + 
  theme(axis.text.x = element_text(angle = 45.0, vjust = 1, hjust = 1)) +
  scale_fill_viridis(trans = "log", labels = scales::number_format(accuracy = 0.1), option = "magma", na.value = "black") +
  labs(fill = "Abundance (%)")
```

![](EMP_figures/plot%20full%20heatmap-1.png)<!-- -->

``` r
# Agglomerate taxa
EMP_merge %>% 
  tax_glom(Rank, NArm = TRUE) %>% 
  transform_sample_counts(function(x) x / sum(x) * 100) ->
  Ps_obj_filt_glom

plot_heatmap(
  Ps_obj_filt_glom,
  # method = "NMDS",
  # distance = "bray",
  trans = log_trans(4),
  sample.order = Grouping_var,
  sample.label = Grouping_var,
  taxa.label = Rank,
  taxa.order = Rank_order$OTU
) + 
  theme_bw(base_size = 18) + 
  theme(axis.text.x = element_text(angle = 45.0, vjust = 1, hjust = 1)) +
  scale_fill_viridis(trans = "log", labels = scales::number_format(accuracy = 0.1), option = "magma", na.value = "black") +
  labs(fill = "Abundance (%)")
```

![](EMP_figures/plot%20merged%20heatmap-1.png)<!-- -->

### Plot a bar chart

``` r
EMP_merge %>% 
  tax_glom(., Rank) %>% 
  transform_sample_counts(., function(x) x / sum(x) * 100) %>% 
  psmelt() %>%  
  mutate_if(is.character, as.factor) %>% 
  mutate(Phylum = fct_relevel(pull(., Rank), rev(unique(pull(Rank_order, Rank))))) %>% 
  ggplot(aes(fill = !!sym(Rank), x = !!sym(Grouping_var), y = log10(Abundance))) +
  geom_col() +
  theme_set(theme_bw()) +
  theme(axis.text.x = element_text(angle = 45.0, vjust = 1, hjust = 1),
        legend.position = "bottom") +
  labs(y = "Log10 rel. abundance") +
  facet_wrap(~get(Rank), ncol = 5) +
  theme(legend.position = "none")
```

![](EMP_figures/plot%20bar-1.png)<!-- -->

### Plot a violin and dot chart

``` r
EMP %>% 
  tax_glom(., Rank) %>% 
  transform_sample_counts(., function(x) x / sum(x) * 100) %>% 
  psmelt() %>%  
  mutate_if(is.character, as.factor) %>% 
  mutate(Taxa = fct_relevel(pull(., Rank), rev(unique(pull(Rank_order, Rank))))) %>%
  mutate(Taxa = fct_other(Taxa, drop = pull(Rank_order, Rank)[Rank_order$`Abundance (mean)` < 0.1], other_level = "Rare")) %>% 
  ggplot(aes(x = Taxa, y = Abundance)) +
  geom_violin(aes(group = interaction(Taxa, !!sym(Grouping_var))),
              scale = "width") +
  geom_point2(aes(colour = Taxa), 
              position = position_jitter(width = 0.2),
             alpha = 1/4,
             stroke = 0, 
             size = 2) +
  theme(axis.text = element_text(angle = 45.0), 
        axis.text.x = element_text(vjust = 1, hjust = 1) ) + 
  # scale_fill_locuszoom() +
  # scale_color_manual(values = pal("d3js")[c(3, 4, 2)]) +
  labs(x = NULL, y = "Abundance (%)", colour = Rank) + 
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) + 
  facet_wrap(~get(Grouping_var), ncol = 2)  +
  theme(legend.position = "none")
```

![](EMP_figures/plot%20violin-1.png)<!-- -->

## Generate heat-tree plots using [Metacoder](https://grunwaldlab.github.io/metacoder_documentation/) ([Foster, Sharpton, and Grünwald 2017](#ref-foster_metacoder_2017))

``` r
# Convert phyloseq to taxmap
EMP %>% 
  # transform_sample_counts(function(x) x / sum(x) ) %>% # not really needed
  parse_phyloseq() %>% 
  filter_taxa(., !taxon_ranks %in% c("Species", "Genus", "Family", "Order")) -> # tips will be Class
  EMP_taxmap

# Potentially drop out low reads
# EMP_taxmap$data$otu_table <- zero_low_counts(EMP_taxmap, "otu_table", min_count = 5)
# no_reads <- rowSums(EMP_taxmap$data$otu_table[, EMP_taxmap$data$sample_data$sample_id]) == 0
# EMP_taxmap <- filter_obs(EMP_taxmap, data = "tax_data", ! no_reads, drop_taxa = TRUE)

# Sum counts per group
EMP_taxmap$data$tax_counts <- calc_taxon_abund(EMP_taxmap, 
                                                      data = 'otu_table', 
                                                      cols = EMP_taxmap$data$sample_data$sample_id, 
                                                      groups = pull(EMP_taxmap$data$sample_data, Grouping_var))

# Generate a plot for each group
my_plots <- lapply(sort(unique(get_variable(EMP, Grouping_var))), function(group) {
  # do something with time point data
  set.seed(2020) # Make each plot layout look the same if using another layout
  EMP_taxmap %>% 
    heat_tree(.,
              node_label = taxon_names,
              node_size = n_obs,
              node_color = EMP_taxmap$data$tax_counts[[group]],
              node_size_axis_label = 'OTU count',
              node_color_axis_label = 'Sequence count',
              layout = "automatic",
              node_size_range = c(0.013, 0.04), 
              title = group)
})

wrap_plots(my_plots, ncol = 3)
```

![](EMP_figures/plot%20metacoder-1.png)<!-- -->

## Colophon

``` r
sessioninfo::session_info() %>%
  details::details(
    summary = 'Current session info',
    open    = TRUE
  )
```

<details open>
<summary>
<span title="Click to Expand"> Current session info </span>
</summary>

``` r
─ Session info ─────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.1.0 (2021-05-18)
 os       Ubuntu 18.04.5 LTS          
 system   x86_64, linux-gnu           
 ui       X11                         
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Prague               
 date     2021-07-23                  

─ Packages ─────────────────────────────────────────────────────────────────────────────
 package          * version    date       lib source                                
 ade4               1.7-17     2021-06-17 [1] CRAN (R 4.1.0)                        
 ape                5.5        2021-04-25 [1] CRAN (R 4.0.3)                        
 assertthat         0.2.1      2019-03-21 [1] CRAN (R 4.0.2)                        
 backports          1.2.1      2020-12-09 [1] CRAN (R 4.0.2)                        
 bayestestR         0.10.0     2021-05-31 [1] CRAN (R 4.1.0)                        
 Biobase            2.52.0     2021-05-19 [1] Bioconductor                          
 BiocGenerics       0.38.0     2021-05-19 [1] Bioconductor                          
 biomformat         1.20.0     2021-05-19 [1] Bioconductor                          
 Biostrings         2.60.1     2021-06-06 [1] Bioconductor                          
 bit                4.0.4      2020-08-04 [1] CRAN (R 4.0.2)                        
 bit64              4.0.5      2020-08-30 [1] CRAN (R 4.0.2)                        
 bitops             1.0-7      2021-04-24 [1] CRAN (R 4.0.3)                        
 broom              0.7.8      2021-06-24 [1] CRAN (R 4.1.0)                        
 cellranger         1.1.0      2016-07-27 [1] CRAN (R 4.0.2)                        
 cli                3.0.1      2021-07-17 [1] CRAN (R 4.1.0)                        
 clipr              0.7.1      2020-10-08 [1] CRAN (R 4.0.2)                        
 cluster            2.1.2      2021-04-17 [1] CRAN (R 4.0.3)                        
 coda               0.19-4     2020-09-30 [1] CRAN (R 4.0.2)                        
 codetools          0.2-18     2020-11-04 [1] CRAN (R 4.0.2)                        
 colorspace         2.0-2      2021-06-24 [1] CRAN (R 4.1.0)                        
 crayon             1.4.1      2021-02-08 [1] CRAN (R 4.0.3)                        
 data.table         1.14.0     2021-02-21 [1] CRAN (R 4.0.3)                        
 DBI                1.1.1      2021-01-15 [1] CRAN (R 4.0.3)                        
 dbplyr             2.1.1      2021-04-06 [1] CRAN (R 4.0.3)                        
 desc               1.3.0      2021-03-05 [1] CRAN (R 4.0.3)                        
 details            0.2.1      2020-01-12 [1] CRAN (R 4.0.2)                        
 digest             0.6.27     2020-10-24 [1] CRAN (R 4.0.2)                        
 dplyr            * 1.0.7      2021-06-18 [1] CRAN (R 4.1.0)                        
 effectsize         0.4.5      2021-05-25 [1] CRAN (R 4.1.0)                        
 ellipsis           0.3.2      2021-04-29 [1] CRAN (R 4.0.3)                        
 emmeans            1.6.2-1    2021-07-08 [1] CRAN (R 4.1.0)                        
 estimability       1.3        2018-02-11 [1] CRAN (R 4.0.2)                        
 evaluate           0.14       2019-05-28 [1] CRAN (R 4.0.2)                        
 extrafont        * 0.17       2014-12-08 [1] CRAN (R 4.1.0)                        
 extrafontdb        1.0        2012-06-11 [1] CRAN (R 4.0.2)                        
 fansi              0.5.0      2021-05-25 [1] CRAN (R 4.0.3)                        
 farver             2.1.0      2021-02-28 [1] CRAN (R 4.0.3)                        
 forcats          * 0.5.1      2021-01-27 [1] CRAN (R 4.0.3)                        
 foreach            1.5.1      2020-10-15 [1] CRAN (R 4.0.2)                        
 fs                 1.5.0      2020-07-31 [1] CRAN (R 4.0.2)                        
 generics           0.1.0      2020-10-31 [1] CRAN (R 4.0.2)                        
 GenomeInfoDb       1.28.1     2021-07-01 [1] Bioconductor                          
 GenomeInfoDbData   1.2.6      2021-05-25 [1] Bioconductor                          
 ggfittext          0.9.1      2021-01-30 [1] CRAN (R 4.0.3)                        
 ggplot2          * 3.3.5      2021-06-25 [1] CRAN (R 4.1.0)                        
 ggridges           0.5.3      2021-01-08 [1] CRAN (R 4.0.2)                        
 ggtext           * 0.1.1      2020-12-17 [1] CRAN (R 4.0.2)                        
 glue               1.4.2      2020-08-27 [1] CRAN (R 4.0.2)                        
 gridExtra          2.3        2017-09-09 [1] CRAN (R 4.0.2)                        
 gridtext           0.1.4      2020-12-10 [1] CRAN (R 4.0.2)                        
 gtable             0.3.0      2019-03-25 [1] CRAN (R 4.0.2)                        
 haven              2.4.1      2021-04-23 [1] CRAN (R 4.0.3)                        
 highr              0.9        2021-04-16 [1] CRAN (R 4.0.3)                        
 hms                1.1.0      2021-05-17 [1] CRAN (R 4.0.3)                        
 htmltools          0.5.1.1    2021-01-22 [1] CRAN (R 4.0.3)                        
 httr               1.4.2      2020-07-20 [1] CRAN (R 4.0.2)                        
 igraph             1.2.6      2020-10-06 [1] CRAN (R 4.0.2)                        
 insight            0.14.2     2021-06-22 [1] CRAN (R 4.1.0)                        
 IRanges            2.26.0     2021-05-19 [1] Bioconductor                          
 iterators          1.0.13     2020-10-15 [1] CRAN (R 4.0.2)                        
 jsonlite           1.7.2      2020-12-09 [1] CRAN (R 4.0.2)                        
 kableExtra       * 1.3.4      2021-02-20 [1] CRAN (R 4.0.3)                        
 knitr              1.33       2021-04-24 [1] CRAN (R 4.0.3)                        
 labeling           0.4.2      2020-10-20 [1] CRAN (R 4.0.2)                        
 lattice            0.20-44    2021-05-02 [1] CRAN (R 4.0.3)                        
 lazyeval           0.2.2      2019-03-15 [1] CRAN (R 4.0.2)                        
 lifecycle          1.0.0      2021-02-15 [1] CRAN (R 4.0.3)                        
 lubridate          1.7.10     2021-02-26 [1] CRAN (R 4.0.3)                        
 magrittr           2.0.1      2020-11-17 [1] CRAN (R 4.0.2)                        
 MASS               7.3-54     2021-05-03 [1] CRAN (R 4.0.3)                        
 Matrix             1.3-4      2021-06-01 [1] CRAN (R 4.1.0)                        
 metacoder        * 0.3.5      2021-07-15 [1] Github (grunwaldlab/metacoder@8146849)
 mgcv               1.8-36     2021-06-01 [1] CRAN (R 4.1.0)                        
 modelr             0.1.8      2020-05-19 [1] CRAN (R 4.0.2)                        
 multcomp           1.4-17     2021-04-29 [1] CRAN (R 4.0.3)                        
 multtest           2.48.0     2021-05-19 [1] Bioconductor                          
 munsell            0.5.0      2018-06-12 [1] CRAN (R 4.0.2)                        
 mvtnorm            1.1-2      2021-06-07 [1] CRAN (R 4.1.0)                        
 nlme               3.1-152    2021-02-04 [1] CRAN (R 4.0.3)                        
 parameters         0.14.0     2021-05-29 [1] CRAN (R 4.1.0)                        
 patchwork        * 1.1.1      2020-12-17 [1] CRAN (R 4.0.2)                        
 permute            0.9-5      2019-03-12 [1] CRAN (R 4.0.2)                        
 phyloseq         * 1.27.6     2021-07-15 [1] Github (joey711/phyloseq@dc35470)     
 pillar             1.6.1      2021-05-16 [1] CRAN (R 4.0.3)                        
 pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.0.2)                        
 plyr               1.8.6      2020-03-03 [1] CRAN (R 4.0.2)                        
 png                0.1-7      2013-12-03 [1] CRAN (R 4.0.2)                        
 purrr            * 0.3.4      2020-04-17 [1] CRAN (R 4.0.2)                        
 R6                 2.5.0      2020-10-28 [1] CRAN (R 4.0.2)                        
 Rcpp               1.0.7      2021-07-07 [1] CRAN (R 4.1.0)                        
 RCurl              1.98-1.3   2021-03-16 [1] CRAN (R 4.0.3)                        
 readr            * 2.0.0      2021-07-20 [1] CRAN (R 4.1.0)                        
 readxl             1.3.1      2019-03-13 [1] CRAN (R 4.0.2)                        
 reprex             2.0.0      2021-04-02 [1] CRAN (R 4.0.3)                        
 reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.0.2)                        
 rhdf5              2.36.0     2021-05-19 [1] Bioconductor                          
 rhdf5filters       1.4.0      2021-05-19 [1] Bioconductor                          
 Rhdf5lib           1.14.2     2021-07-06 [1] Bioconductor                          
 rlang              0.4.11     2021-04-30 [1] CRAN (R 4.0.3)                        
 rmarkdown        * 2.9        2021-06-15 [1] CRAN (R 4.1.0)                        
 rprojroot          2.0.2      2020-11-15 [1] CRAN (R 4.0.2)                        
 rstudioapi         0.13       2020-11-12 [1] CRAN (R 4.0.2)                        
 Rttf2pt1           1.3.9      2021-07-22 [1] CRAN (R 4.1.0)                        
 rvest              1.0.0      2021-03-09 [1] CRAN (R 4.0.3)                        
 S4Vectors          0.30.0     2021-05-19 [1] Bioconductor                          
 sandwich           3.0-1      2021-05-18 [1] CRAN (R 4.0.3)                        
 scales           * 1.1.1      2020-05-11 [1] CRAN (R 4.0.2)                        
 see              * 0.6.4      2021-05-29 [1] CRAN (R 4.1.0)                        
 sessioninfo        1.1.1      2018-11-05 [1] CRAN (R 4.0.2)                        
 speedyseq        * 0.5.3.9001 2020-10-27 [1] Github (mikemc/speedyseq@8daed32)     
 stringi            1.7.3      2021-07-16 [1] CRAN (R 4.1.0)                        
 stringr          * 1.4.0      2019-02-10 [1] CRAN (R 4.0.2)                        
 survival           3.2-11     2021-04-26 [1] CRAN (R 4.0.3)                        
 svglite          * 2.0.0      2021-02-20 [1] CRAN (R 4.1.0)                        
 systemfonts        1.0.2      2021-05-11 [1] CRAN (R 4.0.3)                        
 TH.data            1.0-10     2019-01-21 [1] CRAN (R 4.0.2)                        
 tibble           * 3.1.2      2021-05-16 [1] CRAN (R 4.0.3)                        
 tidyr            * 1.1.3      2021-03-03 [1] CRAN (R 4.0.3)                        
 tidyselect         1.1.1      2021-04-30 [1] CRAN (R 4.0.3)                        
 tidyverse        * 1.3.1      2021-04-15 [1] CRAN (R 4.0.3)                        
 tzdb               0.1.2      2021-07-20 [1] CRAN (R 4.1.0)                        
 utf8               1.2.1      2021-03-12 [1] CRAN (R 4.0.3)                        
 vctrs              0.3.8      2021-04-29 [1] CRAN (R 4.0.3)                        
 vegan              2.5-7      2020-11-28 [1] CRAN (R 4.0.3)                        
 viridis          * 0.6.1      2021-05-11 [1] CRAN (R 4.0.3)                        
 viridisLite      * 0.4.0      2021-04-13 [1] CRAN (R 4.0.3)                        
 vroom              1.5.3      2021-07-14 [1] CRAN (R 4.1.0)                        
 webshot            0.5.2      2019-11-22 [1] CRAN (R 4.0.2)                        
 withr              2.4.2      2021-04-18 [1] CRAN (R 4.0.3)                        
 xfun               0.24       2021-06-15 [1] CRAN (R 4.1.0)                        
 xml2               1.3.2      2020-04-23 [1] CRAN (R 4.0.2)                        
 xtable             1.8-4      2019-04-21 [1] CRAN (R 4.0.2)                        
 XVector            0.32.0     2021-05-19 [1] Bioconductor                          
 yaml               2.2.1      2020-02-01 [1] CRAN (R 4.0.2)                        
 zlibbioc           1.38.0     2021-05-19 [1] Bioconductor                          
 zoo                1.8-9      2021-03-09 [1] CRAN (R 4.0.3)                        

[1] /home/angel/R/library
[2] /usr/local/lib/R/site-library
[3] /usr/lib/R/site-library
[4] /usr/lib/R/library
```

</details>

<br>

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-foster_metacoder_2017" class="csl-entry">

Foster, Zachary S. L., Thomas J. Sharpton, and Niklaus J. Grünwald.
2017. “Metacoder: An R Package for Visualization and Manipulation of
Community Taxonomic Diversity Data.” *PLOS Computational Biology* 13
(2): e1005404. <https://doi.org/10.1371/journal.pcbi.1005404>.

</div>

<div id="ref-thompson_communal_2017" class="csl-entry">

Thompson, Luke R., Jon G. Sanders, Daniel McDonald, Amnon Amir, Joshua
Ladau, Kenneth J. Locey, Robert J. Prill, et al. 2017. “A Communal
Catalogue Reveals Earth’s Multiscale Microbial Diversity.” *Nature* 551
(7681): 457–63. <https://doi.org/10.1038/nature24621>.

</div>

</div>
