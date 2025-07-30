# Phylostratigraphy

To determine the approximate evolutionary age of each orthgroup used in our analyses we used the Phylostratigraphy pipline (https://github.com/AlexGa/Phylostratigraphy). We used the proteomes from a total of 56 species, including the 42 bee species in our study along with proteins from 14 species as a reference set spanning a wide diversity of animals as in Jones et al. 2023 (https://github.com/kocherlab/HalictidCompGen). All taxa and their taxonomy are provided in `taxonomy.txt`.

The results of the phylostratigraphic analysis for each species are in the `phylostratigraphic analysis results/` directory. These species files contain the taxonomic level at which each gene originated represented by a number. We used 18 taxonomic levels with family as the most specific level in our analysis:
3 is Bilateria  
4 is Protostomia  
5 is Ecdysozoa  
6 is Panarthropoda  
7 is Arthropoda  
8 is Mandibulata  
9 is Pancrustacea  
10 is Hexapoda  
11 is Insecta  
12 is Dicondylia  
13 is Pterygota  
14 is Neoptera (winged insects)  
15 is Holometabola  
16 is Hymenoptera  
17 is Apocrita  
18 is Aculeata  
19 is Apoidea  
21 is Family (i.e, Andrenidae, Apidae, Colletidae, Halictidae, Megachilidae, Mellitidae)

[SCRIPT!] Was used to assign each OG to one of 4 age range classes: Bilaterian, younger than Bilaterian but older than Hymenoptera, Hymenoptera to Apoidea, or younger than Apoidea (i.e., family, genus, or species specific). To assign age classifications we required at least 5 species be assigned an age by Phylostratigraphy and that the majority of genes with an assigned age for a given OG be assigned the same age. The results are in `OG_ages.txt`.
