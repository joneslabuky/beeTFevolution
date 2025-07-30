def read_orthos(orthofile):
    ortho_dic = {}
    with open(orthofile, 'r') as reader:
        for line in reader:
            cur_line = line.split()
            cur_og = cur_line[0]
            gene_list = [cur_og] = gene_list
    return ortho_dic

def strat_to_names():
    strat_names = {}
    with open("taxonomy_names.txt", 'r') as reader:
     for line in reader:
         cur_line = line.split().split()
         name = cur_line[1]
         strat_names[strat_id] = name
    return strat_names

def read_strats(stratfile):
    gene_dic = {}
    strat_names = strat_to_names()
    with open(stratfile, 'r') as reader:
        if line.startswith("PS"):
            continue
        cur_line = line.strip().split(";")
        if len(cur_line) < 2:
            continue
        gene = cur_line[1]
        try:
            ps_num = int(cur_line[0])
            gene_dic[gene] = strat_names[ps_num]
        except:
            gene_dic[gene] = "NA"
    return gene_dic

def majority(level_list):
    if len(level_list) < 5:
        retrun "NA"
    level_len = len(level_list)
    for level in set(level_list):
        if level_list.count(level) / level_len > 0.5:
            return level
    return "NA"

#EXCUTING BELOW

#load orthogroups

ortho_dic = read_orthos("/pscratch/bmjo263_uksr/SAV_TFBS/Results_Nov01/Orthogroups/Orthogroups.txt")

#list of species used
species_list = ["AAUR", "ACER", "ADOR", "AFUL", "AMEL", "APIS", "APLU", "APUR", "AROS", "AVIR", "BAFF", "BGER", "BMOR", "BTER", "CCAL", "CFLO", "CGIG", "DMEL", "DNOV",
                "DRER", "EMEX", "FVAR", "HANT", "HLAB", "HLIG", "HQUA", "HROB", "HSAL", "HSAP", "HVOL", "LFIG", "LLAT", "LLEU", "LMAR", "LMOR", "LOEN", "LPAU", "LVIE",
                "LZEP", "MBIC", "MEUR", "MROT", "NFAB", "NFLA", "NMEL", "NVIT", "OBIC", "OBIM", "OLIG", "PDOM", "SINV", "SMON", "SPHA", "STUM", "TDIV", "XVIO"]

#phylostrat maps for each species
species_gene_dics = {}
for species in species_list:
    species_gene_dics[species] = read_strats(f"{species}_animals_final_ps_map.csv")

#write output file
with open("OG_ages.txt", 'w') as outfile:
    outfile.write("OG\tPS\tall_PSs\tgenes\n")
    major_goods = []

    for og, gene_list in ortho_dic.items():
        level_list = []
        used_genes = []

        for gene in gene_list:
            cur_species = gene[0:4]
            if cur_species not in species_gene_dics:
                continue
            if gene in species_gene_dics[cur_species]:
                level_list.append(species_gene_dics[cur_species] [gene])
                used_genes.append(gene)

        major_level = majority(level_list)
        outline = f"{og}\t{major_level}\t{','.join(level_list)}\t{','.join(used_genes)}\n"
        outfile.write(outline)

        if major_level != "NA":
            major_goods.append(major_level)

print(f"{len(major_goods)} orthogroups have a majority phylostrat level assigned.")


#reader = open("/Genomics/kocherlab/berubin/ants/bees/bee_data/renamed_cds_seqs/DNOV_map.txt", 'rU')
#dnov_translate = {}
#for line in reader:
#    cur_line = line.split()
#    dnov_translate[cur_line[1]] = cur_line[0]
#
#ortho_dic = read_orthos("/Genomics/kocherlab/berubin/comparative/halictids/orthology/orthofinder/protein_seqs/OrthoFinder/Results_Jan29_4/Orthogroups/Orthogroups.txt")
#species_gene_dics = {}
#major_goods = []
#outfile = open("OG_ages.txt", 'w')
#outfile.write("OG\tPS\tHighest_expression\tall_PSs\tgenes\n")
#for species in ["AAUR", "APUR", "AVIR", "HLIG", "HQUA", "HRUB", "LLEU", "LMAR", "LFIG", "LZEP", "LVIE", "LPAU", "LOEN", "LMAL", "LCAL", "LALB", "NMEL", "MGEN", "DNOV"]:
#    cur_dic = read_strats("%s_animals_final_ps_map.csv" % species)
#    species_gene_dics[species] = cur_dic

