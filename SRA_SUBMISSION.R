load("output-tables/table_supp.RData", verbose = TRUE)
head(table_wstn_wdepth)

list_files <- read.csv("list_of_files.txt", header = FALSE) %>% 
  mutate(READ = case_when(
    grepl("_R1_", V1) ~ "Read 1",
    grepl("_R2_", V1) ~ "Read 2"
  )) %>% 
  mutate(SAMPLEID = str_remove_all(V1, "_L001_R[:digit:]_001.fastq.gz")) %>% 
  separate(SAMPLEID, into = c("num", "GOM", "STN", "NISKIN", "18s", "dna", "REP"), remove = FALSE) %>% 
  mutate(STN_NISKIN = paste(STN, NISKIN, sep = "_")) %>% 
  select(-num, -dna, -`18s`)
  

## Make MIMARKS for Sequence Reach Archive.
mimarks <- list_files %>% 
  left_join(table_wstn_wdepth) %>% 
  filter(STN_NISKIN != "S0_N0") %>% 
  select(-GOM, -REP) %>% 
  mutate(sample_name = paste("GOM23", STN, NISKIN, "18S", sep = "_")) %>% 
  group_by(sample_name, STN, NISKIN, STN_NISKIN, Station, Niskin, Depth, Latitude, Longitude, Date, Time_UTC, TRANSECT, MLD, DCM, O2_MIN, SAL_MAX, WATER_MASS) %>% 
  summarise(SEQUENCES = paste(V1, collapse = ", ")) %>% 
  mutate(lat_lon = paste(Latitude, "N", Longitude, "W", sep = " ")) %>% 
  add_column(
    sample_title = "Northern Gulf microbial ecology survey 2023",
    # bioproject_accession = "",
    organism = "ncbitaxon:1289542",
    env_broad_scale = "marine biome[ENVO:00000447]",
    env_medium = "subtropical[ENVO:01000205]",
    env_local_scale = "coastal seawater[ENVO:00002150]|marine oligotrophic desert[ENVO:01000073]",
    geo_loc_name = "USA: Gulf of Mexico"
  ) %>% 
  select(
    sample_name, organism, collection_date = Date, depth = Depth, starts_with("env_"), geo_loc_name, lat_lon,
    time_utc = Time_UTC, transect = TRANSECT, 
    MLD, DCM, O2_MIN, SAL_MAX, WATER_MASS) %>% 
  drop_na()


# Metadata - SRA upload document.
sra_metadata_gom <- list_files %>% 
  left_join(table_wstn_wdepth) %>% 
  filter(STN_NISKIN != "S0_N0") %>% 
  filter(STN_NISKIN != "S12_N10") %>% 
  select(-GOM, -REP) %>% 
  mutate(sample_name = paste("GOM23", STN, NISKIN, "18S", sep = "_")) %>% 
  group_by(sample_name) %>% 
  summarise(SEQUENCES = paste(V1, collapse = ", ")) %>% 
  mutate(library_ID = sample_name) %>% 
  add_column(
    title = "Northern Gulf microbial ecology survey 2023",
    library_strategy = "AMPLICON",
    library_source = "METAGENOMIC",
    library_selection = "PCR",
    library_layout = "paired",
    platform = "ILLUMINA",
    instrument_model = "NextSeq 2000",
    design_description = "Amplicon tag-sequencing from the Gulf of Mexico using, 18S rRNA primers to target the microeukaryotepopulation.",
    filetype = "fastq")

# write.csv(mimarks, "output-tables/mimarks.csv")
# write.csv(sra_metadata_gom, "output-tables/sra_metadata_gom.csv")
sra_metadata_gom
seq_names <- sra_metadata_gom %>% 
  select(`Sample Name` = sample_name, SEQUENCES)
bioproj <- read_tsv("output-tables/BioSampleObjects.txt") %>% 
  left_join(seq_names) %>% 
  left_join(mimarks %>% select(`Sample Name` = sample_name, STN, NISKIN, STN_NISKIN))
bioproj

table_supp1 <- table_wfeature %>%
  select(-Pressure, -datetime_cst, -starts_with("SUN"), -HR_OF_DAY, -stn) %>%
  pivot_wider(names_from = ENV_VARIABLE, values_from = value) %>% 
  left_join(tmp_seq_stats) %>%
  left_join(bioproj) %>% 
  select(Station, Niskin, STN_NISKIN, Depth, DEPTH_BIN, Latitude, Longitude, 
         Date, Time_UTC, DAY_NIGHT, TRANSECT, DIST_OUTFLOW, STN_CATEGORY_4, CLUSTER,
         OilRig_Count, MLD, SAL_MAX, O2_MIN, DCM, Temperature, Salinity, 
         DIC, NH4, NO2, NO3, PO4, SIL, Oxygen_CTD, TA, pH, WATER_MASS, Prok_count, Total_ASVs, Total_sequences, SEQUENCES)

write.csv(table_supp1, file = "output-tables/supp-table-allsampleinfo.csv")
