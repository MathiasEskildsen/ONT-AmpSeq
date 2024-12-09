log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")
# List packages required for the analysis
packages <- c("ampvis2")
# Check if packages are not installed, if not, then install them
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], repos = "https://cloud.r-project.org")
}
# Load the required libraries
invisible(lapply(packages, library, character.only = TRUE))

# Load the metadata and OTU table into R
metadata <- read.delim(snakemake@params[["metadata"]], header = TRUE, sep = "\t")
# Load the OTU tables
otu_tables <- read.delim(snakemake@input[[1]], sep="\t", check.names=F, row.names=1)
# Check if the sampleID in the metadata are the same in as the metadata sheet
if (!all(metadata$sampleID %in% colnames(otu_tables))) {
  stop("SampleIDs in metadata are not all found in the OTU table headers")
}

# Load the OTU tables and metadat into a ampvis2 object
ampvis_obj <- amp_load(
    otutable = otu_tables,
    metadata = metadata
)

# Create a heatmap of the OTU table
heat <- amp_heatmap(ampvis_obj,
    group_by = "SampleID",
    tax_aggregate = "Genus",
    tax_add = "Phylum",
    tax_show = 25
)

ggsave(snakemake@output[["Heatmap"]], heat, width = 186, units = "mm")

# Create rarefraction curve
rare <- amp_rarecurve(
    ampvis_obj,
    stepsize = min(max(colSums(ampvis_obj$abund)), 150),
    color_by = "SampleID"    
) + xlim(0,25000) + labs(title="Rarefraction curve") + ylab("Number of observed OTUs") + 
theme(
    axis.text.x = element_text(angle = 90, hjust = 1), panel.border = element_rect(fill=NA, color="Black", size = 0.5, linetype = "solid"),  panel.grid.major = element_blank(), panel.grid.minor  = element_blank(), axis.line = element_line(colour = "black")
    )
ggsave(snakemake@output[["Rarefraction"]], rare, width = 186, units = "mm")

# Create ordination plot
ord <- amp_ordinate(
    ampvis_obj,
    type = "PCA",
    filter_species = 0.1,
    sample_color_by = "SampleID",
    sample_shape_by = "Condition",
    sample_colorframe = T
) + theme(
    axis.text.x = element_text(angle = 90, hjust = 1), panel.border = element_rect(fill=NA, color="Black", size = 0.5, linetype = "solid"),  panel.grid.major = element_blank(), panel.grid.minor  = element_blank(), axis.line = element_line(colour = "black")
    )

ggsave(snakemake@output[["Ordination"]], ord ,width = 186, units = "mm")

# Create Boxplot
df <- amp_alphadiv(ampvis_obj, richness = T)

box <- ggplot(df, aes(x=Condition, y=ObservedOTUs))+
  geom_boxplot(aes(fill = Condition)) +
  geom_jitter(aes(fill = Condition), shape = 21, size = 2, width = 0.2) + 
  facet_grid(.~SampleID, space = "free", scales = "free") +
  scale_y_continuous(limits = c(0,1400)) + 
  theme_bw()+ xlab("Sample ID") +ylab("Observed OTUs") + theme(panel.border =element_rect(fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(snakemake@output[["Boxplot"]], box, width = 186, units = "mm")

save.image(snakemake@output[["r_env"]])