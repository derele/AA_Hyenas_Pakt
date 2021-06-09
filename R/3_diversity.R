### for diversity analysis we use RAREFY!!!

P20k <- prune_samples(sample_sums(PMS)>=10000, PMS)

set.seed(123)
P20kR <- rarefy_even_depth(P20k)
P20kRG <- tax_glom(P20kR, "genus")


pdf("Figures/Genera_per_phylum_Bac.pdf")
ggplot(data.frame(tax_table(PBac)), 
       aes(x=phylum, y=..count..)) +
    geom_bar() +
    coord_polar() +
    scale_y_log10("Number of genera per phylum") +
    theme_bw()
dev.off()

pdf("Figures/Genera_per_phylum_Euk.pdf")
ggplot(data.frame(tax_table(PEuk)), 
       aes(x=phylum, y=..count..)) +
    geom_bar() +
    coord_polar() +
    scale_y_log10("Number of genera per phylum") +
    theme_bw()
dev.off()


pdf("Figures/Reads_per_phylum_Bac.pdf")
ggplot(data.frame(psmelt(PBac)), 
       aes(x=phylum, y=..count.., weight=Abundance)) +
    geom_bar() +
    coord_polar() +
    scale_y_log10("Seqeuncing read counts per phylum") +
    theme_bw()
dev.off()

pdf("Figures/Reads_per_phylum_Euk.pdf")
ggplot(data.frame(psmelt(PEuk)), 
       aes(x=phylum, y=..count.., weight=Abundance)) +
    geom_bar() +
    coord_polar() +
    scale_y_log10("Seqeuncing read counts per phylum") +
    theme_bw()
dev.off()


BacRichness <- estimate_richness(subset_taxa(P20kR, superkingdom%in%"Bacteria"),
                                 measures=c("Observed", "Chao1", "Shannon"))
BacRichness <- merge(BacRichness, sample_data(P20kR), by=0)


ggplot(BacRichness, aes(age_sampling, Observed)) +
    geom_point() +
    geom_line(aes(group=hyena_ID)) +
    stat_smooth(method="lm")

ggplot(BacRichness, aes(age_sampling_cat, Observed)) +
    geom_jitter(width=0.01) +
    geom_line(aes(group=hyena_ID)) 

summary(lm(Observed~age_sampling, BacRichness))


as_tibble(BacRichness) %>%
    group_by(hyena_ID) %>%
    summarize(adu_div = mean(Observed[age_sampling_cat%in%"sampled_as_adult"],
                             na.rm=TRUE),
              juv_div = mean(Observed[age_sampling_cat%in%"sampled_as_juvenile"],
                             na.rm=TRUE)) %>% na.omit() -> pObs

    wilcox.test(pObs$adu_div, pObs$juv_div, paired = TRUE, alternative = "two.sided")

### data:  pObs$adu_div and pObs$juv_div
## V = 136.5, p-value = 0.02784



pdf("Figures/Bacterial_ASVrichness_categorical.pdf")
plot_richness(subset_taxa(P20kR, superkingdom%in%"Bacteria"),
              measures=c("Observed", "Chao1", "Shannon"),
              "age_sampling_cat") +
    geom_boxplot()
dev.off()

pdf("Figures/Bacterial_ASVrichness_continuous.pdf")
plot_richness(subset_taxa(P20kR, superkingdom%in%c("Bacteria")),
              measures=c("Observed", "Chao1", "Shannon"), "age_sampling") +
    stat_smooth(method="lm")
dev.off()


pdf("Figures/Bacterial_Genusrichness_categorical.pdf")
plot_richness(subset_taxa(P20kRG, superkingdom%in%c("Bacteria")),
              measures=c("Observed", "Chao1", "Shannon"),
              "age_sampling_cat") +
    geom_boxplot()
dev.off()

pdf("Figures/Bacterial_Genusrichness_continuous.pdf")
plot_richness(subset_taxa(P20kRG, superkingdom%in%c("Bacteria")),
              measures=c("Observed", "Chao1", "Shannon"), "age_sampling") +
    stat_smooth(method="lm")
dev.off()


pdf("Figures/Euk_Genusrichness_Agecategorical.pdf")
plot_richness(subset_taxa(P20kRG, superkingdom%in%c("Eukaryota")),
              measures=c("Observed", "Chao1", "Shannon"),
              "age_sampling_cat") +
    geom_boxplot()
dev.off()

###

pdf("Figures/Euk_ASVrichness_Juv_RankMom.pdf")
plot_richness(subset_taxa(
    prune_samples(!is.na(sample_data(P20kR)$social_rank_genetic_mum), P20kR),
    superkingdom%in%c("Eukaryota")), 
    measures=c("Observed", "Chao1", "Shannon"),
    "social_rank_genetic_mum") + geom_smooth()
dev.off()

pdf("Figures/Euk_Genusrichness_Juv_RankMom.pdf")
plot_richness(subset_taxa(
    prune_samples(!is.na(sample_data(P20kRG)$social_rank_genetic_mum), P20kRG),
    superkingdom%in%c("Eukaryota")), 
    measures=c("Observed", "Chao1", "Shannon"),
    "social_rank_genetic_mum") + geom_smooth()
dev.off()


pdf("Figures/Euk_ASVrichness_Adu_RankID.pdf")
plot_richness(subset_taxa(
    prune_samples(!is.na(sample_data(P20kR)$social_rank_hyena_ID), P20kR),
    superkingdom%in%c("Eukaryota")), 
    measures=c("Observed", "Chao1", "Shannon"),
    "social_rank_hyena_ID") + geom_smooth()
dev.off()


pdf("Figures/Euk_Genusrichness_Adu_RankID.pdf")
plot_richness(subset_taxa(
    prune_samples(!is.na(sample_data(P20kRG)$social_rank_hyena_ID), P20kRG),
    superkingdom%in%c("Eukaryota")), 
    measures=c("Observed", "Chao1", "Shannon"),
    "social_rank_hyena_ID") + geom_smooth()
dev.off()


