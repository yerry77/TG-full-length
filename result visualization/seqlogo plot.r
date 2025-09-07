# Load required libraries
library(ggseqlogo)
library(dplyr)
library(ggplot2)

# Read the input data
# The input file should contain columns for each isoform's upstream and downstream sequences
df <- read.table("path/to/TG_isoform_sequences_logo_info.tsv", header = TRUE, sep = "\t")

# Extract upstream and downstream sequences
upstream_seqs <- df$upstream_seq_5prime
downstream_seqs <- df$downstream_seq_3prime

# Calculate dinucleotide distribution at critical positions
# Upstream: positions 6–7, Downstream: positions 10–11
table(substr(upstream_seqs, 6, 7))
table(substr(downstream_seqs, 10, 11))

# Generate sequence logos for all isoforms
ggseqlogo(upstream_seqs, seq_type = "dna")
ggseqlogo(downstream_seqs, seq_type = "dna")

# Separate isoforms into known and novel
known_df <- df %>% filter(isoform_type == "known")
novel_df <- df %>% filter(isoform_type == "novel")

# Generate sequence logos for each group
ggseqlogo(known_df$upstream_seq_5prime, seq_type = "dna") + ggtitle("Known Isoforms - Upstream")
ggseqlogo(known_df$downstream_seq_3prime, seq_type = "dna") + ggtitle("Known Isoforms - Downstream")
ggseqlogo(novel_df$upstream_seq_5prime, seq_type = "dna") + ggtitle("Novel Isoforms - Upstream")
ggseqlogo(novel_df$downstream_seq_3prime, seq_type = "dna") + ggtitle("Novel Isoforms - Downstream")

# Set output directory (replace with your own path)
output_dir <- "path/to/output_directory"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Define date prefix for output files
date_prefix <- "date_"

# Save sequence logos as PDF
p1 <- ggseqlogo(known_df$upstream_seq_5prime, seq_type = "dna") +
  ggtitle("Known Isoforms - Upstream Sequence")
ggsave(file.path(output_dir, paste0(date_prefix, "Known_Upstream_Sequence.pdf")),
       plot = p1, width = 3.8, height = 2.4, units = "in")

p2 <- ggseqlogo(known_df$downstream_seq_3prime, seq_type = "dna") +
  ggtitle("Known Isoforms - Downstream Sequence")
ggsave(file.path(output_dir, paste0(date_prefix, "Known_Downstream_Sequence.pdf")),
       plot = p2, width = 3.8, height = 2.4, units = "in")

p3 <- ggseqlogo(novel_df$upstream_seq_5prime, seq_type = "dna") +
  ggtitle("Novel Isoforms - Upstream Sequence")
ggsave(file.path(output_dir, paste0(date_prefix, "Novel_Upstream_Sequence.pdf")),
       plot = p3, width = 3.8, height = 2.4, units = "in")

p4 <- ggseqlogo(novel_df$downstream_seq_3prime, seq_type = "dna") +
  ggtitle("Novel Isoforms - Downstream Sequence")
ggsave(file.path(output_dir, paste0(date_prefix, "Novel_Downstream_Sequence.pdf")),
       plot = p4, width = 3.8, height = 2.4, units = "in")

# Identify canonical splice sites (GT-AG) for known and novel isoforms
known_df <- known_df %>%
  mutate(
    upstream_canonical = grepl("GT", upstream_seq_5prime),
    downstream_canonical = grepl("AG", downstream_seq_3prime)
  )

novel_df <- novel_df %>%
  mutate(
    upstream_canonical = grepl("GT", upstream_seq_5prime),
    downstream_canonical = grepl("AG", downstream_seq_3prime)
  )

# Calculate canonical splice site percentages for each group
known_percentages <- known_df %>%
  summarise(
    upstream_percent = mean(upstream_canonical) * 100,
    downstream_percent = mean(downstream_canonical) * 100
  )

novel_percentages <- novel_df %>%
  summarise(
    upstream_percent = mean(upstream_canonical) * 100,
    downstream_percent = mean(downstream_canonical) * 100
  )

# Prepare data for plotting
known_percentages_long <- data.frame(
  splice_site = c("Upstream", "Downstream"),
  percentage = c(known_percentages$upstream_percent, known_percentages$downstream_percent),
  isoform_type = "Known"
)

novel_percentages_long <- data.frame(
  splice_site = c("Upstream", "Downstream"),
  percentage = c(novel_percentages$upstream_percent, novel_percentages$downstream_percent),
  isoform_type = "Novel"
)

percentages_long <- bind_rows(known_percentages_long, novel_percentages_long)
percentages_long$isoform_type <- factor(percentages_long$isoform_type, levels = c("Novel", "Known"))

# Plot canonical splice site proportions
ggplot(percentages_long, aes(y = splice_site, x = percentage, fill = isoform_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3.5) +
  labs(title = "Prevalence of Canonical Splice Site Motifs",
       x = "Percentage",
       y = "Splice Site") +
  theme_minimal() +
  scale_fill_manual(values = c("Known" = "#FDAF91FF", "Novel" = "#4DBBD5FF")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.background = element_blank(),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black")
  )
