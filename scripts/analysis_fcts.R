


# Load required packages
library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)


#' Read and Prune Phylogenetic Tree
#'
#' @description Loads a phylogenetic tree and removes any tips (samples)
#' that are not in the provided list of sample IDs.
#'
#' @param tree_path Path to the Newick tree file.
#' @param sample_ids Vector of sample IDs to retain in the tree.
#'
#' @return A pruned phylogenetic tree object.
read_and_prune_tree <- function(tree_path, sample_ids, type = c("Newick", "Nexus")) {
  # Match the type argument to determine the function to use, defaulting to "Newick"
  type <- match.arg(type, c("Newick", "Nexus"))
  
  # Read the tree based on the specified type
  if (type == "Newick") {
    tree <- read.tree(tree_path)
  } else if (type == "Nexus") {
    tree <- read.nexus(tree_path)
  } else {
    stop("Unsupported tree type. Use 'Newick' or 'Nexus'.")
  }
  
  # Prune the tree
  tree <- drop.tip(tree, setdiff(tree$tip.label, sample_ids))
  
  return(tree)
}



#' Prepare Annotation Data for Phylogenetic Analysis
#'
#' @description Filters samples used in the phylogeny, merges them with cluster
#' data, and flags potential household transmission events based on index and
#' secondary cases present in the same household and cluster.
#'
#' @param verdi_df Data frame containing metadata (sampleid, hhid, index, etc.).
#' @param df_clusters Data frame with sampleid and cluster assignments.
#'
#' @return Annotated data frame including hh_infected and phy_hhtrans_mem flags.

prepare_annotation_data <- function(verdi_df, df_clusters) {
  df <- verdi_df %>%
    subset(used_in_phylogeny == "YES") %>% 
    mutate(hhid = as.factor(hhid), index = as.factor(index)) %>%
    left_join(df_clusters, by = "sampleid") %>%
    #filter(!is.na(cluster)) %>%
    group_by(hhid, cluster) %>%
    mutate(hh_infected = any(index == 1) & any(index == 0)) %>%
    ungroup() %>%
    mutate(phy_hhtrans_mem = ifelse(index == 0 & hh_infected, TRUE, FALSE))
  return(df)
}

#' Generate Annotated Phylogenetic Tree Plot
#'
#' @description Uses ggtree to visualize a phylogenetic tree with tip annotations
#' indicating household membership and infection status.
#'
#' @param tree Pruned phylogenetic tree object.
#' @param annotation_data Data frame with sampleid, hhid, and index status.
#' @param title Title for the plot.
#'
#' @return A ggplot object representing the annotated tree.

make_ggtree_plot <- function(tree, annotation_data, title = "") {
  levels(annotation_data$index) <- c("Index", "Household Member")
  p <- ggtree(tree) %<+% annotation_data +
    geom_tippoint(aes(color = as.factor(hhid), shape = as.factor(index)),
                  size = 3, na.rm = TRUE) +
    scale_shape_manual(values = c("Index" = 16, "Household Member" = 8),
                       name = "INFECTION STATUS", na.translate = FALSE) +
    scale_color_discrete(name = "HOUSEHOLD INDEX", na.translate = FALSE) +
    theme_tree2() +
    guides(col = guide_legend(nrow = 8)) +
    ggtitle(title)
  return(p)
}


#' Generate Summary Plots for Multiple Clustering Cutoffs
#'
#' @description Evaluates how clustering metrics (e.g., number of clusters,
#' infected households) vary across a range of distance thresholds.
#'
#' @param tree_path Path to the phylogenetic tree file.
#' @param verdi_df Metadata with sample and household information.
#'
#' @return A ggplot object showing clustering metrics vs. distance cutoff.
phylo_summary_plots <- function(tree_path, verdi_df, cutoffs, type = "Newick") {
  phylo_df <- verdi_df %>% filter(used_in_phylogeny == "YES")
  phy <- read_and_prune_tree(tree_path, phylo_df$sampleid, type)
  
  patristic_distances <- cophenetic.phylo(phy)
  distances <- patristic_distances[upper.tri(patristic_distances)]
  hist <- hist(distances, breaks = 30, col = "lightblue", main = "Patristic Distance Distribution")
  
  hc <- hclust(as.dist(patristic_distances), method = "average")

  results <- lapply(cutoffs, function(cutoff) {
    clusters <- cutree(hc, h = cutoff)
    df_clusters <- data.frame(sampleid = names(clusters), cluster = clusters)
    ann <- prepare_annotation_data(verdi_df, df_clusters)
    ann <- ann %>% filter(hh_infected == TRUE) %>% select(sampleid, hhid, index) %>% mutate(cutoff = cutoff,
                                                                                            day = round(cutoff*365))
    list(
      cutoff = cutoff,
      num_clusters = length(unique(clusters)),
      num_hhid = length(unique(ann$hhid)),
      num_hhtrans_mem = sum(ann$index == 0)
    )
  })
  
  ann <- lapply(cutoffs, function(cutoff) {
    clusters <- cutree(hc, h = cutoff)
    df_clusters <- data.frame(sampleid = names(clusters), cluster = clusters)
    ann <- prepare_annotation_data(verdi_df, df_clusters)
    ann <- ann %>% filter(hh_infected == TRUE) %>% select(sampleid, hhid, index) %>% mutate(cutoff = cutoff,
                                                                                            day = round(cutoff*365))
    list(
      ann = ann
    )
  })
  
  # Combine all ann data frames into one
  combined_ann <- do.call(rbind, lapply(ann, function(x) {
    x$ann
  }))
  
  plot_data <- do.call(rbind, lapply(results, as.data.frame))
  print(head(plot_data))
  melted_data <- reshape2::melt(plot_data, id.vars = "cutoff", variable.name = "metric", value.name = "value")
  cutoff_plot <- ggplot(melted_data, aes(x = cutoff, y = value, color = metric)) +
    geom_vline(xintercept = c(15/365, 30/365, 60/365), linetype = "dashed", color = c("black", "green", "red")) +
    geom_line() + geom_point() +
    labs(title = "Clustering Metrics vs. Distance Cutoff (15, 30 and 60 days)", x = "Cutoff", y = "Value", color = "Metric")
  
  return(list(
    hhtrans_table = combined_ann,
    hist = hist,
    cutoff_plot = cutoff_plot))
}


#' Hierarchical Clustering Analysis of Phylogenetic Tree
#'
#' @description Performs clustering on the tree's patristic distances using a
#' specified cutoff. Annotates potential household transmissions and generates
#' summary outputs and a ggtree plot.
#'
#' @param tree_file Path to the phylogenetic tree file.
#' @param verdi_df Data frame with sample metadata.
#' @param threshold Distance cutoff for defining clusters.
#' @param title Title for plot and output.
#'
#' @return A list containing the plot, cluster table, infected household table,
#' summary string, and clustering object.

# Summary plot for phylogenetic analysis
hclust_analysis <- function(tree_file, verdi_df, threshold = 8.4e-5, title, type = "Newick") {
  phylo_df <- verdi_df %>% filter(used_in_phylogeny == "YES")
  pruned_tree <- read_and_prune_tree(tree_file, phylo_df$sampleid, type)
  
  patristic_distances <- cophenetic.phylo(pruned_tree)
  hc <- hclust(as.dist(patristic_distances), method = "average")
  clusters <- cutree(hc, h = threshold)
  df_clusters <- data.frame(sampleid = names(clusters), cluster = clusters)
  
  hhtransmission_df <- prepare_annotation_data(verdi_df, df_clusters)
  annotation_data <- hhtransmission_df %>%
    filter(hh_infected == TRUE) %>%
    select(sampleid, hhid, index)
  
  p <- make_ggtree_plot(pruned_tree, annotation_data,
                        title = paste(title, " WORKFLOW (122 samples, cutoff", threshold))
  print(p)
  message(sprintf("Pour un cutoff de %g, il y a %d clusters", threshold, length(unique(clusters))))
  
  # Create the cluster table
  cluster_table <- hhtransmission_df %>%
    select(pid, hhid, phy_hhtrans_mem, hh_infected) %>%
    rename(cutoff_hhtransmission = phy_hhtrans_mem,
           cutoff_hhtransmission_houselevel = hh_infected)
  
  infected_hh_table <- hhtransmission_df %>%
    group_by(hhid) %>%
    summarise(hh_infected = any(hh_infected)) %>%
    rename(!!title := hh_infected)
  
  
  return(list(
    plot = p,
    cluster_table = cluster_table,
    infected_hh_table = infected_hh_table,
    summary = sprintf("Pour un cutoff de %g, il y a %d clusters", threshold, length(unique(clusters))),
    cluster_object = hc,
    patristic_distances = patristic_distances
  ))
}


#' TreeCluster-based Analysis of Household Transmission
#'
#' @description Uses external cluster assignments from TreeCluster to identify
#' potential household transmission events. Visualizes and summarizes results.
#'
#' @param tree_file Path to the phylogenetic tree file.
#' @param treecluster_file Path to TreeCluster output file (tabular).
#' @param threshold Distance cutoff (used in output text only).
#' @param title Title for plots and summaries.
#'
#' @return A list with the ggtree plot, cluster table, household summary,
#' text summary, clustering object, and distance matrix.

treecluster_analysis <- function(tree_file, treecluster_file, threshold, title, type = "Newick") {
  verdi_df <- read.csv("../datasets/verdi_seq_df.csv")
  phylo_df <- verdi_df %>% filter(used_in_phylogeny == "YES")
  divergence_tree <- read_and_prune_tree(tree_file, phylo_df$sampleid, type)
  
  treecluster <- read.table(treecluster_file, header = TRUE) %>%
    rename(sampleid = SequenceName, cluster = ClusterNumber) %>%
    filter(sampleid %in% phylo_df$sampleid)
  
  singleton_count <- sum(treecluster$cluster == -1)
  unique_cluster_count <- length(unique(treecluster$cluster[treecluster$cluster != -1]))
  
  df_clusters <- treecluster
  df_clusters$cluster[df_clusters$cluster == -1] <- NA
  hhtransmission_df <- prepare_annotation_data(verdi_df, df_clusters) %>%
    mutate(hh_infected = ifelse(is.na(cluster), FALSE, hh_infected),
           phy_hhtrans_mem = ifelse(index == 0 & hh_infected, TRUE, FALSE))
  annotation_data <- hhtransmission_df %>%
    filter(hh_infected == TRUE) %>%
    select(sampleid, hhid, index)
  
  p <- make_ggtree_plot(divergence_tree, annotation_data,
                        title = paste(title," WORKFLOW (122 samples, cutoff", threshold, ",",
                                      unique_cluster_count, "clusters)"))
  print(p)
  message(sprintf("For a cutoff of %g: %d clusters and %d singletons",
                  threshold, unique_cluster_count, singleton_count))
  
  # CLUSTER TABLE RES
  
  # Create the cluster table
  cluster_table <- hhtransmission_df %>%
    select(pid, phy_hhtrans_mem, hh_infected) %>%
    rename(hhtrans_mem = phy_hhtrans_mem,
           hhtrans_house = hh_infected)
  
  infected_hh_table <- hhtransmission_df %>%
    group_by(hhid) %>%
    summarise(hh_infected = any(hh_infected)) %>%
    rename(!!title := hh_infected)
  
  
  
  return(list(
    plot = p,
    cluster_table = cluster_table,
    infected_hh_table = infected_hh_table,
    summary = sprintf("For a cutoff of %g: %d clusters and %d singletons",
                      threshold, unique_cluster_count, singleton_count),
    cluster_object = treecluster))
}

