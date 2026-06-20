# Load the six paper application datasets from the minimal data bundle.

paper_application_data_paths <- function() {
  list(
    moreno_sheep = file.path("data", "moreno_sheep", "edges.csv"),
    strauss_2019b = file.path("data", "Strauss_2019b", "edges.csv"),
    mountain_goats = file.path("data", "mountain_goats", "adjacency_matrix.csv"),
    citations_data = file.path("data", "citations_data", "adjacency_matrix.csv"),
    macaques_data = file.path("data", "macaques_data", "edge_list.tsv"),
    high_school = file.path("data", "high_school", "edges.csv")
  )
}

read_edge_list_adjacency <- function(path, weighted = TRUE, aggregate_counts = FALSE) {
  if (!file.exists(path)) {
    stop("Missing edge-list file: ", path, call. = FALSE)
  }

  edge_list <- if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) {
    utils::read.table(path, header = FALSE, comment.char = "#")
  } else {
    utils::read.csv(path, header = FALSE, comment.char = "#", strip.white = TRUE)
  }
  edge_list <- as.data.frame(edge_list, stringsAsFactors = FALSE)

  if (aggregate_counts) {
    names(edge_list)[1:2] <- c("source", "target")
    edge_list$weight <- 1L
    edge_list <- stats::aggregate(weight ~ source + target, data = edge_list, FUN = sum)
  } else {
    names(edge_list)[1:3] <- c("source", "target", "weight")
    if (!weighted) edge_list$weight <- 1L
  }

  edge_list$source <- as.integer(edge_list$source)
  edge_list$target <- as.integer(edge_list$target)
  edge_list$weight <- as.integer(edge_list$weight)

  if (min(edge_list$source, edge_list$target, na.rm = TRUE) == 0L) {
    edge_list$source <- edge_list$source + 1L
    edge_list$target <- edge_list$target + 1L
  }

  n_nodes <- max(c(edge_list$source, edge_list$target), na.rm = TRUE)
  node_ids <- as.character(seq_len(n_nodes))
  adjacency <- matrix(0L, n_nodes, n_nodes, dimnames = list(node_ids, node_ids))

  for (row_id in seq_len(nrow(edge_list))) {
    i <- edge_list$source[[row_id]]
    j <- edge_list$target[[row_id]]
    w <- edge_list$weight[[row_id]]
    if (w > 0L) adjacency[i, j] <- adjacency[i, j] + w
  }

  diag(adjacency) <- 0L
  adjacency
}

load_application_adjacency <- function(dataset) {
  dataset_paths <- paper_application_data_paths()
  dataset <- match.arg(dataset, choices = names(dataset_paths))

  if (dataset == "mountain_goats" || dataset == "citations_data") {
    adjacency <- as.matrix(
      utils::read.csv(dataset_paths[[dataset]], row.names = 1, check.names = FALSE)
    )
    diag(adjacency) <- 0L
  } else if (dataset == "macaques_data") {
    adjacency <- read_edge_list_adjacency(dataset_paths[[dataset]], weighted = TRUE)
  } else if (dataset == "high_school" || dataset == "moreno_sheep") {
    adjacency <- read_edge_list_adjacency(dataset_paths[[dataset]], weighted = TRUE)
  } else if (dataset == "strauss_2019b") {
    adjacency <- read_edge_list_adjacency(dataset_paths[[dataset]], weighted = FALSE, aggregate_counts = TRUE)
  } else {
    stop("Unknown application dataset: ", dataset, call. = FALSE)
  }

  stopifnot(nrow(adjacency) == ncol(adjacency))
  colnames(adjacency) <- rownames(adjacency)
  adjacency
}
