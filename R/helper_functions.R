run_protein_idr_summary_on_fasta <- function(fasta_file) {
  # Load required packages
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    stop("Please install the 'seqinr' package.")
  }
  if (!requireNamespace("idpr", quietly = TRUE)) {
    stop("Please install the 'idpr' package from GitHub: devtools::install_github('cjf00000/idpr')")
  }

  library(seqinr)
  library(idpr)

  # Read sequences from FASTA file
  sequences <- read.fasta(file = fasta_file, seqtype = "AA", as.string = TRUE, set.attributes = FALSE)

  # Loop over sequences and process each one
  for (i in seq_along(sequences)) {
     seq_name <- names(sequences)[i]
     seq_string <- sequences[[i]]

     cat("\n=====================================\n")
     cat("Processing sequence:", seq_name, "\n")
     cat("=====================================\n\n")

     # Check valid input
     if (!is.character(seq_string) || length(seq_string) != 1) {
       warning(paste("Skipping", seq_name, "- invalid sequence format"))
       next
     }

     # Run summaries
     cat("=== IDProfile Summary ===\n")
     profile <- idprofile(seq_string)
     print(profile)

     cat("\n=== IUPred Anchor Regions ===\n")
     anchor <- iupredAnchor(seq_string)
     print(anchor)

     cat("\n=== IUPred Redox Regions ===\n")
     redox <- iupredRedox(seq_string)
     print(redox)
   }
  return(sequences)
}
