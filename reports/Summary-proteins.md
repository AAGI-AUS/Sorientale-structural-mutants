# Summary of IAA proteins
Dr. Kristina Gagalova

- [Load the protein dataset](#load-the-protein-dataset)
- [Calculate summary of proteins](#calculate-summary-of-proteins)

## Load the protein dataset

``` r
library(seqinr)
library(idpr)

# Read sequences from FASTA file
sequences <- read.fasta(file = "../data/sequences/S.orientale-muts.fasta", seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
```

    Warning in readLines(file): incomplete final line found on
    '../data/sequences/S.orientale-muts.fasta'

## Calculate summary of proteins

``` r
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
      profile <- idprofile(seq_string, proteinName=seq_name)
      print(profile[1:5])
      
    }
```


    =====================================
    Processing sequence: SoIAA2-WT 
    =====================================

    === IDProfile Summary ===
    [[1]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-1.png)


    [[2]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-2.png)


    [[3]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-3.png)


    [[4]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-4.png)


    [[5]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-5.png)



    =====================================
    Processing sequence: SoIAA2-d9 
    =====================================

    === IDProfile Summary ===
    [[1]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-6.png)


    [[2]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-7.png)


    [[3]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-8.png)


    [[4]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-9.png)


    [[5]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-10.png)



    =====================================
    Processing sequence: SoIAA2-d27 
    =====================================

    === IDProfile Summary ===
    [[1]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-11.png)


    [[2]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-12.png)


    [[3]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-13.png)


    [[4]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-14.png)


    [[5]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-15.png)



    =====================================
    Processing sequence: SoIAA2-d33 
    =====================================

    === IDProfile Summary ===
    [[1]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-16.png)


    [[2]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-17.png)


    [[3]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-18.png)


    [[4]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-19.png)


    [[5]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-20.png)



    =====================================
    Processing sequence: SoIAA2-dnd 
    =====================================

    === IDProfile Summary ===
    [[1]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-21.png)


    [[2]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-22.png)


    [[3]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-23.png)


    [[4]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-24.png)


    [[5]]

![](Summary-proteins_files/figure-commonmark/idr-proteins-25.png)
