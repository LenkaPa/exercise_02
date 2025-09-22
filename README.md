# Introduction to (R and) R/Bioconductor and Regular Expressions

## Introduction to (R and) R/Bioconductor
### Task 1
* Load the DNA sequence `fishes.fna.gz` using functions from the `seqinr` package and the `Biostrings` package.
Note the differences between the created variables.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::available()
BiocManager::install("Biostrings")
library(package)

if (!requireNamespace("seqinr", quietly=TRUE)) {
  install.packages("seqinr")
}

library(Biostrings)
library(seqinr)

seq_biostrings <- Biostrings::readDNAStringSet("fishes.fna.gz")
seq_seqinr <- seqinr::read.fasta(file = gzfile("fishes.fna.gz"))
# Differences summary:
# - seqinr: list of character vectors (one character per base), names accessible via names()
# - Biostrings: DNAStringSet object, supports fast vectorized ops (width(), subseq(), writeXStringSet(), etc.)

### Task 2
* Next, focus on the `Biostrings` package. Practice working with loaded data:
seq <- seq_biostrings
    * Check the number of loaded sequences:
        ```R
        length(seq)
        ```
    * Determine the lengths of each sequence:
        ```R
        width(seq[1])
        ```
    * View the sequence names (FASTA headers):
        ```R
        names(seq)
        ```
    * Assign the first sequence including the name to the variable `seq1`:
        ```R
        seq1 <- seq[1]
        ```
    * Assign the first sequence without the name to the variable `seq1_sequence`:
        ```R
        seq1_sequence <- seq[[1]]
        ```
    * Assign the first sequence as a vector of characters to the variable `seq1_string`:
        ```R
        seq1_string <- toString(seq[1])
        ```
    * Learn more about the `XStringSet` class and the `Biostrings` package:
        ```R
        help(XStringSet)
        ```

### Task 3
 * Globally align the two selected sequences using the BLOSUM62 matrix, a gap opening cost of -1 and a gap extension cost of 1.
library(Biostrings)

# load BLOSUM62 substitution matrix
data(BLOSUM62)

# pick two sequences from your loaded DNAStringSet or AAStringSet
seq1 <- seq_biostrings[[1]]
seq2 <- seq_biostrings[[2]]



# If your file is PROTEINS (AAStringSet), you can align directly:
BiocManager::install("pwalign")
alignment <- pwalign::pairwiseAlignment(
  pattern = seq1,
  subject = seq2,
  substitutionMatrix = BLOSUM62,
  gapOpening = -1,
  gapExtension = 1,
  type = "global"
)

# view results
alignment
score(alignment)

## Regular Expressions
### Task 4
* Practice working with regular expressions:
    * Create a list of names, e.g.:
        ```R
        names_list <- c("anna", "jana", "kamil", "norbert", "pavel", "petr", "stanislav", "zuzana")
        ```
    * Search for name `jana`:
        ```R
        grep("jana", names_list, perl = TRUE)
        ```
    * Search for all names containing letter `n` at least once:
        ```R
        grep("n+", names_list, perl = TRUE)
        ```
    * Search for all names containing letters `nn`:
        ```R
        grep("n{2}", names_list, perl = TRUE)
        ```
    * Search for all names starting with `n`:
        ```R
        grep("^n", names_list, perl = TRUE)
        ```
    * Search for names `Anna` or `Jana`:
        ```R
        grep("anna|jana", names_list, perl = TRUE)
        ```
    * Search for names starting with `z` and ending with `a`:
        ```R
        grep("^z.*a$", names_list, perl = TRUE)
        ```

### Task 5
* Load an amplicon sequencing run from 454 Junior machine `fishes.fna.gz`.
* Get a sequence of a sample (avoid conditional statements), that is tagged by forward and reverse MID `ACGAGTGCGT`.
* How many sequences are there in the sample?

# 1. Load sequences and convert to character for simple pattern search
seqs <- readDNAStringSet("fishes.fna.gz")
seq_chars <- as.character(seqs)

# 2. Define MID and its reverse complement
# mid <- DNAString("ACGAGTGCGT")
#mid_rc <- reverseComplement(mid)
mid <- "ACGAGTGCGT"
mid_rc <- as.character(reverseComplement(DNAString(mid)))

# 3. Find seq with MID
# nějaký MID nebo jeho reverzní komplement, pak libovolná sekvence, pak znovu MID nebo jeho RC
pattern <- paste0("(", mid, "|", mid_rc, ").*(", mid, "|", mid_rc, ")")

# 4. Extract these sequences with grepl() -> vrací TRUE/FALSE
sample_idx <- grepl(pattern, seq_chars, perl = TRUE)

# 5. Report how many sequences -> počet TRUE
sum(sample_idx)


### Task 6
* Create a function `Demultiplexer()` for demultiplexing of sequencing data.

* Input:
    * a string with path to fasta file
    * a list of forward MIDs
    * a list of reverse MIDs
    * a list of samples labels

* Output:
    * fasta files that are named after the samples and contain sequences of the sample without MIDs (perform MID trimming)
    * table named `report.txt` containing samples‘ names and  the number of sequences each sample has

* Check the functionality again on the `fishes.fna.gz` file, the list of samples and MIDs can be found in the corresponding table `fishes_MIDs.csv`.

Demultiplexer <- function(fasta_path, forward_mids, reverse_mids, sample_labels,
                          out_dir = "demultiplexed", report_name = "report.txt") {
  # check inputs
  if (length(forward_mids) != length(reverse_mids) ||
      length(forward_mids) != length(sample_labels)) {
    stop("forward_mids, reverse_mids and sample_labels must have the same length.")
  }
  
  # create output directory if not exists
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # load sequences
  seqs <- Biostrings::readDNAStringSet(fasta_path)
  seq_chars <- as.character(seqs)
  seq_names <- names(seqs)
  
  # prepare report and output list
  report <- data.frame(sample = sample_labels,
                       count = integer(length(sample_labels)),
                       stringsAsFactors = FALSE)
  trimmed_seqs <- vector("list", length(sample_labels))
  names(trimmed_seqs) <- sample_labels
  
  for (i in seq_along(sample_labels)) {
    f <- forward_mids[i]
    r <- reverse_mids[i]
    
    # reverse complement of reverse MID
    r_rc <- as.character(reverseComplement(DNAString(r)))
    f_rc <- as.character(reverseComplement(DNAString(f)))
    
    # match sequences that start with forward MID and end with reverse MID
    pattern1 <- paste0("^", f, ".*", r_rc, "$")
    idx1 <- grepl(pattern1, seq_chars, perl = TRUE)
    
    pattern2 <- paste0("^", r, ".*", f_rc, "$")
    idx2 <- grepl(pattern2, seq_chars, perl = TRUE)
    
    idx <- idx1 | idx2
    if (sum(idx) == 0) next
    
    # trim MIDs
    trimmed <- seq_chars[idx]
    trimmed <- sub(f, "", trimmed)
    trimmed <- sub(r_rc, "", trimmed)
    trimmed <- sub(r, "", trimmed)
    trimmed <- sub(f_rc, "", trimmed)
    
    # create DNAStringSet
    trimmed_set <- Biostrings::DNAStringSet(trimmed)
    names(trimmed_set) <- seq_names[idx]
    
    # save to FASTA
    out_fasta <- file.path(out_dir, paste0(sample_labels[i], ".fasta"))
    Biostrings::writeXStringSet(trimmed_set, filepath = out_fasta)
    
    trimmed_seqs[[i]] <- trimmed_set
    report$count[i] <- length(trimmed_set)
    
    message("Sample ", sample_labels[i], ": ", length(trimmed_set), 
            " sequences saved to ", out_fasta)
  }
  
  # save report to file
  report_file <- file.path(out_dir, report_name)
  write.table(report, file = report_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Report saved to ", report_file)
  
  return(list(report = report, trimmed_seqs = trimmed_seqs))
}



# načtení potřebné knihovny
library(Biostrings)

# načtení CSV s MIDs a vzorky
mids <- read.csv("fishes_MIDs.csv", sep = ";", stringsAsFactors = FALSE)
colnames(mids)
str(mids)

# zkontrolovat, zda existují potřebné sloupce
required_cols <- c("SampleID", "FBarcodeSequence", "RBarcodeSequence")
if (!all(required_cols %in% colnames(mids))) {
  stop("CSV musí obsahovat sloupce: ", paste(required_cols, collapse = ", "))
}

forward_mids <- as.character(mids$FBarcodeSequence)
reverse_mids <- as.character(mids$RBarcodeSequence)
sample_labels <- as.character(mids$SampleID)

# spustit demultiplexování
result <- Demultiplexer(
  fasta_path = "fishes.fna.gz",
  forward_mids = forward_mids,
  reverse_mids = reverse_mids,
  sample_labels = sample_labels,
  out_dir = "demultiplexed_fishes",
  report_name = "report.txt"
)

# přístup k reportu
report <- result$report
print(report)

# přístup k oříznutým sekvencím pro jednotlivé vzorky
trimmed_sequences <- result$trimmed_seqs

# např. počet sekvencí pro první vzorek
cat("Počet sekvencí pro první vzorek (", names(trimmed_sequences)[1], "): ",
    length(trimmed_sequences[[1]]), "\n")



## Download files from GitHub
<details>
<summary>Basic Git settings</summary>

>* Configure the Git editor
>    ```bash
>    git config --global core.editor notepad
>    ```
>* Configure your name and email address
>    ```bash
>    git config --global user.name "Zuzana Nova"
>    git config --global user.email z.nova@vut.cz
>    ```
>* Check current settings
>    ```bash
>    git config --global --list
>    ```
>
</details>

* Create a fork on your GitHub account.
  On the GitHub page of this repository find a <kbd>Fork</kbd> button in the upper right corner.

* Clone forked repository from your GitHub page to your computer:
```bash
git clone <fork repository address>
```
* In a local repository, set new remote for a project repository:
```bash
git remote add upstream https://github.com/mpa-prg/exercise_02.git
```

#### Send files to GitHub
Create a new commit and send new changes to your remote repository.
* Add file to a new commit.
```bash
git add <file_name>
```
* Create a new commit, enter commit message, save the file and close it.
```bash
git commit
```
* Send a new commit to your GitHub repository.
```bash
git push origin main
```
