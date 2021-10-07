library(tidyverse)


### Conversion ###


# Get the rows of a df as a list
row_list <- function(df) unname(split(df, seq(nrow(df))))


# Convert decimal to binary, optionally to a given number of bits
intToBin <- function(x, n_bits = NULL) {
  long <- rev(as.integer(intToBits(x)))
  tail(long, if_else(is.null(n_bits), length(long) - match(1, long, length(long)) + 1, n_bits))
}


numToInt <- function(x) { as.integer(unlist(round(x))) }
numToFac <- function(x, levels) { factor(levels[numToInt(x)], levels) }
# as.numeric is guaranteed to return the level number even if the level is an integer itself
facToInt <- function(x) { as.numeric(x) } # Old: setNames(c(1:length(levels(x))), nm = levels(x))[x] # Older: as.numeric(levels(x))[x]
# Simple check:
#all(data == toClassFac(toClassInt(data)))



### Iterable-related ###


# Sort by function, with optional extra arguments to f
sort_by <- function(xs, f, ..., decreasing = F) xs[order(sapply(xs, f, ...), decreasing = decreasing)]
# Sort by i-th element
sort_by_elem <- function(xss, i, decreasing = F) xss[order(sapply(xss, '[[', i), decreasing = decreasing)]
# Sort by arbitrary pluck
sort_by_pluck <- function(xss, ..., decreasing = F) xss[order(sapply(xss, pluck, ...), decreasing = decreasing)]

# Filter by uniqueness of function application output; might want to %>% unname()
unique_by <- function(xs, f, ...) split(xs, sapply(xs, f, ...)) %>% sapply(head, 1)
# Filter by uniqueness of i-th element; might want to %>% unname()
unique_by_elem <- function(xs, i) split(xs, sapply(xs, '[[', i)) %>% sapply(head, 1)
# Filter by uniqueness of arbitraty pluck; might want to %>% unname()
unique_by_pluck <- function(xs, ...) split(xs, sapply(xs, pluck, ...)) %>% sapply(head, 1)


# Powerset, optionally up to a given maximum subset cardinality
# NOTES:
#   - Assumes all elements are different, so input containing identical elements will result in duplicates; add a %>% unique() if problematic
#   - Does not return the empty set; add manually if required
powerset <- function(xs, up_to = NULL) map(1:ifelse(is.null(up_to), length(xs), up_to), ~ combn(xs, ., simplify = F)) %>% unlist(recursive = F)


# Zip any number of items/vectors/lists to the length of the shortest one, i.e. standard zip
lower_zip <- function(...) {
  xss <- list(...)
  min_len <- min(map(xss, length) %>% unlist())
  nice_lists <- map(xss, function(xs) as.vector(xs)[1:min_len])
  pmap(nice_lists, list)
}
# Zip any number of items/vectors/lists to the length of the longest one by repeating the shorter ones
upper_zip <- function(...) {
  xss <- list(...)
  max_len <- max(map(xss, length) %>% unlist())
  nice_lists <- map(xss, function(xs) rep(as.vector(xs), length.out = max_len))
  pmap(nice_lists, list)
}
# Zip any number of items/vectors/lists to the length of the longest one by repeating the last elements of the shorter ones
upper_zip_last <- function(...) {
  xss <- list(...)
  max_len <- max(map(xss, length) %>% unlist())
  nice_lists <- map(xss, function(xs) c(as.vector(xs), rep(xs[length(xs)], length.out = max_len - length(xs))))
  pmap(nice_lists, list)
}



### Other ###


# Partition a data frame into either a given number of random subsets
# or into random subsets of given proportional sizes and names
# Set the argument 'parts' to an integer for the former and to a named list of proportions adding up to 1 for the latter
splitBy <- function(df, parts) {
    labels <- if (length(parts) == 1) cut(seq(nrow(df)), parts, labels = 1:parts)
              else cut(seq(nrow(df)), nrow(df) * cumsum(c(0, parts)), labels = names(parts))
    split(df, sample(labels))
}


# Compute the overlap of two intervals
  # IMPORTANT NOTE: This is not a vectorised function, therefore DO
  #   EITHER: rowwise %>% mutate(... intervalOverlap(...) ...) %>% ungroup
  #   OR: mutate(... Vectorize(intervalOverlap)(...) ...)
  # If this is not done it will return the result of the first evaluation and that will be it for all rows
intervalOverlap <- function(a, b) max(0, min(a[2], b[2]) - max(a[1], b[1]))


