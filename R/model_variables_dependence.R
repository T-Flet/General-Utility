# These functions are used when scheduling model selection steps if variable dependence is taken into account 


# Generate the list of variables with their dependencies, i.e. the other variables which they need to be present to exist
gen_dependent_vars <- function(numeric_vars, factor_vars, max_power = 2, interactions = F) {
  terms <- if (interactions) powerset(rep(numeric_vars, each = max_power) %>% append(factor_vars), max_power) %>% unique()
           else map(numeric_vars, function(v) map(1:max_power, ~ rep(v, .x))) %>% reduce(append)
  terms %>%
    map(function(var) tibble(var) %>% group_by(var)) %>% # Using 'function' instead of ~ already names the column
    map(~ list(id = .x %>% summarise(pow = n()), dp = .x %>% mutate(pow = 1:n()) %>% ungroup()))
}


# Transform a term list (contining tibbles) to its string version
#   NOTE: May want to transform it to a named list of lists of dependencies:
#           setNames(map(RESULT, ~ .x$dp), map(RESULT, ~ .x$id))
terms_to_str <- function(term_deps) {
  df_str <- function(term) (term %>% mutate(string = if_else(pow > 1, paste0('I(', var, '^', pow, ')'), var)))$string %>% paste0(collapse = ':')
  proper <- function(pss) map(pss, ~ bind_rows(.x) %>% group_by(var) %>% summarise(pow = max(pow)) %>% arrange()) %>% unique() %>% head(-1)
  
  term_deps %>% map(function(term) list(
    id = df_str(term$id) %>% paste0(collapse = ':'),
    dp = powerset(row_list(term$dp)) %>% proper() %>% map(df_str) %>% unlist()))
}


# Candidates for removal given the already removed terms
removable_terms <- function(term_deps, removed) {
  remaining <- term_deps %>% discard(~ .x$id %in% removed)
  remaining %>% discard(function(term) any(map(remaining, ~ term$id %in% .$dp) %>% unlist())) %>% map(~ .$id) %>% unlist()
}


# Identify a given model's terms' degrees and base terms; useful to compute further information
#   NOTES:
#     - The result is a list with an unnamed element for each model term;
#         this includes repetitions for factor levels;
#         might want to %>% unique_by_elem('sorted'), which also names the list by the field
#     - The usefulness of the 'sorted' field lies in lining up with gen_dependent_vars outputs
get_terms_details <- function(mod) {
  factors <- keep(colnames(mod$data), ~ is.factor(mod$data[[.x]])) # These may however not be in the formula
  
  tidy(mod)$term %>% discard(~ .x == '(Intercept)') %>% map(~
      list(original = .x, details = map(str_split(.x, ':') %>% unlist(), function(part) {
          inner <- str_match(part, '^I\\((.+)\\^([0-9]+)\\)$') # Full match is [1], base term is [2] and power is [3]
          base_factor <- detect(map(factors, ~ list(f = .x, match = str_starts(part, fixed(.x)))), ~ .x$match)$f
          list(real = ifelse(is.null(base_factor), part, base_factor), # If factor then just its name
               degree = ifelse(is.na(inner[3]), 1, strtoi(inner[3])),
               base = ifelse(is.null(base_factor), ifelse(is.na(inner[2]), part, inner[2]), base_factor))
        })
      ) %>% (function(term) {
        identified <- str_c(map(term$details, ~ .x$real), collapse = ':')
        sorted <- sort_by_elem(term$details, 'base') %>% map(~ .x$real) %>% str_c(collapse = ':')
        list(original = term$original, # the name as it appears in the summary, containing individual factor levels
             identified = identified, # 'original' with factor levels identified, i.e. the term as in a formula
             sorted = sorted, # 'identified' with lexicographically reordered interaction components
             degree = sum(map(term$details, ~ .x$degree) %>% unlist()),
             bases = map(term$details, ~ .x$base) %>% unlist() ) # The base terms involved (>1 for interactions)
      })
    )
}


# Generate the variable dependency list (in string form) from a given model
#   NOTE: the output's fields have the term_details' 'identified' field in place of the 'sorted' one
gen_dependent_vars_from_model <- function(mod, term_details = NULL) {
  term_details <- (if (is.null(term_details)) get_terms_details(mod) else term_details) %>% unique_by_elem('sorted')
  
  present_terms <- map(term_details, ~ .x$identified) # Names are sorted, values are identified
  base_terms <- map(term_details, ~ .x$bases) %>% unlist(recursive = F) %>% unique()
  interaction_terms <- keep(term_details, ~ str_detect(.x$original, ':')) %>% map(~ .x$sorted)
  factors <- intersect(base_terms, keep(colnames(mod$data), ~ is.factor(mod$data[[.x]]))) # I.e. the factors which are present
  
  max_power <- max(map(term_details, ~ .x$degree) %>% unlist())
  interactions <- length(interaction_terms) > 0
  
  gen_dependent_vars(setdiff(base_terms, factors), factors, max_power = max_power, interactions = interactions) %>%
    terms_to_str() %>% keep(~ .x$id %in% names(present_terms)) %>% map(function(tdp) list( # Reorder interaction components
      id = ifelse(tdp$id %in% interaction_terms, present_terms[[tdp$id]], tdp$id),
      dp = map(tdp$dp, ~ ifelse(.x %in% interaction_terms, present_terms[[.x]], .x)) %>% unlist()
    ))
}


# ## Example generation and steps if 'b' is not significant
# AA <- gen_dependent_vars(c('a', 'b'), c('f'), max_power = 3, interactions = T)
# AA[[length(AA)]]
# BB <- terms_to_str(AA)
# BB[[length(BB)]]
# removable_terms(BB, c())
# removable_terms(BB, c('I(b^2):f'))
# removable_terms(BB, c('I(b^2):f', 'a:I(b^2)'))
# removable_terms(BB, c('I(b^2):f', 'a:I(b^2)', 'I(b^3)'))
# removable_terms(BB, c('I(b^2):f', 'a:I(b^2)', 'I(b^3)', 'I(b^2)'))
# removable_terms(BB, c('I(b^2):f', 'a:I(b^2)', 'I(b^3)', 'I(b^2)', 'I(a^2):b'))
# removable_terms(BB, c('I(b^2):f', 'a:I(b^2)', 'I(b^3)', 'I(b^2)', 'I(a^2):b', 'a:b:f'))
# 
# ## Example of generation from a model
# mod <- glm(Petal.Length ~ Species + Sepal.Length + Sepal.Length:Species + I(Sepal.Length^2) + Sepal.Width + Petal.Width + Sepal.Width:I(Petal.Width^2), iris, family = Gamma())
# term_details <- get_terms_details(mod)
# term_details
# gen_dependent_vars_from_model(mod, term_details) %>% map(~ .x$id) %>% unlist()


