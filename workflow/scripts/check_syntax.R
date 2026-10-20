# Parse only; do not execute analyses or read scientific inputs.
here::i_am("workflow/scripts/check_syntax.R")
scripts <- list.files(here::here("workflow", "scripts"), pattern = "[.]R$", full.names = TRUE)
for (script in scripts) {
    invisible(parse(script))
    message("OK: ", basename(script))
}
notebook <- here::here("analysis", "22-scrnaseq-data-analysis", "scripts",
                      "03-score-patients-on-the-basis-of-erap2-module.rmd")
extracted <- tempfile(fileext = ".R")
knitr::purl(notebook, output = extracted, documentation = 0, quiet = TRUE)
invisible(parse(extracted))
unlink(extracted)
message("OK: original S7 notebook R chunks")
