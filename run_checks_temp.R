tryCatch(
    {
        library(devtools)

        cat("Updating documentation...\n")
        devtools::document(".")

        # Run a quick check (checking everything might be slow, but let's try a standard check)
        cat("Running check...\n")
        # Using check(cran = FALSE) to speed it up and avoid strict CRAN policies if not needed immediately
        devtools::check(".", cran = FALSE)
    },
    error = function(e) {
        cat("Error occurred:\n")
        print(e)
        quit(status = 1)
    }
)
