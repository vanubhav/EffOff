tryCatch(
    {
        message("--> Checking reticulate configuration...")
        config <- reticulate::py_config()
        print(config)

        message("--> Installing earthengine-api...")
        # Try installing into the currently active environment found by reticulate
        reticulate::py_install("earthengine-api", pip = TRUE)

        message("--> Verifying installation...")
        ee <- reticulate::import("ee")
        print(ee$String("Success!")$getInfo())

        message("--> SUCCESS: earthengine-api installed and working.")
    },
    error = function(e) {
        message("--> FAILURE: ", e$message)
    }
)
