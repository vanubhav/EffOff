tryCatch(
    {
        library(reticulate)

        # 1. Define a persistent location for the environment (e.g., inside the project)
        # Using a local folder ensures it doesn't get "lost" or mixed up with other system envs
        env_name <- "./effoff_env"

        message("--> Creating persistent virtual environment at: ", env_name)
        virtualenv_create(envname = env_name)

        message("--> Installing earthengine-api and dependencies...")
        virtualenv_install(envname = env_name, packages = c("earthengine-api", "numpy"))

        message("--> Environment setup complete!")

        # Get the absolute path to the python executable
        python_path <- virtualenv_python(env_name)
        abs_python_path <- normalizePath(python_path, winslash = "/")

        message("\nIMPORTANT: To make this setting PERMANENT, please follow these steps:")
        message("1. Run: usethis::edit_r_environ()")
        message("2. Add the following line to the file (if .Renviron opens):")
        message(sprintf("RETICULATE_PYTHON=\"%s\"", abs_python_path))
        message("3. Restart your R session.")
    },
    error = function(e) {
        message("--> FAILURE: ", e$message)
    }
)
