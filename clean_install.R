tryCatch(
    {
        message("--> removing package 'EffOff'...")
        remove.packages("EffOff")
        message("--> Package removed.")
    },
    error = function(e) {
        message("--> Warning during removal: ", e$message)
    }
)

tryCatch(
    {
        message("--> Installing package via devtools...")
        devtools::install("a:/ColumbiaUniversity/AntiGravity_Workspaces/CAMPA/EffOff", upgrade = "never", quick = TRUE)
        message("--> Installation complete.")
    },
    error = function(e) {
        message("--> FAILURE: Installation failed.")
        message(e)
        quit(status = 1)
    }
)
