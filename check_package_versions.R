# check_package_versions.R
#
# Stops unless R and every package of renv.lock are installed in exactly the recorded version
# (for GitHub packages: the recorded commit). Sourced at the top of every R script and notebook
# of the pipeline, because renv only records the versions and does not compare them at run time.
#
# To return to the recorded packages: renv::restore(). R itself is installed system-wide, not
# by renv, and has to be changed there.

local({
  lockfile <- renv::lockfile_read("renv.lock")
  problems <- character(0)

  running_r_version <- as.character(getRversion())
  if (running_r_version != lockfile$R$Version) {
    problems <- c(problems, sprintf("R: running %s, recorded %s", running_r_version, lockfile$R$Version))
  }

  for (recorded in lockfile$Packages) {
    description_file <- system.file("DESCRIPTION", package = recorded$Package)
    if (description_file == "") {
      problems <- c(problems, sprintf("%s: not installed, recorded %s", recorded$Package, recorded$Version))
      next
    }
    installed <- read.dcf(description_file, fields = c("Version", "RemoteSha"))[1, ]

    if (installed[["Version"]] != recorded$Version) {
      problems <- c(problems, sprintf("%s: installed %s, recorded %s",
                                      recorded$Package, installed[["Version"]], recorded$Version))
    } else if (identical(recorded$Source, "GitHub") && !identical(installed[["RemoteSha"]], recorded$RemoteSha)) {
      problems <- c(problems, sprintf("%s: installed commit %s, recorded commit %s",
                                      recorded$Package, installed[["RemoteSha"]], recorded$RemoteSha))
    }
  }

  if (length(problems) > 0) {
    stop("The installed versions differ from renv.lock:\n  ", paste(problems, collapse = "\n  "), call. = FALSE)
  }
})
