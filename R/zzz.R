# .onLoad <- function(libname, pkgname)
# {
#   library.dynam("admix", pkgname, libname)
# }

admixStartupMessage <- function()
{
  msg <- c(paste0(
    "This is package admix, version ",
    utils::packageVersion("admix")),
    "\n-------------------------------\n",
    "Type 'citation(\"admix\")' for citing this R package in publications.",
    "\n-------------------------------\n",
    "This work was partly conducted within the Research Chair DIALog under the aegis of the Risk Foundation, an initiative by CNP Assurances.\n")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # unlock .admix variable allowing its modification
  #unlockBinding(".admix", asNamespace("admix"))
  # startup message
  msg <- admixStartupMessage()
  if(!interactive())
    msg[1] <- paste("Package 'admix' version", packageVersion("admix"))
  base::packageStartupMessage(msg)
  base::invisible()
}
