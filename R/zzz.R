damdaStartupMessage <- function()
{
  msg <- c(paste0(
    "damda --- version ", utils::packageVersion("damda")),
    "\nType 'citation(\"damda\")' for citing this R package in publications."
    )
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- damdaStartupMessage()
  packageStartupMessage(msg)      
  invisible()
}