maximMessage <- function(code) {
   message <- switch(code,
                     "1" = "gradient close to zero",
                     "2" = "successive function values within tolerance limit",
                     "3" = paste("Last step could not find a value above the",
                     "current.\nBoundary of parameter space?",
                     " \nConsider switching to a more robust optimisation method temporarily."),
                     "4" = "Iteration limit exceeded.",
                     "5" = "Infinite value",                        
                     "6" = "Infinite gradient",                        
                     "7" = "Infinite Hessian",                        
                     "100" = "Initial value out of range.",
                     paste("Code", code))
   return(message)
}
