
showMaxControl <- function(object) {
   cat("A 'MaxControl' object with slots:\n")
   for(s in slotNames(object)) {
      cat(s, "=", slot(object, s), "\n")
   }
}

setMethod("show", "MaxControl", showMaxControl)
