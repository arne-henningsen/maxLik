
showMaxControl <- function(object) {
   for(s in slotNames(object)) {
      cat(s, "=", slot(object, s), "\n")
   }
}

setMethod("show", "MaxControl", showMaxControl)
