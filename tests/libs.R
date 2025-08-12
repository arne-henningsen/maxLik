cat("Installed to:\n")
cat(".Library:", .Library, "\n")
cat(".Library.site:", .Library.site, "\n")
cat(".libPaths:\n")
print(.libPaths())

options(width=9999)
installed.packages()

cat("Installed packages:\n")
list.files(.libPaths()[1])

