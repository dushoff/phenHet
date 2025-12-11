library(shellpipes)
library(Rcpp)

sourceCpp(matchFile(exts=c("cpp", "Cpp")))

dlls <- getLoadedDLLs()
sc <- grep("sourceCpp", names(dlls))
stopifnot(length(sc)==1)

dll_path <- dlls[[sc]][["path"]]

print(dll_path)

stopifnot(file.copy(dll_path, targetname(ext=".so"), overwrite = TRUE))

saveEnvironment()
