
.First.lib<-function(libname,pkgname){
    if (TRUE)
        library.dynam("netCDF",pkgname,libname)
    else
        stop("You don't have netCDF installed.")
}
 
