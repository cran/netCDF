##
## NetCDF code for  R by Thomas Lumley (c) 1999-2000
##
##  Based on code for S-PLUS by Gordon Maclean and others
#               Copyright (C) 1989,90,91,92,93,94 by UCAR

.First.lib<-function(libname,pkgname){
  library.dynam("netCDF",pkgname,libname)
}

is.open.netCDF<-function(x){
    .Call("is_open_netcdf",as.integer(x$id))>=0
}

if (!exists("close"))
    close<-function(x,...) UseMethod("close")

open.netCDF<-function(filename,verbose=F){
    if (!file.exists(filename))
        stop(paste("File",filename," not found"))
    rval<-.C("open_netcdf",filename,id=as.integer(verbose))
    rval$filename<-filename
    class(rval)<-"netCDF"
    return(rval);
}

close.netCDF<-function(x){
    if (.Call("is_open_netcdf",as.integer(x$id))<0)
        return(FALSE)
    else
        .C("close_netcdf",as.integer(x$id))
    return(TRUE)
}

print.netCDF<-function(x){
    cat("netCDF file",x$filename,"is ")
    if (!is.open.netCDF(x))
        cat("closed\n")
    else
        cat("open\n")
    x
}

names.netCDF<-function(x){
    rval<-.Call("do_netcdf_varnames",as.integer(x$id))
    names(rval)<-c("names","short.names","long.names","ids")
    rval
}

dim.netCDF<-function(x){
    .Call("do_dim_netcdf",as.integer(x$id))
}

summary.netCDF<-function(x){
    if  (.Call("is_open_netcdf",as.integer(x$id))<0)
        return(x)
    nn<-names(x)
    info<-.Call("do_inquire_varid",as.integer(x$id),as.integer(nn$ids))
    names(info)<-c("types","dims","dimids")
    info$names<-nn$names
    rval<-list(x,info)
    att<-.Call("do_get_global_atts",as.integer(x$id))
    
    for(i in seq(along=att)){
        attr(rval,names(att)[i])<-att[[i]]
    }
    class(rval)<-"summary.netCDF"
    rval
}

print.summary.netCDF<-function(x){
    print(x[[1]])
    print(x[[2]])
}

read.netCDF<-function(x,name=NULL,id=NULL,start=NULL,count=NULL,byrow=T,attr=T){
    
    if (is.character(x)){
        y<-open.netCDF(x)  ##x is a filename
        x<-y
        on.exit(close(x))  
    }

    if (!is.open.netCDF(x))    ## if not open then reopen it.
        x<-open.netCDF(x$filename)
    
    if (is.null(name) & is.null(id))
        name<-names(x)[[1]]
    if (!is.null(name) & !is.null(id))
        stop("Can't give both variable names and ids")
    
    if (is.null(id)){
        if (!is.character(name))
            stop("Names must be names")
        id<- .Call("do_varname_2_id",as.integer(x$id),name)
        if (any(is.na(id)))
            stop("Can't find those variables")
    }

    if (!is.numeric(id) | any(is.na(id)))
        stop("id must be numeric")

    if (is.null(start) | is.null(count))
        info<-.Call("do_inquire_varid",as.integer(x$id),as.integer(id))

    rval<-vector("list",length(id))
    for(i in 1:length(id)){
        if (is.list(start))
            s<-start[[i]]
        else if (!is.null(start))
            s<-start
        else
            s<-0*info[[2]][[i]]
        if (is.list(count))
            n<-count[[i]]
        else if (!is.null(count))
            n<-count
        else
            n<-info[[2]][[i]]-s

        if (any(n==0)) { ##empty variables
            rval[[i]]<-array(1,dim=rev(n))
            next
        }
        
        rval[[i]]<-.Call("do_read_netcdf",as.integer(x$id),as.integer(id[i]),as.integer(rev(s)),as.integer(rev(n)))
        if (byrow){
            if (length(dim(rval[[i]])) == 2)
                rval[[i]]<-t(rval[[i]])
            else if (length(dim(rval[[i]]))>2)
                rval[[i]]<-aperm(rval[[i]],length(dim(rval[[i]])):1)	# Transpose >2 dims
        }
        if (attr){
            att<-.Call("do_get_attribute",as.integer(x$id),as.integer(id[i]))[[1]]
            for(j in seq(along=att)){
                attr(rval[[i]],names(att)[j])<-att[[j]]
            }
        }
    }

    if(is.null(name))
        names(rval)<-as.character(id)
    else
        names(rval)<-name

    if (length(rval)==1)
        rval<-rval[[1]]

    rval
}


