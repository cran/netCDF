
/* 
 * C routines callable by R to read netcdf files.
 *
 */

/*
 * Copyright (c) 1999 Thomas Lumley
 * S+ version by Gordon Maclean
 * Initial port to R by Steven Oncley
 * 
 * Portions of this code copyright 1998 University Corporation for
 * Atmospheric Research (UCAR)/ Atmospheric Technology Division 
 *
 * For complete license terms see file LICENSE
 *
 * The user is granted the right, without any fee or cost, to use, copy, 
 * modify and distribute this software, and any derivative works thereof, for 
 * any purpose whatsoever, provided that this entire notice appears in all 
 * copies of the software.  UCAR is not obligated to provide the user with 
 * any updates or support of this software.
 *
 *
 * THIS SOFTWARE IS PROVIDED "AS IS" AND ANY IMPLIED WARRANTIES ARE DISCLAIMED.
 */

/** TODO:                    **/
/**       update to netCDF 3 **/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "netcdf.h"

#include <Rinternals.h>
#include <Rdefines.h>
#include <S.h>


SEXP is_open_netcdf(SEXP ncid){
    SEXP rval;
    int ec,oldncopt;
    
    oldncopt=ncopts;
    ncopts=0;

    ec=ncinquire(INTEGER(ncid)[0],(int *)0, (int *)0,(int *)0, (int *)0);

    ncopts=oldncopt;

    PROTECT(rval=allocVector(INTSXP,1));
    INTEGER(rval)[0]=ec;
    UNPROTECT(1);
    return rval;
}
	

void open_netcdf(
  char **filename,	/* NetCDF filename, passed by S function */
  int *ncid_r)		/* input: verbose error option
			   output:NetCDF file id */
{

  int  ncid;			/* netCDF id */
  /* is this a horrible way to set options or what?*/
  ncopts=0;
  if (*ncid_r)
      ncopts = NC_VERBOSE;

  /* open cdf file */
  if ((ncid = ncopen(R_ExpandFileName(filename[0]), NC_NOWRITE)) < 0) 
      warning("netCDF error opening file");
  *ncid_r = ncid;

  return;

}
void close_netcdf(
  int *ncid_p)		/* NetCDF file id */
{
  ncclose(*ncid_p);
  return;
}

SEXP do_varname_2_id(SEXP nc_id,SEXP namelist){
    SEXP  rval;
    
    int ncid,i,nn;

    /* variable ids */
    int  var_id;

    ncid=INTEGER(nc_id)[0];
    
    nn=length(namelist);
    PROTECT(rval=allocVector(INTSXP,nn));

    for (i=0;i<nn;i++){
	var_id=ncvarid(ncid,CHAR(STRING(namelist)[i]));
	if (var_id<0)
	    INTEGER(rval)[i]=NA_INTEGER;
	else 
	    INTEGER(rval)[i]=var_id;
    }
    UNPROTECT(1);
  return rval;
}



SEXP do_inquire_varid(SEXP ncid, SEXP var_id){

    SEXP dimidlist,dimlist,vartypes,rval;
    int i,errorcode,j;
    nc_type vtype;

    /* variable shapes */
    int dims[MAX_VAR_DIMS];
    long dimval;
    
    int ndims;
    
    int nvars=length(var_id);
    
    PROTECT(dimlist=allocVector(VECSXP,nvars));
    PROTECT(dimidlist=allocVector(VECSXP,nvars));
    PROTECT(vartypes=allocVector(STRSXP,nvars));

    for (i=0; i<nvars;i++){
	errorcode=ncvarinq(INTEGER(ncid)[0],INTEGER(var_id)[i],(char *)0,
			   &vtype, &ndims,dims,(int *)0);
	if (errorcode<0) {
	    warning("Nonexistent netCDF variable ID");
	    VECTOR(dimidlist)[i]=allocVector(INTSXP,1);
	    INTEGER(VECTOR(dimidlist)[i])[0]=NA_INTEGER;
	    VECTOR(dimlist)[i]=allocVector(INTSXP,1);
	    INTEGER(VECTOR(dimlist)[i])[0]=NA_INTEGER;
	    STRING(vartypes)[i]=NA_STRING;
	    continue;
	} else {
	    VECTOR(dimidlist)[i]=allocVector(INTSXP,ndims);
	    for (j=0;j<ndims;j++){
		INTEGER(VECTOR(dimidlist)[i])[j]=dims[j];
	    }
	    switch (vtype) {
	    case NC_LONG:
		STRING(vartypes)[i]=mkChar("integer");
		break;
	    case NC_FLOAT:
		STRING(vartypes)[i]=mkChar("single");
		break;
	    case NC_DOUBLE:
		STRING(vartypes)[i]=mkChar("numeric");
		break;
	    case NC_CHAR:
		STRING(vartypes)[i]=mkChar("char");
		break;
	    case NC_SHORT:
		STRING(vartypes)[i]=mkChar("short integer");
		break;
	    case NC_BYTE:
		STRING(vartypes)[i]=mkChar("byte");
		break;
	    default:
		warning("Unknown variable type");
		STRING(vartypes)[i]=NA_STRING;
	    }
	}	    
	VECTOR(dimlist)[i]=allocVector(INTSXP,ndims);
	for (j=0;j<ndims;j++){
	    errorcode=ncdiminq(INTEGER(ncid)[0],dims[j],(char*)0, &dimval);
	    INTEGER(VECTOR(dimlist)[i])[j]=dimval;
	}

    }
    PROTECT(rval=allocVector(VECSXP,3));
    VECTOR(rval)[0]=vartypes;
    VECTOR(rval)[1]=dimlist;
    VECTOR(rval)[2]=dimidlist;
    UNPROTECT(4);
    return rval;
}


SEXP do_read_netcdf(SEXP ncid, SEXP var_id, SEXP start,SEXP count){

    SEXP rval;
    
    int ndims,ndims0,ec;
    nc_type vtype;
    int i;
    int size;
    long *startl,*countl;

    float *df; short int* ds; char *dc;
    
    size=1;
    ndims=length(count);
    if (ndims!=length(start)) 
	error("Start and count do not match");
    for(i=0;i<ndims;i++){
	size *=INTEGER(count)[i];
    }

    ncvarinq(INTEGER(ncid)[0],INTEGER(var_id)[0],(char *)0,&vtype,&ndims0,(int *)0,(int *)0);
    if (ndims0!=ndims)
	error("Dimensions do not match file");
    
    startl=(long *) Calloc(ndims,long);
    countl=(long *) Calloc(ndims,long);
    for(i=0;i<ndims;i++){  /* S dimensions are backwards */
	startl[ndims-i-1]=(long) (INTEGER(start)[i]);
	countl[ndims-i-1]=(long) (INTEGER(count)[i]);
    }

    rval=R_UnboundValue; /*-Wall*/
    switch (vtype) {
    case NC_LONG:
	PROTECT(rval=allocVector(INTSXP,size));
	ec=ncvarget(INTEGER(ncid)[0],INTEGER(var_id)[0],
		 startl,countl, INTEGER(rval));
	break;
    case NC_FLOAT:
	PROTECT(rval=allocVector(REALSXP,size));
	df=Calloc(size,float);
	ec=ncvarget(INTEGER(ncid)[0],INTEGER(var_id)[0],
		 startl,countl, df);
	for (i = 0; i < size; i++) REAL(rval)[i] = (double) df[i];
	Free(df);
	break;
    case NC_DOUBLE:
	PROTECT(rval=allocVector(REALSXP,size));
	ec=ncvarget(INTEGER(ncid)[0],INTEGER(var_id)[0],
		 startl,countl, REAL(rval));
	break;
    case NC_CHAR:
	error("Can't handle NC_CHAR yet");
	break;
    case NC_SHORT:
	PROTECT(rval=allocVector(INTSXP,size));
	ds=Calloc(size,short int);
	ec=ncvarget(INTEGER(ncid)[0],INTEGER(var_id)[0],
		 startl,countl, ds);
	for (i = 0; i < size; i++) INTEGER(rval)[i] = (int) ds[i];
	Free(ds);
	break;
    case NC_BYTE: /** this may not be the Right Thing to do with bytes **/
	PROTECT(rval=allocVector(INTSXP,size));
	dc=Calloc(size,char);
	ec=ncvarget(INTEGER(ncid)[0],INTEGER(var_id)[0],
		 startl,countl, dc);
	for (i = 0; i < size; i++) INTEGER(rval)[i] = (int) dc[i];
	Free(dc);
	break;
    default:
	error("Unknown variable type");
    }
    
    if (ec<0) {
	rval=R_NilValue;
	error("Error in reading");
    }

    /* fix up dimensions */
    setAttrib(rval,R_DimSymbol,count);

    UNPROTECT(1);
    return rval;
}

SEXP do_get_global_atts(SEXP ncid){
    SEXP attname, attributes;
    int j,k,errorcode,var_id,attlen,natts;
    nc_type atype;
    float *df;short int *ds; char *dc;
    char name[MAX_NC_NAME];

    errorcode=ncinquire(INTEGER(ncid)[0],(int *)0,(int *)0,&natts,(int *)0);

    if ((errorcode<0) || (natts==0))
	return R_NilValue;

    PROTECT(attributes=allocVector(VECSXP,natts));
    PROTECT(attname=allocVector(STRSXP,natts));
 
    for(k=0;k<natts;k++){
	    errorcode=ncattname(INTEGER(ncid)[0],NC_GLOBAL,k,name);
	    if (errorcode<0){
		STRING(attname)[k]=NA_STRING;
		VECTOR(attributes)[k]=R_NilValue;
		continue;
	    }
	    STRING(attname)[k]=mkChar(name);
	    errorcode=ncattinq(INTEGER(ncid)[0],NC_GLOBAL,
			   name,&atype,&attlen);

	    if (errorcode<0){
		VECTOR(attributes)[k]=R_NilValue;
		continue;
	    }
	    
	    switch (atype) {
	    case NC_LONG:
		VECTOR(attributes)[k]=allocVector(INTSXP,attlen);
		if (ncattget(INTEGER(ncid)[0],var_id,name,
			     (void *)(INTEGER(VECTOR(attributes)[k]))) <0)
		    VECTOR(attributes)[k]=R_NilValue;
		break;
	    case NC_FLOAT:
		df=Calloc(attlen,float);
		if (ncattget(INTEGER(ncid)[0],var_id,
			     name,(void *)df) <0)
		    VECTOR(attributes)[k]=R_NilValue;
		else {
		    VECTOR(attributes)[k]=allocVector(REALSXP,attlen);
		    for(j=0;j<attlen;j++)
			REAL(VECTOR(attributes)[k])[j]=(double) df[j];
		}
		Free(df);
		break;
	    case NC_DOUBLE:
		VECTOR(attributes)[k]=allocVector(REALSXP,attlen);
		if (ncattget(INTEGER(ncid)[0],var_id,name,
			     (void *)(REAL(VECTOR(attributes)[k])))<0)
		    VECTOR(attributes)[k]=R_NilValue;
		break;
	    case NC_CHAR:
		dc=Calloc(attlen+1,char);
		if (ncattget(INTEGER(ncid)[0],var_id,name,
			     (void *)dc)<0)
		    VECTOR(attributes)[k]=R_NilValue;
		else {
		dc[attlen]=0; /* may not be zero-terminated */
		VECTOR(attributes)[k]=allocVector(STRSXP,1);
		STRING(VECTOR(attributes)[k])[0]=mkChar(dc);
		}
		Free(dc);
		break;
	    case NC_SHORT:
		ds=Calloc(attlen,short int);
		if (ncattget(INTEGER(ncid)[0],var_id,name,
			(void *)ds)<0)
		    VECTOR(attributes)[k]=R_NilValue;
		else {
		    VECTOR(attributes)[k]=allocVector(INTSXP,attlen);
		    for(j=0;j<attlen;j++)
			INTEGER(VECTOR(attributes)[k])[j]=(int) ds[j];
		}
		Free(ds);
		break;
	    case NC_BYTE:
		VECTOR(attributes)[k]=R_NilValue;
		warning("Can't handle NC_BYTE attributes");
		break;
	    default:
		warning("Unknown attribute type");
		VECTOR(attributes)[k]=R_NilValue;
	    }

    }
    setAttrib(attributes,R_NamesSymbol,attname);
    UNPROTECT(2);
    
    return(attributes);
}

SEXP do_get_attribute(SEXP ncid, SEXP var_ids){

    SEXP attname, attributes,v;
    int i,j,k,nvars,errorcode,var_id,attlen,natts;
    nc_type atype;
    float *df;short int *ds; char *dc;
    char name[MAX_NC_NAME];

    nvars=length(var_ids);
    PROTECT(v=allocVector(VECSXP,nvars));
    for(i=0;i<nvars;i++){
	var_id=INTEGER(var_ids)[i];
	errorcode=ncvarinq(INTEGER(ncid)[0],var_id,(char *)0,(nc_type *)0,
			   (int*)0,(int *)0,&natts);
	if ((errorcode<0) || (natts==0)){
	    VECTOR(v)[i]=R_NilValue;
	    continue;
	}
	PROTECT(attname=allocVector(STRSXP,natts));
	PROTECT(attributes=allocVector(VECSXP,natts));
	for(k=0;k<natts;k++){
	    errorcode=ncattname(INTEGER(ncid)[0],var_id,k,name);
	    if (errorcode<0){
		STRING(attname)[k]=NA_STRING;
		VECTOR(attributes)[k]=R_NilValue;
		continue;
	    }
	    STRING(attname)[k]=mkChar(name);
	    errorcode=ncattinq(INTEGER(ncid)[0],var_id,
			   name,&atype,&attlen);

	    if (errorcode<0){
		VECTOR(attributes)[k]=R_NilValue;
		continue;
	    }
	    
	    
	    switch (atype) {
	    case NC_LONG:
		VECTOR(attributes)[k]=allocVector(INTSXP,attlen);
		if (ncattget(INTEGER(ncid)[0],var_id,name,
			     (void *)(INTEGER(VECTOR(attributes)[k]))) <0)
		    VECTOR(attributes)[k]=R_NilValue;
		break;
	    case NC_FLOAT:
		df=Calloc(attlen,float);
		if (ncattget(INTEGER(ncid)[0],var_id,
			     name,(void *)df) <0)
		    VECTOR(attributes)[k]=R_NilValue;
		else {
		    VECTOR(attributes)[k]=allocVector(REALSXP,attlen);
		    for(j=0;j<attlen;j++)
			REAL(VECTOR(attributes)[k])[j]=(double) df[j];
		}
		Free(df);
		break;
	    case NC_DOUBLE:
		VECTOR(attributes)[k]=allocVector(REALSXP,attlen);
		if (ncattget(INTEGER(ncid)[0],var_id,name,
			     (void *)(REAL(VECTOR(attributes)[k])))<0)
		    VECTOR(attributes)[k]=R_NilValue;
		break;
	    case NC_CHAR:
		dc=Calloc(attlen+1,char);
		if (ncattget(INTEGER(ncid)[0],var_id,name,
			     (void *)dc)<0)
		    VECTOR(attributes)[k]=R_NilValue;
		else {
		dc[attlen]=0; /* may not be zero-terminated */
		VECTOR(attributes)[k]=allocVector(STRSXP,1);
		STRING(VECTOR(attributes)[k])[0]=mkChar(dc);
		}
		Free(dc);
		break;
	    case NC_SHORT:
		ds=Calloc(attlen,short int);
		if (ncattget(INTEGER(ncid)[0],var_id,name,
			(void *)ds)<0)
		    VECTOR(attributes)[k]=R_NilValue;
		else {
		    VECTOR(attributes)[k]=allocVector(INTSXP,attlen);
		    for(j=0;j<attlen;j++)
			INTEGER(VECTOR(attributes)[k])[j]=(int) ds[j];
		}
		Free(ds);
		break;
	    case NC_BYTE:
		VECTOR(attributes)[k]=R_NilValue;
		warning("Can't handle NC_BYTE attributes");
		break;
	    default:
		warning("Unknown attribute type");
		VECTOR(attributes)[k]=R_NilValue;
	    }
	}
	setAttrib(attributes,R_NamesSymbol,attname);
	VECTOR(v)[i]=attributes;
	UNPROTECT(2);
    }
    UNPROTECT(1); /*v*/
    return v;
}




SEXP do_netcdf_varnames(SEXP ncid){

    SEXP varnames,short_names,long_names,ids,rval;
    int i,l,nvars,ec,vid,oldncopt;
    char *vname;
    char name[MAX_NC_NAME];
    nc_type atttype;
    
    oldncopt=ncopts;
    ncopts=0;  /* don't want errors when name attributes are missing */

    ncinquire(INTEGER(ncid)[0],(int *)0, &nvars,(int *)0, (int *)0);
    
    
    PROTECT(varnames=allocVector(STRSXP,nvars));
    PROTECT(short_names=allocVector(STRSXP,nvars));
    PROTECT(long_names=allocVector(STRSXP,nvars));
    PROTECT(ids=allocVector(INTSXP,nvars));

    for (i=0;i<nvars;i++){
	vid=i;
	INTEGER(ids)[i]=vid;
	ec=ncattinq(INTEGER(ncid)[0],vid,"short_name",(nc_type *) &atttype,&l);    
	if ((ec<0) || (l<=0) || (atttype!=NC_CHAR))
	    STRING(short_names)[i]=NA_STRING;
	else{  
	    vname=Calloc(l+1,char);
	    ncattget(INTEGER(ncid)[0],vid,"short_name",(void *) vname);
	    vname[l]='\0';
	    STRING(short_names)[i]=mkChar(vname);
	    Free(vname);
	}
	ec=ncattinq(INTEGER(ncid)[0],vid,"long_name",(nc_type *) &atttype,&l);    
	if ((ec<0) || (l<=0) || (atttype!=NC_CHAR))
	    STRING(long_names)[i]=NA_STRING;
	else{  
	    vname=Calloc(l+1,char);
	    ncattget(INTEGER(ncid)[0],vid,"long_name",(void *) vname);
	    vname[l]='\0';
	    STRING(long_names)[i]=mkChar(vname);
	    Free(vname);
	}
	ncvarinq(INTEGER(ncid)[0],vid,name,(nc_type*)0,(int *)0,(int *)0,(int *)0);
	STRING(varnames)[i]=mkChar(name);
    }

    ncopts=oldncopt;
    
    PROTECT(rval=allocVector(VECSXP,4));
    VECTOR(rval)[0]=varnames;
    VECTOR(rval)[1]=short_names;
    VECTOR(rval)[2]=long_names;
    VECTOR(rval)[3]=ids;
    UNPROTECT(4);

    UNPROTECT(1); /*rval*/
    return rval;
}


SEXP do_dim_netcdf (SEXP ncid){

    int i,ndims,errorcode;
    long size;
    SEXP dims,dimnames;
    char name[MAX_NC_NAME];

    errorcode=ncinquire(INTEGER(ncid)[0],&ndims, (int *)0,(int *)0, (int *)0);
    
    if (errorcode<0) 
	return R_NilValue;

    PROTECT(dims=allocVector(INTSXP,ndims));
    PROTECT(dimnames=allocVector(STRSXP,ndims));

    for (i=0;i<ndims;i++){ 
	errorcode=ncdiminq(INTEGER(ncid)[0],i,name,&size);
	if (errorcode<0) {
	    INTEGER(dims)[i]=errorcode;
	    STRING(dimnames)[i]=NA_STRING;
	    continue;
	}
	INTEGER(dims)[i]=(int) size;
	STRING(dimnames)[i]=mkChar(name);
    }
    setAttrib(dims,R_NamesSymbol,dimnames);
    UNPROTECT(2);
    return dims;
}

