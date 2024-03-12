##Geographical Temporal Quantile Regression=name
##GWmodel=group
##Layer=vector
##Dependent_Variable=Field Layer
##QgsProcessingParameterField|explana_variables|Expanatory Variables|None|Layer|-1|True|False
##GeoId=Field Layer
##QgsProcessingParameterNumber|minhIn|min H|1
##QgsProcessingParameterNumber|maxhIn|Max H|1
##QgsProcessingParameterNumber|tauIn|Tau|1|0.5|False|0|1
##kernel=selection global;gaussian;exponential;bisquare global
##QgsProcessingParameterBoolean|diskmIn|Great Circle Distance|False
##QgsProcessingParameterBoolean|locallinIn|Local Linear
##Output=output vector

# libraries
library(sp)
library(sf)
library(quantreg)
library(MatrixModels) ## 缺少
library(MASS)
library(parallel)
library(snow)
library(stringr)

#---------------- Verify User's Input ---------------

if (tauIn > 1) stop("tau must be less than 1")
if (minhIn > maxhIn) stop("Min_H must less than Max_H")

#----------------- define functions -----------------

cvgwqr <- function(formula,tau=0.5,cordx,cordy,minh,maxh,kernel,indata,distkm=F,locallin=F){

    time1=Sys.time()
    #Set CV function
    cvh<-function(h){
        mat<-model.frame(formula, data = indata)
        y=mat[, 1]
        ox=as.matrix(model.matrix(formula, data = indata))[,-1]  #designed x matrix excluding intercept (1's)
        colnames(ox)<-names(mat)[-1]
        nobs=length(y)
        x=cbind(rep(1,nobs),ox)
        nvar=ncol(x)

        beta=c(); muhat=c()
        for (i in 1:nobs){
            xstart=cordx[i]
            ystart=cordy[i]
            xd=cordx-xstart
            yd=cordy-ystart
            d<-sqrt(xd**2+yd**2)
            if (distkm == T) d = geosphere::distCosine(cbind(cordx,cordy),cbind(xstart,ystart),r=6371)  #great circle distances
            rh=round(h)
            sd=sort(d)
            disth=sd[rh]
            cid=which(d!=0)
            cy=y[cid]
            cx=ox[cid,]
            dd=d[cid]

        
            if (kernel=="global") {
                w=rep(1,nobs-1)
            }
            if (kernel=="gaussian"){
                w=exp((-0.5)*((dd/disth)**2))
                w=as.numeric(w)
            }# /*gaussian kernel*/
            
            if (kernel=="exponential"){
                w=exp(-dd/disth)
                w=as.numeric(w)
            } #/*exponential kernel*/
            if (kernel=="bisquare"){
                #w=rep(0,nr)
                w=(1-(dd/disth)^2)^2
                index=which(dd>disth)
                w[index]=0
            } #/*bisquare nearest neighbor*/

            newdata = data.frame(cx)
            if (locallin==T) {  #for local linear estimation
                nvar1=3*nvar
                xdz=xd*ox; ydz=yd*ox
                allx = cbind(ox,xd,yd,xdz,ydz)
                callx = allx[cid,]
                newdata = data.frame(callx) }
                fit <-rq(cy~.,tau=tau,weights=w,method="fn",data=newdata)
                fit1 <- summary(fit)
                bmat <- coef(fit1)[1:nvar,1]
                #print(bmat)

                beta=rbind(beta,bmat)
                yhat = x[i,]%*%as.matrix(bmat[1:nvar])
                muhat =rbind(muhat,yhat)
        }
        error=y-muhat
        wsum=sum((error*(tau-ifelse(error<0,1,0))))
        out=list(y=y,x=x,beta=beta,muhat=muhat,error=error,wsum=wsum)  
        return(out)
    }
    # /*Golden Section Serach*/
    eps=1;
    r=(sqrt(5)-1)/2
    a0=minh
    b0=maxh
    x1=r*a0+ (1-r)*b0
    x2=(1-r)*a0+ r*b0
    fx1=cvh(x1)$wsum
    #print(cbind(fx1))
    fx2=cvh(x2)$wsum
    
    it=0;
    while ((b0-a0) > eps){
        it=it+1
        if (fx1<fx2){
          b0=x2
          x2=x1
          fx2=fx1
          x1=r*a0+(1-r)*b0
          fx1=cvh(x1)$wsum
        }else{
          a0=x1
          x1=x2
          fx1=fx2
          x2=(1-r)*a0+r*b0
          fx2=cvh(x2)$wsum
        }
    }
    h=ifelse(fx1<=fx2,a0,b0) #*將fx1小於fx2的b0改成a0;
    print(paste(time1,Sys.time()))
    print(paste(a0,b0))
    names(h)=c("lambda")
    return(h)
}


gwqr <- function(formula,tau=0.5,geoid,cordx,cordy,h,kernel,indata,distkm=F,locallin=F){
  mat<-model.frame(formula, data = indata)
  y=mat[, 1]
  x=as.matrix(model.matrix(formula, data = indata))  #designed x matrix including intercept (1's)
  ox=x[,-1]  #designed x matrix excluding intercept (1's)
  nvar=ncol(x)
  nobs=length(y)
 
  beta=c(); muhat=c(); std=c(); tvalue=c()
  for (i in 1:nobs){
    xstart=cordx[i]
    ystart=cordy[i]
    xd=cordx-xstart
    yd=cordy-ystart
    d<-sqrt(xd**2+yd**2)
    if (distkm == T) d = geosphere::distCosine(cbind(cordx,cordy),cbind(xstart,ystart),r=6371)  #great circle distances
    rh=round(h)
    sd=sort(d)
    disth=sd[rh]
    # kernel function
    if (kernel=="global") w=rep(1,nobs)
    if (kernel=="gaussian"){
      w=exp((-0.5)*((d/disth)**2))
      w=as.numeric(w)}# /*gaussian kernel*/
    
    if (kernel=="exponential"){
      w=exp(-d/disth)
      w=as.numeric(w)} #/*exponential kernel*/
    if (kernel=="bisquare"){
      #w=rep(0,nobs)
      w=(1-(d/disth)^2)^2
      index=which(d>disth)
      w[index]=0}   #/*bisquare nearest neighbor*/
    
    newdata = data.frame(ox)
    if (locallin==T) {      #for local linear estimation
       nvar1=3*nvar
       xdz=xd*ox; ydz=yd*ox
       allx = cbind(ox,xd,yd,xdz,ydz)
       newdata = data.frame(allx) }
    fit <-rq(y~.,tau=tau,weights=w,method="fn",data=newdata)
    fit1 <- suppressWarnings(summary(fit,se="nid"))
    bmat <- coef(fit1)[1:nvar,1]
    #print(bmat)
    semat <- coef(fit1)[1:nvar,2]
    #print(semat)
    tmat <- coef(fit1)[1:nvar,3]
  
    std=rbind(std,semat)
    tvalue=rbind(tvalue,tmat)
    yhat = x[i,]%*%as.matrix(bmat[1:nvar])
    beta=rbind(beta,bmat)
    muhat=rbind(muhat,yhat)
    #print(muhat)
  }
  error=y-muhat
  #loss function
  wsum=sum((error*(tau-ifelse(error<0,1,0))))
  out=list(y=y,x=x,geoid=geoid,cordx=cordx,cordy=cordy,beta=beta,std=std,tvalue=tvalue,muhat=muhat,error=error,wsum=wsum)  
  return(out)
}


#--------------------- Input ---------------------
kernel_list <- c("global", "gaussian", "exponential", "bisquare")
kernel_choice <- kernel_list[kernel + 1]

# rename selection
# define shorten_names function
shorten_names <- function(names) {
  name_dict <- list()
  new_names <- character(length(names))
  
  for (i in seq_along(names)) {
    name <- names[i]
    short_name <- substr(name, 1, 10)
    
    if (!is.null(name_dict[[short_name]])) {
      count <- name_dict[[short_name]] + 1
    } else {
      count <- 1
    }
    
    if (count > 1) {
      short_name <- substr(short_name, 1, 8)
      short_name <- paste0(short_name, "_", count - 1)
    }
    
    name_dict[[short_name]] <- count
    new_names[i] <- short_name
  }
  
  dupes <- duplicated(new_names)
  if (any(dupes)) {
    for (i in which(dupes)) {
      base_name <- substr(new_names[i], 1, 8)
      count <- 2
      new_name <- paste0(base_name, "_", count)
      
      while (new_name %in% new_names) {
        count <- count + 1
        new_name <- paste0(base_name, "_", count)
      }
      
      new_names[i] <- new_name
    }
  }
  
  return(new_names)
}
explana_variables = shorten_names(explana_variables)
Dependent_Variable = shorten_names(Dependent_Variable)

# create formula by user's input 
formula_str = paste0(Dependent_Variable, " ~ ")
for (x in 1:length(explana_variables)) {
  if (x != length(explana_variables)) {
    formula_str = paste0(formula_str, explana_variables[x], " + ")
  } else {
    formula_str = paste0(formula_str, explana_variables[x])
  }
}
formula = formula(formula_str)

data = na.omit(Layer)
data = as.data.frame(data) # convert sf_data.frame to dataframe
coords_from_sfdf = st_coordinates(Layer)
coord_X = coords_from_sfdf[, "X"]
coord_Y = coords_from_sfdf[, "Y"]

#--------------------- model ---------------------

bw <- cvgwqr(formula, tau = tauIn, cordx = coord_X, cordy = coord_Y,
    minh = minhIn, maxh = maxhIn, kernel = kernel_choice, indata = data,
    distkm = diskmIn, locallin = locallinIn)
model <- gwqr(formula, tau = tauIn, geoid = data[GeoId],
    cordx = coord_X, cordy = coord_Y, h = round(bw), kernel = kernel_choice,
    indata = data, distkm = diskmIn, locallin = locallinIn)

# -------------------- result --------------------

res_df = data.frame(geoId = model$geoid, Y = model$y)
res_df = cbind(res_df, as.data.frame(model$x[,2:ncol(model$x)]))
temp = data.frame(CoordX = model$cordx, CoordY = model$cordy)
res_df = cbind(res_df, temp)
var = c("beta", "std", "tvalue", "muhat", "error")
for (v in var) {
  temp = model[[v]]
  name = ifelse(is.null(colnames(temp)), v, paste0(v, "_"))
  colnames(temp) = paste0(name, colnames(temp))
  res_df = cbind(res_df, temp)
}  


LayerCRS = st_crs(Layer)$epsg
sp_res_df = SpatialPointsDataFrame(cbind(model$cordx, model$cordy), res_df) # convert to SpatialPointsDataFrame 
Output=st_as_sf(sp_res_df)