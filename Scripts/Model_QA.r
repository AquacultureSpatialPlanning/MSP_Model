Raw_Impacts_CW <- do.call('cbind',readMat('~/MSP_Model/Scripts/CrowT0v1/Raw_Impacts.mat')$Raw.Impacts)

## Run Crow Code Version and compare results
# setwd(paste0(wkdir,'/MSP_Model/Scripts/CrowT0v1'))
if(readline('Run CW version? ') == T){
  system2('/Applications/MATLAB_R2016b.app/bin/matlab',
    args = c('-nodesktop','-noFigureWindows','-nosplash','-r',
    "run\\(\\'~/MSP_Model/Scripts/CrowT0v1/TOA_AquaMSP_CrowCode_v3NaN.m\\'\\);","exit"))
}else{
  # Loading CW results
  print('loading CW results')
  CW_Variables <- readMat('~/MSP_Model/Scripts/CrowT0v1/TOA_data.mat')
  print(str(CW_Variables))
  print('Primary model variables include the following........')
  CW_vars <- c('Raw_Impacts','R.bar.n','R.n.i.p','V.n.i.p','max.p.X.n.i.p','aMatrix','Policy_i_a')
  invisible(lapply(CW_vars, FUN = print))
  CW.list <- lapply(CW_vars[c(-1,-7)], FUN = function(x){CW_Variables[[x]]})
  print(str(CW.list))

  # Loading JS results
  print('loading JS results')
  load('~/MSP_Model/Output/Data/JS_MSP_Output.Rdata')
  print(str(JS_MSP.list))
}
# Case study as outlined by CW on 4/20
i.test <- 138
a.test <- 55988
test.vars <- function(a,b,str){
  print(paste0(str,': '))
  invisible(lapply(paste0(unique(ifelse(unique(a == b) == T,'YES','NO'))), FUN = print))
}
print(paste('Evaluate test case......'))
print(paste('Test site # = ',i.test))
print(paste('Alpha scenario = ',a.test))
# Test 1 compare alpha values
a.JS <- JS_MSP.list[['a']][a.test,-6]
a.CW <- CW_Variables[['aMatrix']][a.test,]
test.vars(a.JS,a.CW,'Are alpha values identical?')
# Test 2 compare Raw Impacts
Raw_Impacts.JS <- JS_MSP.list[['Raw_Impacts']] %>% slice(i.test) %>% glimpse()
Raw_Impacts.CW <- setNames(data.frame(CW_Variables[['NUM1']],stringsAsFactors = F),c(names(Raw_Impacts.JS))) %>% slice(i.test) %>% glimpse()
test.vars(Raw_Impacts.JS,Raw_Impacts.CW,'Are the sector values for the test scenario identical?')
# Test 3 compare the R_n_i_p
print('JS Response data.....')
print(str(JS_MSP.list[['R_n_i_p']]))
print('CW Response data.....')
print(dim(CW_Variables[['R.n.i.p']]))

R_n_i_p.JS <- lapply(JS_MSP.list[['R_n_i_p']],"[", i.test,)
print(str(R_n_i_p.JS))
# Convert CW output to a comparable f
R_n_i_p.CW <- setNames(lapply(1:ncol(CW_Variables[['R.n.i.p']][i.test, ,]), FUN = function(x){
    setNames(data.frame(t(CW_Variables[['R.n.i.p']][i.test, ,x]),stringsAsFactors = F), c(names(R_n_i_p.JS[[1]])[1:(which(names(R_n_i_p.JS[[1]]) %in% names(select(R_n_i_p.JS[[1]],contains('Viewshed_'))))[1] - 1)],
      'Viewshed',
      names(R_n_i_p.JS[[1]])[which(names(R_n_i_p.JS[[1]]) %in% names(select(R_n_i_p.JS[[1]],contains('Viewshed_'))))[2] + 1], names(R_n_i_p.JS[[1]])[length(names(R_n_i_p.JS[[1]]))]))
  }),names(R_n_i_p.JS))
print(str(R_n_i_p.CW))
# Come up with a better test later
# Test 4 compare X_n_i_p
X_n_i_p.JS <- unname(sapply(JS_MSP.list[['V_n_i_p']],"[", i.test, ))
row.names(X_n_i_p.JS) <- c();
print(X_n_i_p.JS)
X_n_i_p.CW <- CW_Variables[['X.n.i.p']][i.test, ,]
print(X_n_i_p.CW)
# Test 1 compare Alpha sets
# Because of this, columns 5 and 6 in JS a-matrix are identical to one another. These need to therefore be removed in order to do an accurate comparison.
# Note: JS version uses 8 sectors (in a quasi sort of way). This is due to how the model is set up with the two viewsheds which are calculated seperately from one another.
a.JS <- JS_MSP.list[['a']][,-6]
print(paste('JS alpha matrix = ',dim(a.JS)))
a.CW <- CW_Variables[['aMatrix']]
print(paste('CW alpha matrix = ',dim(a.CW)))
# Take the absolute difference between the two
# apply(a.JS,MARGIN = 1,FUN = function(x){as.character((factor(x)})))
# apply(a.CW,MARGIN = 1,FUN = sum)
# as.character(factor(a.JS))
# Take the
# Test 1 Compare Raw_Impacts, i.e. the sector inputs to the two versions
Raw_Impacts.CW <- CW_Variables[['NUM1']]
