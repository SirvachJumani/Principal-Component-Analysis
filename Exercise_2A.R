##########################################
##                                      ##       
##       Multivariate Methods           ##
##                                      ##
##       Exercise Sheet 2A:             ##
##      Principal Component             ##
##            Analysis                  ## 
##                                      ##
##########################################



#-------------------------------------------------------------------------------------------------
##  Exercise  2.1A        
#-------------------------------------------------------------------------------------------------

## Exercise 2.1A a)

load("UKZ.Rda")
View(UKZ)

pairs(UKZ)
help("pairs")


# covariance-PCA:

?prcomp
prcomp(UKZ,center=TRUE,scale=FALSE,retx=FALSE)      # scale = FALSE for covariance-PCA
pca.cov <- prcomp(UKZ)               
summary(pca.cov) 

UKZ.c <- scale(UKZ, center=TRUE, scale=FALSE) # centering of the columns 

?svd
UKZ.c.svd <- svd(UKZ.c, nu = nrow(UKZ.c), nv = ncol(UKZ.c)) # singular value decomposition

UKZ.c.svd$d
dim(UKZ.c.svd$d)


UKZ.c.svd$d^2         # squared standard deviactions  
pca.cov$sdev^2
UKZ.c.svd$d^2/(nrow(UKZ.c) - 1)

# or:

UKZ.c.svd$d
pca.cov$sdev
1/sqrt(nrow(UKZ.c)-1) * UKZ.c.svd$d


dim(UKZ.c.svd$u)
length(UKZ.c.svd$d)
dim(UKZ.c.svd$v)

# Construct a diagonal matrix frm the singular values  

UKZ.c.svd.dMatrix <- matrix(0, nrow=100, ncol=7) 
UKZ.c.svd.dMatrix[1:7,] <- diag(UKZ.c.svd$d)
UKZ.c.svd.dMatrix


UKZ.c.svd$u%*%UKZ.c.svd.dMatrix%*%t(UKZ.c.svd$v)
head(UKZ.c.svd$u%*%UKZ.c.svd.dMatrix%*%t(UKZ.c.svd$v)) 
head(UKZ.c)

## Exercise 2.1A b)

apply(UKZ.c, 2 ,var) 

# variances a very different 

# correlaiton-PCA:
 
prcomp(UKZ, retx=FALSE, center=TRUE, scale=TRUE)  # scale = TRUE for  correlation-PCA

# calculation of the PCA with singular decompositon by hand 

UKZ.s <- scale(UKZ, center=TRUE, scale=TRUE)
UKZ.s.svd <- svd(UKZ.s)
UKZ.s.svd$d^2/(nrow(UKZ.s) - 1)
sqrt(UKZ.s.svd$d^2/(nrow(UKZ.s) - 1))

## Exercise 2.1A c)

# Variant 1: T = USigma'

Score1 <- UKZ.s.svd$u %*% diag(UKZ.s.svd$d)

# Variant 2: T = X(V')'

Score2 <- UKZ.s %*% UKZ.s.svd$v

# Variant 3: use of  prcomp

prcomp(UKZ, center=TRUE, scale=TRUE, retx=TRUE)$x

# comparison:

sum(abs(Score1-Score2))
sum(abs(Score1 - prcomp(UKZ, center=TRUE, scale=TRUE, retx=TRUE)$x))


#-------------------------------------------------------------------------------------------------
##  Exercise  2.2A
#-------------------------------------------------------------------------------------------------

## Exercise 2.2A a)

load("SPData.Rda")

head(SPData.raw)
str(SPData.raw)

SPData <- SPData.raw[-c(1:3, ncol(SPData.raw))]
rownames(SPData) <- SPData.raw[,2]

SPData <- SPData[complete.cases(SPData)==TRUE,]

SPData <- as.matrix(SPData)

head(SPData)
str(SPData)
View(SPData)

## Exercise 2.2A b)

apply(SPData, 2 ,var) 
SPData.pca <- prcomp(SPData, retx=TRUE, center=TRUE, scale = TRUE) 
summary(SPData.pca)

# Calculation of the PCA with singular value decompositon by hand 

SPData.s <- scale(SPData, center=TRUE, scale=TRUE)
SPData.s.svd <- svd(SPData.s)
sdev.s <- 1/sqrt(nrow(SPData.s)-1) * SPData.s.svd$d

head(SPData.s.svd$u%*%diag(SPData.s.svd$d)%*%t(SPData.s.svd$v))
head(SPData.s)

round(sdev.s - SPData.pca$sdev, 6)

## Exercise 2.2A c)

V1 <- SPData.s.svd$v
V2 <- SPData.pca$rotation

V1 - V2

T1 <- SPData.s %*% SPData.s.svd$v       # XV*'

T2 <- SPData.s %*% SPData.pca$rotation

T3 <- SPData.pca$x

table(T1 - T2)
table(T2 - T3)

## Exercise 2.2A d)

V.2 <- round(V1^2, 2)
V.2

# Which cummunalities will result if all PC are used?

rowSums(V.2)
rowSums(V1^2)

#  Which cummunalities will result if the first three PC are used 

rowSums(V.2[,1:3])
# rowSums(V1[,1:3]^2)


# the variable  Price.earnings should be till well explained with a reduced model
# Which PC have to be considered?

View(SPData)

V.2

## Exercise 2.2A e)

# Calculation of the eigenvalues from the covariance matrix of the scaled and centered data matrix 

SPData.s.cov <- cov(SPData.s)
SPData.s.cov.eig <- eigen(SPData.s.cov)
SPData.s.cov.eig

SPData.s.cov.eig$values

# eigen(cov(scale(SPData)))$values

# or:

prcomp(SPData, retx=TRUE, center=TRUE, scale=TRUE)$sdev^2

SPData.pca$sdev^2 

# Calculation of the relative proportion of the variance explained by the corresponding PC. 
# cumulating these proportions

# What is the variance of scaled and centred variables?

eigenVec <- SPData.pca$sdev^2

relprop <- round(eigenVec/sum(eigenVec),4)*100 
relprop

cumprop <- round(cumsum(eigenVec)/sum(eigenVec),4)*100
cumprop

# Matrix for the results:

rbind(eigenvalue=eigenVec, ProportionVar=relprop, expltotalVar=cumprop)

## Exercise 2.2A f)

# Which principal components are chosed based on the Kaiser-Guttman criteria?

eigenVec > 1

# draw a scree-plot 
# how is this plot interpret?

screeplot(SPData.pca, type="barplot")
screeplot(SPData.pca, type="lines")

# alternative:

plot(SPData.pca)
plot(SPData.pca, type="l")

# Note:

plot(SPData)
screeplot(SPData) # only worked for output objects of prcomp()


# you want to explain at least 85% of the total variance 
# How many PC have to be chosen?

cumsum(eigenVec)/sum(eigenVec)        # here: sum(eigenVec) = total variance 

cumsum(eigenVec)/sum(eigenVec)>0.85 

# Chose all PC which explain more than the mean vairance of the single components 

average <- mean(eigenVec) 

eigenVec > average

sum(as.numeric(eigenVec > average))

 





