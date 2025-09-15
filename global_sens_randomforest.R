# Global Sensitivity Analysis of Commander et al. oyster predator-prey model
# Using procedures described in Harper et al. 2011 Ecological Applications

#install.packages('rpart')
#install.packages('randomForest')

library('rpart')
library('randomForest') #Liaw and Wiener 

#A. Liaw and M. Wiener (2002). Classification and Regression by #randomForest. R News 2(3), 18--22.

# load in data.  22Dec runs had 1e5 simulations
D1 = read.csv('Global_sens_8Sept2025.csv',head=FALSE)

# The first column is the average population density

# Columns 2:17 are the parameter values
for (i in 2:17){
  D1[,i]= (mean(D1[,i])-D1[,i])/sd(D1[,i])
}

#D1 = D1[1:2000,]



# D.impt gives the relative importance of each parameter

#Insert variable names to make it easier to interpret
Parm_names = c('Prey.k','Prey.Linf','Prey.Mj','Prey.M',
               'Prey.Mat','Prey.Fec.coeff',
               'Prey.DDa','Prey.DDb','Prey.b2',
               'Prey.b4','Prey.b5a','Prey.b4a','pred.M',
               'pred.fec.coeff','pred.LEP','pred.aP')

Parm_names2 = c('X',Parm_names)

colnames(D1) = Parm_names2


# Perform randomForests
D1.rf = randomForest(X~.,data=D1,ntree=1000,importance=TRUE)
D1.impt = importance(D1.rf,type=1) # type 1 is mean decrease in accuracy


rownames(D1.impt) = Parm_names



barplot(abs(D1.impt)/sum(abs(D1.impt)),beside=TRUE,names.arg=rownames(D1.impt))


partialPlot(D1.rf, D1, Prey.b2, xlab="Larval salinity parameter", ylab="Oyster population density", lwd=4, col="green")



