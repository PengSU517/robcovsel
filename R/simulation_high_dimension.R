library(doParallel)
library(robcovsel)
registerDoParallel(cores=20)
getDoParWorkers()


###### data generation settings
{
  ms = 1:10
  ps = 200 # ps = 200 in high dimensional settings
  es = c(0.02, 0.05, 0.1)
  rs = c(0.5, 0.7, 0.9)
  gammas = seq(0, 10, 2)
}


{
  m = 1
  p = 10
  e = c(0.02)
  r = 0.5
  gamma = c(5)
}


##### methods
{
  lassof = function(x,y){covlasso(x,y,cor.method = "pearson",scale.method = "sd", ada = F)$betahat_opt}
  alassof = function(x,y){covlasso(x,y,cor.method = "pearson",scale.method = "sd", ada = T)$betahat_opt}
  sparseShootingSf = function(x,y){shootings::sparseshooting(x,y)$coef[-1]}
  sltsf = function(x,y){fit = robustHD::sparseLTS(x,y);fit$coefficients[-1]}
  lgrf = function(x,y){covlasso(x,y,cor.method = "gaussrank", ada = F)$betahat_opt}
  lrpf = function(x,y){covlasso(x, y, ada = F)$betahat_opt}
  algrf = function(x,y){covlasso(x,y,cor.method = "gaussrank", ada = T)$betahat_opt}
  alrpf = function(x,y){covlasso(x, y, ada = T)$betahat_opt}

  mtds = list(Lasso = lassof,
              aLasso = alassof,
              sLTS = sltsf,
              sShootingS = sparseShootingSf,
              LGR = lgrf,
              LRP = lrpf,
              ALGR = algrf,
              ALRP = alrpf
  )
}

#### measurements
{
  msef<-function(betahat, beta){mean(((betahat-beta))^2)}
  mse5f<-function(betahat, beta){mean(((betahat-beta)[1:5])^2)}
  tpf<-function(betahat, beta){sum((as.logical(betahat)==as.logical(beta))[1:5])}
  tnf<-function(betahat,beta){sum((as.logical(betahat)==as.logical(beta))[-(1:5)])}
}



{
  systemtime = system.time({
    result <- foreach(m = ms,
                      .packages = c("robustHD","robcovsel", "shootings"))%:%
      foreach(p = ps)%:%
      foreach(e = es)%:%
      foreach(r = rs)%:%
      foreach(gamma = gammas)%dopar% {

        seed = m
        set.seed(seed = seed)

        beta = c(1,2,1,2,1,rep(0,p-5))
        dataset = robcovsel::genevar(n = 100, p = p, e = e, r = r, beta = beta, gamma = gamma)
        x = dataset$x
        y = dataset$y

        rst = list(lasso = NULL, slts = NULL, sss = NULL, algr = NULL, alrp = NULL)
        for (mtd in 1:length(mtds)) {
          try({
            betahat = mtds[[mtd]](x,y)
            mse = msef(betahat, beta)
            mse5 = mse5f(betahat, beta)
            tp = tpf(betahat, beta)
            tn = tnf(betahat, beta)
            rst[[mtd]] = c(m = m, n = 100, p = p, e = e, r = r, gamma = gamma, mtd = names(mtds)[mtd], seed =seed,
                           mse = mse, mse5 = mse5, tp = tp, fp = (p-5)-tn)
          }, TRUE
          )
          if(is.null(rst[[mtd]])){rst[[mtd]]=rep(NA, 11)}
        }
        rst
      }
  })



  save(result, file = "result_simu_p200.RData")
}
