## required library
library(deSolve) ## ode
library(pracma)  ## lsqnonlin
library(data.table)
library(magrittr)
library(parallel)
#setwd("~/ShinyApps/covid-seoul")


## Data 
info.sl <- fread("datasl_utf8.csv")[, `확진일` := lubridate::ymd(`확진일`)][`확진일` >= ymd("2020-11-01")][]

info.n <- fread("report.txt")[-c(1, 2), c(2, 4)][, `인구` := as.integer(gsub(",", "", `인구`))][]
for (v in info.n[["자치구"]][-1]){
  info.sl[["지역"]][grep(v, info.sl[["지역"]])] <- v
}


inci <- rbind(info.sl[, .(`지역` = "합계", .N), by = "확진일"],
              info.sl[`지역` %in% info.n[["자치구"]][-1], .N, by = c("확진일", "지역")]
              )[order(`지역`, `확진일`)]


#sigu <- "강남구"

## Run: All + sigu
res <- mclapply(info.n[["자치구"]], function(sigu){
  
  ## load data
  date <- inci[`지역` == sigu, `확진일`]
  incidata <- inci[`지역` == sigu, N]
  incidata.g <- cumsum(incidata)
  
  
  ## make discrete time array, which has same length with the loaded data
  times.g   <- seq(0, length(incidata.g) - 1, by = 1)
  
  
  tt <- seq(min(times.g), max(times.g), len = max(times.g) +1)
  B  <- bs(tt, knots = seq(min(times.g), max(times.g), len = ceil(length(date)/7)), Boundary.knots = c(min(times.g), max(times.g)+1)) # boundaryknots: ode goes beyond t0, t1
  
  ## main ODE form
  covODE <- function(time, state, modelparams) {
    with(as.list(c(state, modelparams)), {
      beta <- exp(predict(B, time) %*% modelparams[-c(1, 2)])
      N <- S + E + I + R
      dS <- (-1 * beta)*S*I/N 
      dE <- beta*S*I/N - theta*E
      dI <- theta*E - gamma*I
      dR <- gamma*I
      dCummrep <- gamma*I
      #이건 모수추정 과정에서 누적데이터에 상응하는 모델의 적분형태를 위한 부분 입니다.
      return(list(c(dS, dE, dI, dR, dCummrep)))
    })
  }
  
  ## set known model parameters !!(these can be adjusted by user)!!
  theta.g <- 1/3.5
  # 1/theta 값은 증상발현까지 잠복기간 입니다.
  gamma.g <- 1/6.8
  # 1/alpha 값은 증상발생환자의 입원까지의 평균 기간 입니다.
  
  ## setting fixed parameters
  fixparams.g  <- c(theta = theta.g, gamma = gamma.g)
  
  ## setting model initial conditions !!(these can be adjusted by user)!!
  init.g    <- c(S = info.n[`자치구` == sigu, `인구`] - incidata[1] - incidata[1] - incidata[1], E = incidata[1], I = incidata[1], R = incidata[1], Cummrep = incidata[1])
  # 모델 초기조건 입니다. 현재 이 셋팅은 대구 데이터에 맞춘 셋팅 입니다.
  
  ## parameter (beta) estimation
  ## Poisson likelihood
  ftemp2 <- function(x) {
    modelparams <- c(fixparams.g, x)
    tempoutput <- as.data.frame(ode(y = init.g, times = times.g, func = covODE, parms = modelparams))
    #fit  = ode.model[, "i"]
    llik <- -sum(dpois(incidata[-1], lambda = diff(tempoutput$Cummrep), log = TRUE)) 
    
    return(llik)
  }
  
  
  ## Least square
  ftemp <- function(x) {
    modelparams = c(fixparams.g, x)
    tempoutput = as.data.frame(ode(y = init.g, times = times.g, func = covODE, parms = modelparams))
    return(abs(tempoutput$Cummrep-incidata.g))
  }
  
  #inibeta = rep(0.5, ceil(length(times.g)/7))    ## initial assumption, !!(these can be adjusted by user)!!
  #initbeta   = sort(runif(ncol(B), -1, 1))
  initbeta <- rep(0, ncol(B))
  
  fitresult <- lsqnonlin(ftemp, initbeta) ## least square fit
  #fitresult <- optim(initbeta, ftemp2, method = "L-BFGS-B")  ## poisson fit
  betas <- fitresult[[1]]   
  
  #beta = betas[1]
  #beta_F = betas[2]
  # beta는 전파율이며 beta_F는 행동변화를 야기하는 행동변화 전파율 입니다.
  
  ## solve the model
  #params = c(fixparams.g, beta = beta, beta_F = beta_F)
  params = c(fixparams.g, betas)
  output <- ode(y = init.g, times = times.g, func = covODE, parms = params)
  output <- data.table(output)
  output[, `:=`(time = date, obs.cum = incidata.g, beta = exp(B %*% betas))]
  names(output)[6] <- "exp.cum"
  return(output)
}, mc.cores = 13)


saveRDS(res, "res.RDS")


library(showtext);library(ggplot2);library(gridExtra)
showtext_opts(dpi = 90)                                                ## same to rmd chunk
font_add("NanumGothic", "/usr/share/fonts/truetype/nanum/NanumGothic.ttf")
showtext_auto()

plot.cum <- lapply(1:length(info.n[["자치구"]]), function(v){
  ggplot(melt(res[[v]][, .(time, exp.cum, obs.cum)], id = "time", variable.name = "Type"), aes(time, value, color = Type)) + geom_line() + xlab("") + 
    scale_x_date(date_labels = "%b %d", date_breaks = "1 week") + theme_bw() + ggtitle(info.n[["자치구"]][v])
})

plot.beta <- lapply(1:length(info.n[["자치구"]]), function(v){
  ggplot(res[[v]], aes(time, beta)) + geom_line() + ylim(c(0, 1)) + scale_x_date(date_labels = "%b %d", date_breaks = "1 week") + 
    theme_bw() + xlab("") + ylab(expression(beta)) + ggtitle(info.n[["자치구"]][v])
})

 
allplot.cum <- marrangeGrob(plot.cum, nrow = 13, ncol = 2)
allplot.beta <- marrangeGrob(plot.beta, nrow = 13, ncol = 2)

ggsave("beta.jpg", allplot.beta, width = 5, height = 20)


## plot
ggplot(melt(output[, .(time, exp.cum, obs.cum)], id = "time", variable.name = "Type"), aes(time, value, color = Type)) + geom_line() + xlab("") + scale_x_date(date_labels = "%b %d", date_breaks = "1 week") + theme_bw()

#plot(date, incidata.g)        ## data
#lines(date, output$Cummrep)   ## model

ggplot(output, aes(time, beta)) + geom_line() + ylim(c(0, 1)) + scale_x_date(date_labels = "%b %d", date_breaks = "1 week") + theme_bw() + xlab("") + ylab(expression(beta))

#plot(date, exp(B %*% betas), col = 2, lty = 2, ylim = c(0, 1))
