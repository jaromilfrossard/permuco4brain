context("Test permutation method")
library(permuco4brain)
library(permuco)
library(fields)
library(dplyr)
library(tidyr)
library(purrr)
library(Matrix)
library(MASS)


### create fake data

df_channel <- data.frame( channel = c("A","B","C","D"), x = c(0,1,1,1.2),y= c(1,0,1,1.2),z =c(0,0,0,0),stringsAsFactors = F)

df_corr<- expand.grid(channel = c("A","B","C","D"),sample = c(1:100),stringsAsFactors = F)%>%
  left_join(df_channel,by="channel")%>%
  arrange(channel)%>%
  nest(space_time = c(channel,sample,x,y,z))

#add covariance function
df_corr<-df_corr%>%
  mutate(cov_mat = purrr::map(space_time,function(dfi){
    mat_xy <- dfi%>%
      dplyr::select(x,y,z)%>%
      dist()%>%
      as.matrix()
    mat_t <- dfi%>%
      transmute(sample = sample/max(sample))%>%
      dist()%>%
      as.matrix()

    Exponential(1/3*mat_xy^2+190*mat_t^2,range = 1)
  }))


#correalation matrix
df_corr<- df_corr%>%
  mutate(corr_mat = purrr::map(cov_mat,function(mati){
    (Matrix::nearPD(mati,corr=T)$mat)%>%as.matrix()
  }))



#add_design
design <-
  expand_grid(A=c("a1","a2","b2"),B = c("b1","b2"),n = 1:6)%>%
  dplyr::select(-n)%>%
  filter(1:n()%in%c(1:(n()-2)))%>%
  mutate(x1 = runif(n(),0,12)%>%round(),
         x1 = x1-mean(x1))

signal = MASS::mvrnorm(n=nrow(design), mu =matrix(rep(0,(df_corr%>%pull(corr_mat))[[1]]%>%ncol()),ncol=1),
                       Sigma = (df_corr%>%pull(corr_mat))[[1]])%>%
  array(dim = c(nrow(design),ncol(.)/nrow(df_channel),nrow(df_channel)))

dimnames(signal)[[3]]<-c("A","B","C","D")


gi<- position_to_graph(df_channel,delta=1.2,name="channel")

fl_f <- brainperm(signal~x1*A*B,data=design,graph = gi,np =2, method = "freedman_lane",ncores=1,multcomp = c("clustermass","troendle"))
m_f <- brainperm(signal~x1*A*B,data=design,graph = gi,np =2, method = "manly",ncores=1,multcomp = c("clustermass","troendle"))

fl_t <- brainperm(signal~x1*A*B,data=design,graph = gi,np =2, method = "freedman_lane",ncores=1,multcomp = c("clustermass","troendle"),
                  test ="t")
m_t <- brainperm(signal~x1*A*B,data=design,graph = gi,np =2, method = "manly",ncores=1,multcomp = c("clustermass","troendle"),
                 test = "t")


test_that("Equal statistic freedman_lane and kennedy", {



  max_diff <- max(abs(unlist(lapply(fl_f$multiple_comparison,function(effi)effi$uncorrected$statistic))-
    unlist(lapply(m_f$multiple_comparison,function(effi)effi$uncorrected$statistic))))

  expect_true(max_diff<1e-12)
})
