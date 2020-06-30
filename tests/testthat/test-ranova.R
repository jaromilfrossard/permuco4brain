context("ranova")
library(permuco4brain)
library(permuco)
library(fields)
library(dplyr)
library(tidyr)
library(purrr)
library(Matrix)
library(MASS)


set.seed(42)

### create fake data

df_channel <- data.frame( channel = c("A","B","C","D"), x = c(0,1,1,1.2),y= c(1,0,1.2,1),z =c(0,0,0,0),stringsAsFactors = F)

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
  expand_grid(wt=c("w1","w2","w2"), id =1:12)%>%
  arrange(id)%>%
  mutate(bw = if_else(id>5,"b1","b2"),
         id = paste0("p",id))

#add error id
signal<-
  design%>%
  nest(data=-id)%>%
  mutate(signal_id =
           map(id,function(idi){0.5*MASS::mvrnorm(n=1, mu =matrix(rep(0,(df_corr%>%pull(corr_mat))[[1]]%>%ncol()),ncol=1),
              Sigma = (df_corr%>%pull(corr_mat))[[1]])}))%>%
  unnest(data)

#add error res

signal<-
  signal%>%
  mutate(signal =
           map(id,function(idi){0.5*MASS::mvrnorm(n=1, mu =matrix(rep(0,(df_corr%>%pull(corr_mat))[[1]]%>%ncol()),ncol=1),
                                                  Sigma = (df_corr%>%pull(corr_mat))[[1]])}))%>%
  transmute(signal = map2(signal,signal_id,function(si1,si2){si1+si2}))

#create matrix

signal<-
  signal%>%
  pull(signal)%>%
  do.call(what="rbind")%>%
  array(dim = c(nrow(design),ncol(.)/nrow(df_channel),nrow(df_channel)))

dimnames(signal)[[3]]<-c("A","B","C","D")


gi<- position_to_graph(df_channel,delta=1.2,name="channel")


rd_cm <- brainperm(signal~bw*wt + Error(id/wt),data=design,graph = gi,np =2,
                     method = "Rd_kheradPajouh_renaud",ncores=1,multcomp = c("clustermass"))

rde_t <- brainperm(signal~bw*wt + Error(id/wt),data=design,graph = gi,np =2,
                   method = "Rde_kheradPajouh_renaud",ncores=1,multcomp = c("troendle"))



test_that("Equal statistic rd and rde", {
  max_diff <- max(abs(unlist(lapply(rd_cm $multiple_comparison,function(effi)effi$uncorrected$statistic))-
                        unlist(lapply(rde_t$multiple_comparison,function(effi)effi$uncorrected$statistic))))

  expect_true(max_diff<1e-12)
})


