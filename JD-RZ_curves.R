library(ggplot2);theme_set(theme_bw())
library(tidyr)
library(dplyr)

sigma_negbinom <- function(S, delta, kappa){
  out <- S^kappa+(S^kappa-1)/(kappa*delta)
  return(out)
}

sigma_gamma <- function(S, delta, kappa){
  out <- S^kappa
  return(out)
}

sigma_poisson <- function(S, delta, kappa){
  out <- 1 + log(S)/delta
  return(out)
}

sigma_geom <- function(S, delta, kappa){
  out <- S + (S-1)/delta
  return(out)
}

dat <- (expand.grid(  S=seq(from=0.001, to=1, by=0.001)
                    , delta=c(0.1, 1, 10, 1000)
                    , kappa=c(0.01, 0.1, 0.5, 1, 10, 100)
                    )
       %>% as_tibble()
       %>% mutate(NegBinom = sigma_negbinom(S, delta, kappa))
       %>% mutate(Gamma = sigma_gamma(S, delta, kappa))
       %>% mutate(Poisson = sigma_poisson(S, delta, kappa))
       %>% mutate(Geom = sigma_geom(S, delta, kappa))
)

(dat 
  |> pivot_longer(  cols = 4:7
                  , names_to="Distrib"
                  , values_to = "Sigma"
                  )
  )-> dat_long

# brkvec <- c(10^(-3:-1), 0.5)

dat_d01  <- filter(dat_long,delta==0.1)
dat_d1   <- filter(dat_long,delta==1)
dat_d10  <- filter(dat_long,delta==10)
dat_d1000 <- filter(dat_long,delta==1000)

fig_d01 <- (
  ggplot(dat_d01,aes(S,Sigma,col=Distrib,group=kappa))
  + geom_point(size=0.5,alpha=0.5)
  + facet_wrap(~kappa,scale="free",labeller = label_both)
  + geom_hline(yintercept=1)
  # + scale_y_log10()
  + scale_x_reverse()
  + ylim(0,1)
  #+ scale_colour_viridis_c(trans="log10", breaks = brkvec)
  + ggtitle("Sigma-S curve with delta=0.1")
)
print(fig_d01)
ggsave("Sigma-S_Curve_delta01.png",plot=fig_d01, path = "./docs/pix", width=3200,height=1800,units="px")


fig_d1 <- (
  ggplot(dat_d1,aes(S,Sigma,col=Distrib,group=kappa))
  + geom_point(size=0.5,alpha=0.5)
  + facet_wrap(~kappa,scale="free",labeller = label_both)
  + geom_hline(yintercept=1)
  # + scale_y_log10()
  + scale_x_reverse()
  + ylim(0,1)
  #+ scale_colour_viridis_c(trans="log10", breaks = brkvec)
  + ggtitle("Sigma-S curve with delta=1")
)
print(fig_d1)
ggsave("Sigma-S_Curve_delta1.png",plot=fig_d1, path = "./docs/pix", width=3200,height=1800,units="px")

fig_d10 <- (
  ggplot(dat_d10,aes(S,Sigma,col=Distrib,group=kappa))
  + geom_point(size=0.5,alpha=0.5)
  + facet_wrap(~kappa,scale="free",labeller = label_both)
  + geom_hline(yintercept=1)
  # + scale_y_log10()
  + scale_x_reverse()
  + ylim(0,1)
  #+ scale_colour_viridis_c(trans="log10", breaks = brkvec)
  + ggtitle("Sigma-S curve with delta=10")
)
print(fig_d10)
ggsave("Sigma-S_Curve_delta10.png",plot=fig_d10, path = "./docs/pix", width=3200,height=1800,units="px")



fig_d1000 <- (
  ggplot(dat_d1000,aes(S,Sigma,col=Distrib,group=kappa))
  + geom_point(size=0.5,alpha=0.5)
  + facet_wrap(~kappa,scale="free",labeller = label_both)
  + geom_hline(yintercept=1)
  # + scale_y_log10()
  + scale_x_reverse()
  + ylim(0,1)
  #+ scale_colour_viridis_c(trans="log10", breaks = brkvec)
  + ggtitle("Sigma-S curve with delta=1000")
)
print(fig_d1000)
ggsave("Sigma-S_Curve_delta1000.png",plot=fig_d1000, path = "./docs/pix", width=3200,height=1800,units="px")
