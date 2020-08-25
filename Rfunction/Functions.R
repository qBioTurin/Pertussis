b_rate<-function(file, x = NULL)
{
    load(file)
    b_rates <- birth_rates
    return(unlist(b_rates,use.names = FALSE))
}

c_rate<-function(file, x = NULL)
{
    c_rates <- readRDS(file)
    return(matrix(unlist(c_rates), nrow = 3))
}

d_rate<-function(file, x = NULL)
{
    load(file)
    d_rates <- death_rates
    return(matrix(unlist(d_rates), nrow = 3))
}

v_rate<-function(file, x = NULL)
{
    load(file)
    ExpMin <- function( lambda, mu, xi, phases, probability ) {
        ( lambda / ( lambda + mu + xi ) ) ^ phases - probability
    }
    
    v_rates<-sapply(1:length(t(vaccination_coverage)), function(x){
        uniroot( f =ExpMin , interval = c(0,100), mu = death_rates[1,x], xi = 0.0030303, probability = t(vaccination_coverage)[x] , phases = 3 ,tol = 1e-8 )$root
    } )
    return(v_rates)
}

v_rate_wif<-function(file, x = NULL)
{
    load(file)
    ExpMin <- function( lambda, mu, xi, phases, probability ) {
        ( lambda / ( lambda + mu + xi ) ) ^ phases - probability
    }
    vaccination_coverage[c(33:43)] <- 0.8
    v_rates<-sapply(1:length(t(vaccination_coverage)), function(x){
        uniroot( f =ExpMin , interval = c(0,100), mu = death_rates[1,x], xi = 0.0030303, probability = t(vaccination_coverage)[x] , phases = 3 ,tol = 1e-8 )$root
    } )
    return(v_rates)
}

probability <- function(file, x = NULL)
{
    load(file)
    if( is.null(x) ){
        x <- c(runif(n=1,min=0,max=0.1), # 0.01
               runif(n=1,min=0,max=0.1), # 0.005
               runif(n=1,min=0,max=0.1), # 0.01
               0)
    }
    else{
        x <- c(x[c(1:3)],0.3)
        # x <- c(0.0008,0.00218,0.002985)
    }
    return(matrix(x, ncol = 1))
}

probability_wif <- function(file, x = NULL)
{
    load(file)
    if( is.null(x) ){
        x <- runif(n = length(probabilities), min=0, max=0.25)
        x[length(x)] = 0
    }
    else{
        x <- c(x[c(1:3)],0.2)
    }
    return(matrix(x, ncol = 1))
}

initial_marking <- function(file, x = NULL)
{
    load(file)
    if( is.null(x))
    {
        n_variables <- 3
        x <- runif(n_variables, min=0, max=1)
        n_a1<-sum(yini[c("S_a1", "R_a1_nv_l4")])
        x_a1<-c(runif(1,
                      min=x[1]*0.99*n_a1,
                      max=x[1]*n_a1),
                runif(1,
                      min=(1-x[1])*0.99*n_a1,
                      max=(1-x[1])*n_a1)
        )
        yini[c("S_a1", "R_a1_nv_l4")] <- x_a1
        n_a2 <- sum(yini[c("S_a2", "R_a2_nv_l1", "R_a2_nv_l2", "R_a2_nv_l3", "R_a2_nv_l4")])
        x_a2 <- c(runif(1,
                        min=x[2]*0.99*n_a2,
                        max=x[2]*n_a2),
                  runif(4,
                        min=((1-x[2])/4)*0.99*n_a2,
                        max=((1-x[2])/4)*n_a2)
        )
        yini[c("S_a2", "R_a2_nv_l1", "R_a2_nv_l2", "R_a2_nv_l3", "R_a2_nv_l4")] <- x_a2
        n_a3 <- sum(yini[c("S_a3", "R_a3_nv_l1", "R_a3_nv_l2", "R_a3_nv_l3", "R_a3_nv_l4")])
        x_a3 <- c(runif(1,
                        min=x[3]*0.99*n_a3,
                        max=x[3]*n_a3),
                  runif(4,
                        min=((1-x[3])/4)*0.99*n_a3,
                        max=((1-x[3])/4)*n_a3)
        )
        yini[c("S_a3", "R_a3_nv_l1", "R_a3_nv_l2", "R_a3_nv_l3", "R_a3_nv_l4")]<- x_a3
        return(matrix(unlist(yini), ncol = 1))
        
    }else{
        n_variables <- 12
        x <- x[c((length(x)-n_variables+1):length(x))]
        yini[c("S_a1", "R_a1_nv_l4")] <- sum(yini[c("S_a1", "R_a1_nv_l4")]) * (x[c(1:2)]/sum(x[c(1:2)]))
        yini[c("S_a2", "R_a2_nv_l1", "R_a2_nv_l2", "R_a2_nv_l3", "R_a2_nv_l4")] <- sum(yini[c("S_a2", "R_a2_nv_l1", "R_a2_nv_l2", "R_a2_nv_l3", "R_a2_nv_l4")]) * (x[c(3:7)]/sum(x[c(3:7)]))
        yini[c("S_a3", "R_a3_nv_l1", "R_a3_nv_l2", "R_a3_nv_l3", "R_a3_nv_l4")]<-  sum(yini[c("S_a3", "R_a3_nv_l1", "R_a3_nv_l2", "R_a3_nv_l3", "R_a3_nv_l4")]) * (x[c(8:12)]/sum(x[c(8:12)]))
        return(matrix(unlist(yini), ncol = 1))
    }
}