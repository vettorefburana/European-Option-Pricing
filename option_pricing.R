

library(ggplot2)


### Compute the payoff of a short european call option ############
S_t <- seq(5, 20, 1) 
K = 15
premium = 1

intrinsic_value <- S_t - K - premium
payoff_short_call <- pmin(premium,-intrinsic_value)

payoff <- rowSums(cbind(payoff_short_call))

results <- data.frame(cbind(S_t,payoff))

ggplot(results, aes(x=S_t)) + 
  geom_line(aes(y=payoff, color = 'Payoff')) +
  scale_colour_manual("", 
                      breaks = c("Payoff"),
                      values = c("darkgreen")) + ylab("Payoff")+
  ggtitle("Call Payoff") 

## Compute the payoff of a european put option in a single-period binomial tree ##############
u_prob = function(r, sigma, delta_t = 1) {
  
  u = exp(sigma*sqrt(delta_t))
  d = exp(-sigma*sqrt(delta_t))
  
  return(data.frame( cbind(u,d) ) )
}

q_prob = function(u, d, r, sigma, delta_t = 1) {
  
  q = ( (1+r)^delta_t - d )/( u - d )
  
  return(q)
}

put_eur = function(S0, K, sigma, r, t0 = 0, T = 1) {
  
  delta_t = T - t0
  
  probs = u_prob(r = r, sigma = sigma)
  u = probs$u
  d = probs$d
  
  qu = q_prob(r = r, sigma = sigma, u = u, d = d)
  qd = 1 - qu
 
  phi_u = max(K - S0*u, 0)
  phi_d = max(K - S0*d, 0)
  
  B = (1+r)^(delta_t)
  
  x = 1/B*( ( u*phi_d - d*phi_u )/(u - d) )
  y = 1/S0*( (phi_u - phi_d)/(u - d) ) 
  
  price = x*B + y*S0
  
  results = list()
  results$price = price
  results$q = c(qu, qd)
  results$ptf_weights = cbind(x, y)
  
  return(results)
}

## Compute the payoff of a call option using put-call parity #####
call_eur = function(S0, K, sigma, r, delta_t = 1){
  
  output = put_eur(sigma = sigma, r = r, S0 = S0, K = K)
  call = - K*(1-r)^delta_t  + S0 + output$price
  
  return(call)
  
}

## Verify that the european call option satisfies Merton's constraints #####
merton = function(S0, K, sigma, r, delta_t = 1){
  # returns true if constraint is satisfied
  
  call = call_eur(sigma = sigma, r = r, S0 = S0, K = K) 
  
  payoff = max(S0 - K*(1-r)^delta_t, 0)
  
  return( (call >= payoff ) )
}

merton(sigma = 0.1, r = 0.1, S0 = 10, K = 10) 

### Examine the behaviour of the payoff functions ##########
put_eur(sigma = 0.2, r = 0.1, S0 = 10, K = 10)

call_eur(sigma = 0.2, r = 0.1, S0 = 10, K = 10)

# Evolution of the put price with respect to the underlying price 
stock_price = c(seq(from = 1, to = 15, by = 0.1))

option_price = sapply(stock_price, put_eur, sigma = 0.2, r = 0.1, K = 10)

data <- data.frame(stock_price, option = unlist( option_price["price", ] ) )

ggplot(data, aes(x=stock_price, y=option)) +
  geom_line()

# Evolution of the call price with respect to the underlying price 
option_price = sapply(stock_price, call_eur, sigma = 0.2, r = 0.1, K = 10)

data <- data.frame(stock_price, option_price )

ggplot(data, aes(x=stock_price, y=option_price)) +
  geom_line()

## Compute the payoff of a european put option in a multi-period binomial tree using backward induction ######
final_payoff = function(S0, K, u, d){
  
  Puu = max(K-S0*u^2,0)
  Pud = max(K-S0*u*d,0)
  Pdd = max(K-S0*d^2,0)
  
  return( list(Puu = Puu, Pud = Pud, Pdd = Pdd) )
  
}

binomial_BI = function(S0,K, u, d, r, T, t0 = 0, N = 2){
  
  deltaT = (T-t0)/N # deltaT è il tempo totale diviso il numero di periodi
  
  qu = q_prob(r = r, sigma = sigma, u = u, d = d)
  qd = 1 - qu
  
  P_final = final_payoff(S0 = S0, K = K, u = u, d = d)
  
  Puu = P_final$Puu
  Pud = P_final$Pud
  Pdd = P_final$Pdd
  
  discount = 1/(1+r)^(deltaT)
  
  Pu = discount*(Puu*qu+Pud*qd)
  Pd = discount*(Pud*qu+Pdd*qd)
  P0 = discount*(Pu*qu+Pd*qd)
  
  return(P0)  
  
}

# Check the behaviour of the payoff function
r_annual = 0.04
r_sixmonths = (1 + (r_annual / 2))^2 - 1
binomial_BI(S0 = 10, T = 1, K = 10, u = 1.25, d = 0.8, r = r_sixmonths)

## Compute the payoff of a put option using a direct formula ######
binomial_DIR = function(K, S0, u, d, r, T, t0 = 0, N = 2){
  
  deltaT = (T-t0)/N 
  
  qu = q_prob(r = r, sigma = sigma, u = u, d = d)
  qd = 1 - qu
  
  P_final = final_payoff(S0 = S0, K = K, u = u, d = d)
  
  Puu = P_final$Puu
  Pud = P_final$Pud
  Pdd = P_final$Pdd
  
  P0 = 1/(1+r)^(T-t0)*(qu^2*Puu+2*qu*(1-qu)*Pud+(1-qu)^2*Pdd)
  
  return(P0)
}

binomial_DIR(S0 = 10, T = 1, K = 10, u = 1.25, d = 0.8, r = r_sixmonths)

# Evolution of the option price with respect to the underlying price
stock_price = c(seq(from = 0, to = 15, by = 0.1))
put_price= sapply(stock_price, binomial_DIR, T = 1, K = 10, u = 1.25, d = 0.8, r = 0.04)

data <- data.frame(stock_price,put_price)

ggplot(data, aes(x=stock_price, y=put_price)) +
  geom_line()

