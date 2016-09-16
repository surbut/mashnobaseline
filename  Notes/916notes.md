##Notes from 9.16 Meeting

* Instead of comparing ***mashnobaseline*** with ash and bmalite, let's try comparing ***mashnobaseline*** with old mash, where we simply center $\hat{C}$ and then model $\hat{C} - \mu$ as before we modeled $\hat{B}$ in exchangeable effects , with 
$$\hat{C} - \mu = \beta + E$$ where $$\beta \sim \sum_{p} \pi_{p} U_{k}$$ and $E \sim N(0,diag(S^2)$.
* First, I will apply mash to my existing simulated data, and then I will create a simpler set of simulations in which effects are either up(down) in a few subgroups to see if ***mashnobaseline*** is better at capturing the ***noncompetitive situation***, i.e., ensuring that just because average expression is higher in a few subgroups, it is not necessarily negative in the others.
 

 
