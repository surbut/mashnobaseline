In **JOL**, we compute just one covariance matrix using ED with MNB and with MASH only. Critically, here we **cheat** and use the **true** associations for creating the covariance matrices and estimating the weights, 
$\hat{pi}$. We also examine the Likelihoods using just the truth matrix.

We find that the likelihood with MNB is -9553.341 and with MASH is -11,512 K. 

The likelihood with MNB and the truth is -9556.
The likelihood with MASH and the truth is -11925.

however, the RMSE with MNB is not really better than the RMSE with MASH (0.855) but it is better than using the truth:

Rmse.mash
## [1] 0.4698918
rmse.mnb
## [1] 0.4698918
rmse.truth
## [1] 0.3588735
rmse.mashtruth
## [1] 0.3756168




In *MASH Median*, we show that applying Mash to the median centered estimates gets an RMSE that is close to using *mnb* on the truth (0.4094) vs 0.3595518.





