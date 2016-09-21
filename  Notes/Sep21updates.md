##Implementation

* Added simple sparse factor model where there are $2^R$ possible combinations from which Uk is simulated
* We want to see if MashNobaseline is better at at capturing the 1 1 0 0 0 case than Mash - for example, producing estimates that are closer to 1 1 0 0 0 than 0.2 0.2 -0.2 -0.2 -0.2 etc.
* How do we evaluate that? 
	* First we run  mashnobaseline with input summary statistics as *w = (LChat)* and the list of arrays *LVL*, then mash with *LChat* and *Se*, and we can compare RMSE.
	* I want to run on 8, 20 and 40 tissue cases to see if mashnobaseline improves as number of subgroups goes up.
	