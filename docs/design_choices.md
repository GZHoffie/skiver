# Design choices of skiver

In this file, I document some key design choices of skiver as well as some known issues, which can probably be useful for future reference and can be possible directions for future improvements.

## Definition of sequencing error rate

Usually, the sequencing error rate is defined using alignment results (see [Heng Li's blog](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity) on this topic),  which can vary a lot with different alignment scoring scheme. In light of this,
I define the hazard rate to be the probability that **the next base disagrees with the consensus (reference genome)**. The new definition is not dependent on ways of alignment, and thus is more consistent.

Unfortunately, the new definition can also be ambiguous since the base disagreeing with the consensus doesn't necessarily means that the sequencing error happens at this position. On one hand, the consensus can be wrong if the coverage is not high enough. On the other hand, let's say an insertion happens in squence `ACGT`, giving the read `ACGGT`, it's hard to say whether which `G` is the insertion error. According to our definition, we will treat the last `G` as the insertion error. The same thing happens for a deletion in a homopolymer region.

This ambiguity seems to have little effect when estimating Weibull parameters for hazard/survival rate estimation, but might cause a problem in quality score calibration, in which case, we may wrongly accuse a base with certain Q-score to be erroneous, causing error.

## Why survival analysis?

The first version of skiver assumes that the error is distributed uniformly at random, so the number of errors inside the value of each (*k*, *v*)-mer should follow a Poisson distribution. However, in both simulated and real data, this assumption largely underestimates the error rate.

One key observation is that the probability of an error happening in the value is not unconditioned - it is conditioned on the key agree with the consensus (i.e. being free of sequencing errors). It is thus natural to borrow concepts from survival analysis, where the definition of hazard rate is similiar to the conditional probability described above. 


## Finding parameters for the discrete Weibull

We assume that $T$, the time until first occurrence of sequencing error, follows a discrete Weibull distribution. To find its parameters, we apply complementary log-log transform to the estimated hazard rate and perform a linear regression.

One problem with this is that **it doesn't work with simulated perfect reads**, as the clog-log transformation of a zero hazard rate goes to -infinity, making the regression impossible. As of now, I set a lower limit to clip the hazard rate ($10^{-4}$) to keep the regression stable even for perfect reads. As a result, for simulated perfect reads, skiver will report an error rate of $10^{-4}$.

Another problem is that we can only observe hazard rate $h(t)$ for large $t$, $k+1\leq t\leq k+v$. As a result, to estimate $h(0)$, the per-base error rate, we need to extrapolate the observed $h(t)$. This is usually unstable, and the resulting $h(0)$ can vary a lot with different parameter choices. Fortunately, the resulting hazard rate/survival rate estimation remains relatively robust in my current experiments.


## Outlier filter

Due to repetitions/presence of multiple alleles, a key can be associated with multiple values of high count. Increasing $k$ is an option, but it makes the regression less stable (see above).

The first version of outlier filter simply takes the non-zero hazard rate at each $t$ and filter out base on the IQR. It works on simple simulated/real data, but for more complicated metagenomic samples, the non-zero hazard rate has higher variance and thus a larger IQR, making the method not so effective.

The newer version of outlier filter takes the keys with zero hazard rate into consideration. It works iteratively - estimate hazard rates, exclude outliers based on a Binomial distribution, and estimate hazard rates again, until they converges. This is the best solution I can come up with so far, and it produces reasonable results even for complicated metagenomic samples.

However, this filter is not perfect. I observe that it can be too strict for keys with high coverage, and throwing them away has a large impact on the estimation of hazard rates. At the same time it can be too lenient for keys with low coverage, making the results less accurate.

