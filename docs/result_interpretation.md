# Interpreting skiver reports

After running `skiver analyze`, it produces files with the specified prefix,

|File name|Content|
|---------|-------|
|`[prefix].summary_error_rate.csv`|Estimated error rate and coverage|
|`[prefix].kvmer.csv`|Info of the sketched (*k*, *v*)-mers|
|`[prefix].hazard_rate.csv`|Estimated hazard rate at different *t*|
|`[prefix].summary_error_spectrum.csv`|Frequency of each error type|
|`[prefix].summary_error_spectrum_dependence_on_t.csv`|Dependence of the frequency of each error type on *t*|
|`[prefix].summary_phred.csv`|Empirical error rate of the Phred scores|
|`[prefix].summary_read_position.csv`|Empirical error rate of each position in the read|
|`[prefix].summary_gc_content.csv`|Empirical error rate of the reads with a certain GC-content|
|`[prefix].summary_read_length.csv`|Empirical error rate of the reads in read-length intervals|

Below, we explain the fields reported in each csv file.

### `[prefix].summary_error_rate.csv`

[Here](./example/SRR7498042.summary_error_rate.csv) is an example of the error rate report file. 

|Fields|Interpretation|
|------|--------------|
|`per_base_error_rate`|Estimated per-base error rate. In other words, the probability of a random base is erroneous. </br> This value is estimated by $\hat{\varepsilon}_{\text{perbase}}:=1-\exp(-\hat{\lambda})$.|
|`mean_hazard_rate`|Estimated mean hazard rate: $1$ minus the geometric mean of $(1-h(t))$ over $t$. This value can be used to estimate the survival rate (i.e. the chance of a k-mer being free of sequencing error) by raising this number to the power of $k$.|
|`lambda`|Estimated $\lambda$ of the discrete Weibull distribution of the survival model. A larger $\lambda$ means a higher per-base error rate.|
|`beta`|Estimated $\beta$ of the discrete Weibull distribution of the survival model. Smaller $\beta$ means that the errors are likely clustered together, and $\beta$ closer to 1 means that the errors are close to be randomly distributed.|
|`key_median_coverage`|Median coverage of the keys of the sketched $(k,v)$-mers that passes the outlier filter.|
|`true_median_coverage`|Estimated true coverage of the keys (if they are free of sequencing errors) of the sketched $(k,v)$-mers that passes the outlier filter. The true coverage is estimated by taking the observed median coverage above, and divide by $\hat{S}(k)=\exp(-\hat{\lambda} k^{-\hat{\beta}})$.|

Each value is paired with a `5-95th_percentile` confidence interval, which is estimated using bootstrapping experiments.

### `[prefix].kvmer.csv`

[Here](./example/SRR7498042.kvmer.csv) is an example of the $(k,v)$-mer summary file. 

|Fields|Interpretation|
|------|--------------|
|`key`|The key sequence in the $(k,v)$-mer.|
|`consensus_value`|The most frequently appearing value that follows the key.|
|`passes_filter`|A boolean value that is `true` if the key passes the outlier filter and is included in the analysis, and `false` otherwise.|
|`homopolymer_length`|The longest length of homopolymer that is recorded in the consensus value, can be used to study the dependence of error rate on homopolymer length.|
|`consensus_count`|Number of time the `consensus_value` appear in the reads that that follows the key.|
|`neighbor_count`|Number of time a value that has an edit distance of 1 to the `consensus_value` appear in the reads that that follows the key.|
|`total_count`|Number of time the `key` appear in the reads.| 
|`AC`, `AG`, ..., `T_`|Number of time a value that has an edit distance of 1 with a specific edit type (`AC` means substitution of `A` to `C`, `T_` means a deletion of base `T`, etc.) to the `consensus_value` appear in the reads.|
|`consensus_count_up_to_v${i}`|Number of values that agree with the consensus value up to the `i`-th base. The hazard rate is calculated by these columns.|



### `[prefix].hazard_rate.csv`

[Here](./example/SRR7498042.hazard_rate.csv) is an example of the estimated hazard rate summary.

|Fields|Interpretation|
|------|--------------|
|`num_candidates`|Total number of values that agree with the consensus value up to T=`t-1`|
|`num_survival`|Total number of values that agree with the consensus value up to T=`t`|
|`hazard_ratio`|Estimated hazard rate $\hat{h}(t)$, which is calculated by `num_survival` divided by `num_candidates`.|
|`5th_percentile`, `95th_percentile`|Confidence interval of the hazard rate estimated by bootstrapping experiments.|

### `[prefix].summary_error_spectrum.csv`

[Here](./example/SRR7498042.summary_error_spectrum.csv) is an example of the error spectrum.

|Fields|Interpretation|
|------|--------------|
|`operation`|Edit operation, including substitution (e.g. `A>C`), insertion (e.g. `->A`), and deletion (e.g. `A>-`).|
|`prev_base`|Base preceding the edit operation.|
|`next_base`|Base after the edit operation.|
|`total`|Number of times the error type happen in the read set, including both the forward strand of the read and the reverse complement of the read.|
|`forward`|Number of times the error type happen in the read set, including only the forward strand of the read.|

### `[prefix].summary_error_spectrum_dependence_on_t.csv`

[Here](./example/SRR7498042.summary_error_spectrum_dependence_on_t.csv) is an example of the error spectrum dependence on t.

|Fields|Interpretation|
|------|--------------|
|`operation`, `prev_base`, `next_base`, `total`|Same as `summary_error_spectrum.csv`.|
|`freq_at_t${i}`|Number of times the error type happen at the position t=`i` in the sketched $(k,v)$-mers.|

### `[prefix].summary_phred.csv`

[Here](./example/SRR7498042.summary_phred.csv) is an example of the Phred score summary file.

|Fields|Interpretation|
|------|--------------|
|`qscore`|The reported Phred score in the fastq file|
|`empirical_qscore`|The observed error rate of the base that has the Phred score.</br> **Note**: the error here includes substitution, insertion, and deletion. As long as the base differs from the consensus value, it is regarded as an error.|
|`per_base_error_rate`|The per base error rate that is estimated using the clog-log regression.|
|`per_base_error_rate_5-95th_percentile`|The confidence interval of per-base error rate estimated using bootstrapping.|
|`num_correct`|Number of bases in the values that has the reported phred score that are correct (agree with the consensus).|
|`num_error`|Number of bases in the values that has the reported phred score that are incorrect (disagree with the consensus).|

### `[prefix].summary_gc_content.csv`

[Here](./example/SRR7498042.summary_gc_content.csv) is an example of the GC-content summary file.

|Fields|Interpretation|
|------|--------------|
|`gc_content_min`, `gc_content_max_exclusive`|Represent a range of reads whose GC-content is >= `gc_content_min` and < `gc_content_max_exclusive`|

The other fields have the identical meaning as that of `summary_phred.csv`.

### `[prefix].summary_read_length.csv`

This file estimates the dependence of error rate on read length using read-length bins. The default bin width is 10kb and can be changed with `--read-length-step`.

|Fields|Interpretation|
|------|--------------|
|`read_length_min`, `read_length_max_exclusive`|Represent a range of reads whose trimmed read length is >= `read_length_min` and < `read_length_max_exclusive`, e.g. with the default setting: `0,10000`, `10000,20000`, and so on.|

The other fields have the identical meaning as that of `summary_phred.csv`.

### `[prefix].summary_read_position.csv`

[Here](./example/SRR7498042.read_position.csv) is an example of the error rate dependence on position in the read.

> [!NOTE]  
> The error rate is calculated by `num_error / (num_correct + num_error)`, instead of the clog-log regression used for Phred score and GC-content. This is because for long reads, it is hard to have enough bases in the sketch to obtain accurate estimate of hazard rates for each position in the read for each *t*. The error rate is usually lower than the per-base error rate, but the general trend would be similar.

|Fields|Interpretation|
|------|--------------|
|`index`|Position from the two ends of the read.|
|`from_start`|Whether we are counting `index` from the start of the read (true) or the end of the read (false).|
|`num_correct`|Number of bases that are at the position that are correct (agree with the consensus).|
|`num_error`|Number of bases that are at the position that are incorrect (disagree with the consensus).|
|`error_rate`|Percentage of bases that are at the position that are incorrect (disagree with the consensus).|
