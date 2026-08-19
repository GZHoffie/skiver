# FASTA/FASTQ processing performance

Skiver uses `seq_io_parallel` for uncompressed FASTA/FASTQ input. Each worker
builds a thread-local `(k,v)`-mer set, and the sets are merged after parsing
completes. Gzip-compressed FASTA/FASTQ input is decompressed with
`rapidgzip-core` and parsed with Needletail when more than one thread is
available. With one thread, gzip input uses Needletail's sequential decoder.

## Benchmark setup

The following end-to-end `skiver sketch` benchmarks were run on
a machine with 16 logical CPUs and 62 GiB RAM. Release builds were used, the
input was warmed in the operating-system page cache, and every run used
`-c 10000`. Timings include parsing, `(k,v)`-mer extraction, merging, and
writing the sketch. Each table reports one full-file run.

The baseline was the release build from master commit `dd57abb` (v0.3.1).

### Uncompressed FASTQ

Input: `SRR13128014.fastq` (34 GiB).

|Implementation|Threads|Elapsed (s)|Speedup vs master|Time reduction|Peak RSS (MiB)|
|---|---:|---:|---:|---:|---:|
|Master / Needletail|1|106.62|1.00x|0.0%|294.8|
|`seq_io_parallel`|1|99.08|1.08x|7.1%|294.7|
|`seq_io_parallel`|2|50.99|2.09x|52.2%|394.9|
|`seq_io_parallel`|4|26.67|4.00x|75.0%|459.7|
|`seq_io_parallel`|8|15.24|7.00x|85.7%|499.7|
|`seq_io_parallel`|16|14.52|7.34x|86.4%|519.9|

The master and 16-thread sketches contained the same 3,521,210 observations
across 47,731 keys. The complete deserialized structures were equal after
canonicalizing observation order.

### Gzip-compressed FASTQ

Input: `ERR3152366.fastq.gz` (16 GiB compressed).

|Implementation|Threads|Elapsed (s)|Speedup vs master|Time reduction|Peak RSS (MiB)|
|---|---:|---:|---:|---:|---:|
|Master / Needletail|1|278.86|1.00x|0.0%|709.7|
|Needletail sequential fallback|1|280.96|0.99x|-0.8%|710.1|
|rapidgzip + Needletail|2|98.98|2.82x|64.5%|721.8|
|rapidgzip + Needletail|4|98.71|2.83x|64.6%|804.8|
|rapidgzip + Needletail|8|55.95|4.98x|79.9%|982.9|
|rapidgzip + Needletail|16|43.49|6.41x|84.4%|1282.2|

For gzip input, approximately three quarters of the requested thread budget is
assigned to decompression and the remainder to record processing. Consequently,
the 2- and 4-thread configurations both have one processing worker and perform
similarly on this workload.

The master and 16-thread sketches contained the same 3,331,044 observations
across 650,881 keys. The complete deserialized structures were equal after
canonicalizing observation order.

These numbers are specific to the machine, inputs, and subsampling rate above.
Performance and peak memory scale with read length, compression ratio, selected
thread count, and the number of retained `(k,v)`-mers.
