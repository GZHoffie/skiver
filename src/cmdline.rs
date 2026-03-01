use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about = "Skiver: Alignment-free estimation of sequencing error rates and spectra using (k,v)-mer sketches", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// Sketch the given sequencing files into kv-mer sketches.
    #[clap(display_order = 1)]
    Sketch(SketchArgs),

    /// Analyze a given sequencing file.
    #[clap(display_order = 2)]
    Analyze(AnalyzeArgs),

    /// For testing only: Try mapping the reads to reference genomes, and check how many k-mers are error-free.
    #[clap(display_order = 3)]
    Map(MapArgs),
}

#[derive(Args, Default)]
pub struct SketchArgs {
    #[clap(multiple=true, help_heading = "INPUT", help = "fasta/fastq files; gzip optional.")]
    pub files: Vec<String>,

    #[clap(short, default_value_t = 21, help_heading = "ALGORITHM", help ="Length of keys.")]
    pub k: u8,

    #[clap(short, default_value_t = 13, help_heading = "ALGORITHM", help ="Length of values.")]
    pub v: u8,

    #[clap(short, default_value_t = 1000, help_heading = "ALGORITHM", help = "Subsampling rate.")]
    pub c: usize,

    #[clap(short = 'f', default_value_t = 0, help_heading = "INPUT", help = "Number of bases to trim from the start of each read.")]
    pub trim_front: usize,

    #[clap(short = 'b', default_value_t = 0, help_heading = "INPUT", help = "Number of bases to trim from the end of each read.")]
    pub trim_back: usize,

    //#[clap(short, default_value_t = 4, help_heading = "ALGORITHM", help = "Number of threads.")]
    //pub threads: usize,

    #[clap(short, default_value_t = String::new(), help_heading = "OUTPUT", help = "Output file.")]
    pub output_path: String,

    #[clap(long, help_heading = "ALGORITHM", help = "Use the forward strand of the reads only. Default: use both forward and reverse strands of the reads.")]
    pub forward_only: bool,
}


#[derive(Args, Default, Clone)]
pub struct AnalyzeArgs {
    #[clap(multiple=true, help_heading = "INPUT", help = "fasta/fastq files; gzip optional.")]
    pub files: Vec<String>,

    #[clap(short = 'k', default_value_t = 21, help_heading = "ALGORITHM", help ="Length of keys.")]
    pub k: u8,

    #[clap(short = 'v', default_value_t = 13, help_heading = "ALGORITHM", help ="Length of values.")]
    pub v: u8,

    #[clap(short = 'c', default_value_t = 1000, help_heading = "ALGORITHM", help = "Subsampling rate.")]
    pub c: usize,

    #[clap(short = 'l', long = "lower-bound", default_value_t = 5, help_heading = "ALGORITHM", help = "Lower bound for the number of times the consensus appears in the read for it to be considered in the profiling.")]
    pub lower_bound: u32,

    #[clap(long, help_heading = "ALGORITHM", help = "Use the forward strand of the reads only. Default: use both forward and reverse strands of the reads.")]
    pub forward_only: bool,

    #[clap(long = "use-all", help_heading = "ALGORITHM", help = "Not excluding the outliers.")]
    pub use_all: bool,

    #[clap(short = 'e', long = "outlier-threshold", default_value_t = 3.0, help_heading = "ALGORITHM", help = "The multiplier to the IQR for defining outliers.")]
    pub outlier_threshold: f32,

    #[clap(long = "num-experiments", default_value_t = 100, help_heading = "ALGORITHM", help = "Number of experiments in bootstrapping for estimating the parameters.")]
    pub num_experiments: u32,

    #[clap(long = "bootstrap-sample-rate", default_value_t = 0.1, help_heading = "ALGORITHM", help = "Proportion of data points to sample per experiment in bootstrapping.")]
    pub bootstrap_sample_rate: f32,

    #[clap(short = 'r', long = "reference", help_heading = "ALGORITHM", help = "Reference genomes.")]
    pub reference: Option<String>,

    #[clap(short = 'f', long = "trim-front", default_value_t = 0, help_heading = "INPUT", help = "Number of bases to trim from the start of each read.")]
    pub trim_front: usize,

    #[clap(short = 'b', long = "trim-back", default_value_t = 0, help_heading = "INPUT", help = "Number of bases to trim from the end of each read.")]
    pub trim_back: usize,

    #[clap(long, default_value_t = 2, help_heading = "ALGORITHM", help = "Number of estimated hazard ratios to ignore from the largest v.")]
    pub ignore_last_hazard_ratios: usize,

    //#[clap(short = 't', long = "threads", default_value_t = 4, help_heading = "ALGORITHM", help = "Number of threads.")]
    //pub threads: usize,

    #[clap(short = 'o', long = "verbose-output", help_heading = "OUTPUT", help = "Output file.")]
    pub output_path: Option<String>,

    #[clap(long, default_value_t = String::from("sum_ratio"), hidden = true, help = "One of 'slope', 'linear_fit', 'ratio_mean', 'sum_ratio'.")]
    pub estimation_method: String,

    #[clap(long, default_value_t = String::from("weibull"), help = "Model used to fit the hazard rates vs. t. Should be one of 'constant' (assuming that the hazard rate is constant over t), 'weibull' (assuming T follows a discrete Weibull distribution).")]
    pub hazard_model: String,

    #[clap(long, help_heading = "OUTPUT", help = "Output the estimated hazard ratio and their confidence intervals as a csv file.")]
    pub hazard_rate: Option<String>,
}

#[derive(Args, Default)]
pub struct MapArgs {
    #[clap(multiple=true, help_heading = "INPUT", help = "fasta/fastq files; gzip optional.")]
    pub files: Vec<String>,

    #[clap(short, default_value_t = 21, help_heading = "ALGORITHM", help ="Length of keys.")]
    pub k: u8,

    #[clap(short, default_value_t = 1000, help_heading = "ALGORITHM", help = "Subsampling rate.")]
    pub c: usize,

    #[clap(short, default_value_t = 100, help_heading = "ALGORITHM", help = "Read sampling rate.")]
    pub sample_rate: usize,

    #[clap(short, default_value_t = 5, help_heading = "ALGORITHM", help = "Lower bound for the number of times the consensus appears in the read for it to be considered in the profiling.")]
    pub lower_bound: u32,

    #[clap(long, help_heading = "ALGORITHM", help = "Use the forward strand of the reads only. Default: use both forward and reverse strands of the reads.")]
    pub forward_only: bool,

    #[clap(short, help_heading = "ALGORITHM", help = "Reference genomes.")]
    pub reference: String,

    #[clap(short = 'f', default_value_t = 0, help_heading = "INPUT", help = "Number of bases to trim from the start of each read.")]
    pub trim_front: usize,

    #[clap(short = 'b', default_value_t = 0, help_heading = "INPUT", help = "Number of bases to trim from the end of each read.")]
    pub trim_back: usize,

    //#[clap(short, default_value_t = 4, help_heading = "ALGORITHM", help = "Number of threads.")]
    //pub threads: usize,

    #[clap(short, help_heading = "OUTPUT", help = "Verbose output per-read k-mer hit information to stdout.")]
    pub print_verbose: bool,
}