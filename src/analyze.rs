use crate::kvmer::*;
use crate::utils::*;
use crate::inference::*;
use crate::cmdline::AnalyzeArgs;

use simple_logger::SimpleLogger;
use log::{info, warn, error};
use glob::glob;
use needletail::parse_fastx_file;
use std::fs::{self, OpenOptions};
use std::path::Path;

const OUTPUT_SUFFIXES: [&str; 15] = [
    "hazard_rate.csv",
    "kvmer.csv",
    "summary_error_rate.csv",
    "survival_rate.csv",
    "summary_error_spectrum.csv",
    "summary_error_spectrum_dependence_on_t.csv",
    "summary_read_position.csv",
    "summary_phred.csv",
    "summary_gc_content.csv",
    "plot_spectrum.pdf",
    "plot_coverage.pdf",
    "plot_hazard_survival.pdf",
    "plot_qscore_calibration.pdf",
    "plot_gc_content.pdf",
    "plot_read_position.pdf",
];

fn survival_rate_to_csv(lambda: f32, beta: f32) -> String {
    let mut result = String::from("t,survival_rate\n");
    for t in 1..=100 {
        let survival_rate = (-(lambda as f64) * (t as f64).powf(beta as f64)).exp();
        result.push_str(&format!("{},{:.6}\n", t, survival_rate));
    }
    result
}

fn validate_output_prefix(prefix: &str) -> Result<(Option<String>, Vec<String>), String> {
    if prefix.is_empty() {
        return Err("Output prefix cannot be empty.".to_string());
    }

    let prefix_path = Path::new(prefix);
    let file_name = prefix_path.file_name()
        .ok_or_else(|| format!("Output prefix '{}' does not contain a file name.", prefix))?;
    if file_name.is_empty() {
        return Err(format!("Output prefix '{}' does not contain a file name.", prefix));
    }

    let parent = prefix_path.parent().filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let created_dir = if !parent.exists() {
        fs::create_dir_all(parent)
            .map_err(|e| format!("Could not create output directory '{}': {}", parent.display(), e))?;
        Some(parent.display().to_string())
    } else {
        if !parent.is_dir() {
            return Err(format!("Output directory '{}' is not a directory.", parent.display()));
        }
        None
    };

    let probe = parent.join(format!(".skiver-write-test-{}", std::process::id()));
    OpenOptions::new().write(true).create_new(true).open(&probe)
        .map_err(|e| format!("Output prefix '{}' is not writable: {}", prefix, e))?;
    fs::remove_file(&probe)
        .map_err(|e| format!("Could not remove output validation file '{}': {}", probe.display(), e))?;

    let existing = OUTPUT_SUFFIXES.iter()
        .map(|suffix| format!("{}.{}", prefix, suffix))
        .filter(|path| Path::new(path).exists())
        .collect();
    Ok((created_dir, existing))
}

pub fn analyze(args: AnalyzeArgs) {
    SimpleLogger::new().with_level(log::LevelFilter::Info).init().unwrap();

    let (created_dir, existing_outputs) = match validate_output_prefix(&args.output_prefix) {
        Ok(result) => result,
        Err(message) => {
            error!("{}", message);
            std::process::exit(1);
        }
    };
    if let Some(directory) = created_dir {
        info!("Created output directory '{}'.", directory);
    }
    if !existing_outputs.is_empty() {
        warn!(
            "Output files with prefix '{}' already exist and will be overwritten.",
            args.output_prefix
        );
    }

    // [TODO] Multithreaded version is under development.
    //rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

    //info!("Using {} threads for analysis.", args.threads);

    let mut kvmer_set = KVmerSet::new(args.k, args.v, !args.forward_only);

    // Expand globs and try each input as FASTA/FASTQ, then as a sketch.
    let mut raw_files: Vec<String> = Vec::new();
    let mut sketch_file_count = 0;
    let mut has_invalid_files = false;
    for file in &args.files {
        let entries = match glob(file) {
            Ok(entries) => entries,
            Err(_) => {
                error!("{} is not a valid FASTA/FASTQ or kv-mer sketch file; skipping.", file);
                has_invalid_files = true;
                continue;
            }
        };
        let mut matched = false;
        for entry in entries {
            matched = true;
            match entry {
                Ok(path) => {
                    if let Some(file_str) = path.to_str() {
                        if parse_fastx_file(file_str).is_ok() {
                            raw_files.push(file_str.to_string());
                        } else {
                            match kvmer_set.load(file_str) {
                                Ok(()) => sketch_file_count += 1,
                                Err(message) => {
                                    error!("{}; skipping.", message);
                                    has_invalid_files = true;
                                }
                            }
                        }
                    } else {
                        error!("{} is not a valid FASTA/FASTQ or kv-mer sketch file; skipping.", path.display());
                        has_invalid_files = true;
                    }
                }
                Err(e) => {
                    error!("{} is not a valid FASTA/FASTQ or kv-mer sketch file; skipping.", e.path().display());
                    has_invalid_files = true;
                }
            }
        }
        if !matched {
            error!("{} is not a valid FASTA/FASTQ or kv-mer sketch file; skipping.", file);
            has_invalid_files = true;
        }
    }

    if (!raw_files.is_empty() && sketch_file_count > 0)
        || sketch_file_count > 1
    {
        error!("{}", "The current version of skiver analyze only supports either exactly one kv-mer sketch file or one or more FASTA/FASTQ files (gzip optional) as input.");
        std::process::exit(1);
    }
    if raw_files.is_empty() && sketch_file_count == 0 {
        if !has_invalid_files {
            error!("{}", "The current version of skiver analyze only supports either exactly one kv-mer sketch file or one or more FASTA/FASTQ files (gzip optional) as input.");
        }
        std::process::exit(1);
    }

    let c = args.c.unwrap_or_else(|| {
        let raw_refs: Vec<&str> = raw_files.iter().map(|s| s.as_str()).collect();
        let (auto_c, est_file_size) = estimate_c_from_raw_files(&raw_refs);
        info!("Total weighted input sequence file size: {:.2} GB", est_file_size as f64 / (1024.0 * 1024.0 * 1024.0));
        info!("Auto-determined subsampling rate: -c {}", auto_c);
        auto_c
    });

    // Read query files
    info!("Processing query files...");
    let threads = if args.threads == 0 { std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1) } else { args.threads };
    info!("Using {} FASTA/FASTQ worker thread(s).", threads);
    for file_str in &raw_files {
        kvmer_set.add_file_to_kvmer_set_with_threads(file_str, c, args.trim_front, args.trim_back, threads);
    }
    info!("Finished processing query files.");

    let analyzer = ErrorAnalyzer::new(args.clone());

    
    let stats: KVmerStats;
    if let Some(reference) = &args.reference {
        if args.lower_bound.is_none() {
            info!("Reference is provided. Using default lower bound of 0.");
        }
        let lower_bound = args.lower_bound.unwrap_or(0);

        let mut reference_kvmer_set = KVmerSet::new(args.k, args.v, true);
        reference_kvmer_set.add_file_to_kvmer_set_with_threads(reference, c, args.trim_front, args.trim_back, threads);
        info!("Loaded reference file: {}", reference);

        stats = kvmer_set.get_stats_with_reference(lower_bound, &reference_kvmer_set);
    } else {
        let lower_bound = args.lower_bound.unwrap_or(10);
        //println!("Error rate: {}", kvmer_set.get_stats(args.threshold));
        stats = kvmer_set.get_stats(lower_bound);
    }
    // if reference is set, the filter should be disabled
    // [FIXME] enable --use-all by default
    if args.reference.is_some() && !args.use_all {
        warn!("If reference is provided, --use-all is recommended.");
    }

    let spectrum = analyzer.analyze(&stats);
    let analysis_output = format!("{}\n{}", header_str(!args.forward_only), spectrum_to_str(&spectrum, !args.forward_only));

    fs::write(format!("{}.summary_error_rate.csv", args.output_prefix), &analysis_output).unwrap();
    fs::write(
        format!("{}.survival_rate.csv", args.output_prefix),
        survival_rate_to_csv(spectrum.estimated_lambda.0, spectrum.estimated_beta.0),
    ).unwrap();
    crate::plot::generate(&args.output_prefix, &crate::plot::PlotOptions::default());
    info!("Output written to prefix {}.", args.output_prefix);
}
