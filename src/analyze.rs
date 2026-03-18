use crate::kvmer::*;
use crate::utils::*;
use crate::inference::*;
use crate::cmdline::AnalyzeArgs;

use clap::error;
use simple_logger::SimpleLogger;
use log::{info, warn, error};
use glob::glob;
use std::fs;

pub fn analyze(args: AnalyzeArgs) {
    SimpleLogger::new().with_level(log::LevelFilter::Info).init().unwrap();
    // [TODO] Multithreaded version is under development.
    //rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

    //info!("Using {} threads for analysis.", args.threads);

    // Expand globs and categorize files before processing so we can auto-determine -c
    let mut raw_files: Vec<String> = Vec::new();
    let mut sketch_files: Vec<String> = Vec::new();
    for file in &args.files {
        for entry in glob(file).expect("Failed to read glob pattern") {
            match entry {
                Ok(path) => {
                    let file_str = path.to_str().unwrap().to_string();
                    if is_fastx_file(&file_str) {
                        raw_files.push(file_str);
                    } else if is_sketch_file(&file_str) {
                        sketch_files.push(file_str);
                    } else {
                        warn!("File format not recognized for file: {}. Skipping.", file_str);
                    }
                }
                Err(e) => warn!("Error reading file: {:?}", e),
            }
        }
    }

    let c = args.c.unwrap_or_else(|| {
        let raw_refs: Vec<&str> = raw_files.iter().map(|s| s.as_str()).collect();
        let (auto_c, est_file_size) = estimate_c_from_raw_files(&raw_refs);
        info!("Total estimated input sequence file size (decompressed): {:.2} GB", est_file_size as f64 / (1024.0 * 1024.0 * 1024.0));
        info!("Auto-determined subsampling rate: -c {}", auto_c);
        auto_c
    });

    let mut kvmer_set = KVmerSet::new(args.k, args.v, !args.forward_only);

    // Read query files
    info!("Processing query files...");
    for file_str in &raw_files {
        kvmer_set.add_file_to_kvmer_set(file_str, c, args.trim_front, args.trim_back);
    }
    for file_str in &sketch_files {
        kvmer_set.load(file_str);
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
        reference_kvmer_set.add_file_to_kvmer_set(reference, c, args.trim_front, args.trim_back);
        info!("Loaded reference file: {}", reference);

        stats = kvmer_set.get_stats_with_reference(lower_bound, &reference_kvmer_set, args.first_base_only);
    } else {
        let lower_bound = args.lower_bound.unwrap_or(10);
        //println!("Error rate: {}", kvmer_set.get_stats(args.threshold));
        stats = kvmer_set.get_stats(lower_bound, args.first_base_only);
    }
    // if reference is set, the filter should be disabled
    // [FIXME] enable --use-all by default
    if args.reference.is_some() && !args.use_all {
        warn!("If reference is provided, --use-all is recommended.");
    }

    let spectrum = analyzer.analyze(&stats);
    let analysis_output = format!("{}\n{}", header_str(!args.forward_only), spectrum_to_str(&spectrum, !args.forward_only));

    if let Some(prefix) = &args.output_prefix {
        fs::write(format!("{}.summary_error_rate.csv", prefix), &analysis_output).unwrap();
        info!("Output written to prefix {}.", prefix);
    } else {
        error!("No output prefix provided. Use -o or --output-prefix to specify the output file prefix for the analysis results.");
    }
}