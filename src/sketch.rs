

use crate::kvmer::*;
use crate::utils::estimate_c_from_raw_files;
use simple_logger::SimpleLogger;
use log::{info, warn, error};
use crate::cmdline::SketchArgs;
use std::fs::{self, OpenOptions};
use std::path::Path;
//use rayon::prelude::*;

fn validate_output_path(output_path: &str) -> Result<(Option<String>, bool), String> {
    if output_path.is_empty() {
        return Err("Output file path is empty.".to_string());
    }

    let output = Path::new(output_path);
    if output.file_name().is_none() {
        return Err(format!("Output path '{}' does not contain a file name.", output_path));
    }

    let parent = output.parent().filter(|path| !path.as_os_str().is_empty())
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
        .map_err(|e| format!("Output path '{}' is not writable: {}", output_path, e))?;
    fs::remove_file(&probe)
        .map_err(|e| format!("Could not remove output validation file '{}': {}", probe.display(), e))?;

    Ok((created_dir, output.exists()))
}

pub fn sketch(args: SketchArgs) {
    SimpleLogger::new().with_level(log::LevelFilter::Info).init().unwrap();

    let (created_dir, output_exists) = match validate_output_path(&args.output_path) {
        Ok(result) => result,
        Err(message) => {
            error!("{}", message);
            std::process::exit(1);
        }
    };
    if let Some(directory) = created_dir {
        info!("Created output directory '{}'.", directory);
    }
    if output_exists {
        warn!("Output file '{}' already exists and will be overwritten.", args.output_path);
    }

    //rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

    let c = args.c.unwrap_or_else(|| {
        let raw_refs: Vec<&str> = args.files.iter().map(|s| s.as_str()).collect();
        let (auto_c, est_file_size) = estimate_c_from_raw_files(&raw_refs);
        info!("Total estimated input sequence file size (decompressed): {:.2} GB", est_file_size as f64 / (1024.0 * 1024.0 * 1024.0));
        info!("Auto-determined subsampling rate: -c {}", auto_c);
        auto_c
    });

    info!("Processing query files...");

    let mut kvmer_set = KVmerSet::new(args.k, args.v, !args.forward_only);
    let threads = if args.threads == 0 { std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1) } else { args.threads };
    info!("Using {} FASTA/FASTQ worker thread(s).", threads);
    for file in &args.files {
        kvmer_set.add_file_to_kvmer_set_with_threads(file, c, args.trim_front, args.trim_back, threads);
    }
    info!("Finished processing query files.");

    kvmer_set.dump(&args.output_path);
    info!("Sketch saved to {}", args.output_path);
}
