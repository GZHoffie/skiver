use simple_logger::SimpleLogger;
use log::info;
use glob::glob;

use crate::cmdline::CalibrateArgs;
use crate::kvmer::KVmerSet;
use crate::utils::{is_fastx_file, is_sketch_file};

pub fn calibrate(args: CalibrateArgs) {
    SimpleLogger::new().with_level(log::LevelFilter::Info).init().unwrap();

    let mut kvmer_set = KVmerSet::new(args.k, args.v, !args.forward_only);

    info!("Processing input files...");
    for file in &args.files {
        for entry in glob(file).expect("Failed to read glob pattern") {
            match entry {
                Ok(path) => {
                    let file_str = path.to_str().unwrap();
                    if is_fastx_file(file_str) {
                        kvmer_set.add_file_to_kvmer_set(file_str, args.c, args.trim_front, args.trim_back);
                    } else if is_sketch_file(file_str) {
                        kvmer_set.load(file_str);
                    } else {
                        log::warn!("File format not recognized: {}. Skipping.", file_str);
                    }
                }
                Err(e) => log::warn!("Error reading file: {:?}", e),
            }
        }
    }
    info!("Finished processing input files.");

    let stats = if let Some(reference) = &args.reference {
        let mut ref_kvmer_set = KVmerSet::new(args.k, args.v, true);
        ref_kvmer_set.add_file_to_kvmer_set(reference, args.c, args.trim_front, args.trim_back);
        info!("Loaded reference: {}", reference);
        kvmer_set.get_stats_with_reference(args.lower_bound, &ref_kvmer_set)
    } else {
        kvmer_set.get_stats(args.lower_bound)
    };

    // Collect all observed qscores from both maps.
    let mut qscores: Vec<u8> = stats.qscore_correct.keys()
        .chain(stats.qscore_error.keys())
        .cloned()
        .collect();
    qscores.sort_unstable();
    qscores.dedup();

    println!("qscore,num_correct,num_error,error_rate");
    for q in qscores {
        let correct = *stats.qscore_correct.get(&q).unwrap_or(&0);
        let error   = *stats.qscore_error.get(&q).unwrap_or(&0);
        let total   = correct + error;
        let error_rate = if total > 0 { error as f64 / total as f64 } else { 0.0 };
        println!("{},{},{},{:.6}", q, correct, error, error_rate);
    }
}
