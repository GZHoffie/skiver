use crate::kvmer::*;
use crate::utils::*;
use crate::types::*;
use crate::inference::*;
use crate::cmdline::AnalyzeArgs;

use simple_logger::SimpleLogger;
use log::{info, warn};
use glob::glob;

pub fn analyze(args: AnalyzeArgs) {
    SimpleLogger::new().with_level(log::LevelFilter::Info).init().unwrap();
    // [TODO] Multithreaded version is under development.
    //rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

    //info!("Using {} threads for analysis.", args.threads);

    let mut kvmer_set = KVmerSet::new(args.k, args.v, !args.forward_only);
    
    // Read query files
    info!("Processing query files...");
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
                        warn!("File format not recognized for file: {}. Skipping.", file_str);
                    }
                }
                Err(e) => warn!("Error reading file: {:?}", e),
            }
        }
    }
    info!("Finished processing query files.");

    let analyzer = ErrorAnalyzer::new(args.clone());

    
    let stats: KVmerStats;
    if let Some(reference) = &args.reference {
        let lower_bound = args.lower_bound.unwrap_or(0);

        let mut reference_kvmer_set = KVmerSet::new(args.k, args.v, true);
        reference_kvmer_set.add_file_to_kvmer_set(reference, args.c, args.trim_front, args.trim_back);
        info!("Loaded reference file: {}", reference);

        stats = kvmer_set.get_stats_with_reference(lower_bound, &reference_kvmer_set);
    } else {
        let lower_bound = args.lower_bound.unwrap_or(10);
        //println!("Error rate: {}", kvmer_set.get_stats(args.threshold));
        stats = kvmer_set.get_stats(lower_bound);
    }
    if let Some(output_path) = &args.output_path {
        kvmer_set.output_stats(output_path, &stats, true, true);
    }
    let spectrum = analyzer.analyze(&stats);

    // output to stdout a csv file
    // [TODO] allow output a separate line per file
    println!("{}", header_str(!args.forward_only));
    let spectrum_str = spectrum_to_str(&spectrum, !args.forward_only);
    println!("{}", spectrum_str);

    
    
    /*
    let mut vmer_set = VmerSet::new(args.v, false);
    for file in &args.files {
        vmer_set.add_file_first_pass(file, args.c);
    }
    let relevant_values = vmer_set._get_relevant_values(args.threshold);
    let mut value_counts: HashMap<u64, u32> = HashMap::new();
    for file in &args.files {
        vmer_set.add_file_second_pass(file, &mut value_counts, &relevant_values);
    }
    vmer_set.add_value_counts(&value_counts);
    let stats = vmer_set.get_stats(args.threshold);

    let spectrum = error_profile(&stats, false);
    output_error_spectrum(&spectrum, args.v);
    */
}