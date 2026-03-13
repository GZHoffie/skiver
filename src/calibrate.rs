use simple_logger::SimpleLogger;
use log::info;
use glob::glob;

use crate::cmdline::{AnalyzeArgs, CalibrateArgs};
use crate::inference::ErrorAnalyzer;
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

    let analyze_args = AnalyzeArgs {
        k: args.k,
        v: args.v,
        outlier_threshold: args.outlier_threshold,
        num_experiments: args.num_experiments,
        use_all: args.use_all,
        ignore_last_hazard_ratios: 2,
        estimation_method: "sum_ratio".to_string(),
        hazard_model: "weibull".to_string(),
        ..Default::default()
    };
    let analyzer = ErrorAnalyzer::new(analyze_args);

    let results = analyzer.calibrate_qscores(&stats);

    println!("qscore,empirical_qscore,num_correct,num_error,error_rate,5th_percentile,95th_percentile");
    for (q, correct, error, error_rate, lower, upper) in results {
        let empirical_q = if error_rate > 0.0 { -10.0 * error_rate.log10() } else { f64::INFINITY };
        println!("{},{:.4},{},{},{:.6},{:.6},{:.6}", q, empirical_q, correct, error, error_rate, lower, upper);
    }
}
