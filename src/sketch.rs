

use crate::kvmer::*;
use simple_logger::SimpleLogger;
use log::info;
use crate::cmdline::SketchArgs;
//use rayon::prelude::*;


pub fn sketch(args: SketchArgs) {
    SimpleLogger::new().with_level(log::LevelFilter::Info).init().unwrap();

    // check if the output file is valid
    if args.output_path.is_empty() {
        panic!("Output file path is empty.");
    }
    let output = std::path::Path::new(&args.output_path);
    let parent = output.parent().unwrap_or(std::path::Path::new("."));
    if !parent.exists() {
        panic!("Output directory '{}' does not exist.", parent.display());
    }
    // verify we can actually create/write the output file before doing any work
    std::fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(false)
        .open(output)
        .unwrap_or_else(|e| panic!("Cannot write to output path '{}': {}", args.output_path, e));
    std::fs::remove_file(output).ok();

    //rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
    info!("Processing query files...");
    
    let mut kvmer_set = KVmerSet::new(args.k, args.v, !args.forward_only);
    for file in &args.files {
        kvmer_set.add_file_to_kvmer_set(file, args.c, args.trim_front, args.trim_back);
    }
    info!("Finished processing query files.");

    kvmer_set.dump(&args.output_path);
    info!("Sketch saved to {}", args.output_path);
}