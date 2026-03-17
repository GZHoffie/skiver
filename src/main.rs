use clap::Parser;


use skiver::cmdline::*;
use skiver::analyze;
use skiver::sketch;
use skiver::mapping;


//Use this allocator when statically compiling
//instead of the default
//because the musl statically compiled binary
//uses a bad default allocator which makes the
//binary take 60% longer!!! Only affects
//static compilation though. 
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
//    set_hook(Box::new(|info| {
//        if let Some(s) = info.payload().downcast_ref::<String>() {
//            log::error!("{}", s);
//        }
//    }));

    let cli = Cli::parse();
    match cli.mode {
        Mode::Sketch(sketch_args) => sketch::sketch(sketch_args),
        Mode::Analyze(analyze_args) => analyze::analyze(analyze_args),
        Mode::Map(map_args) => mapping::map(map_args),
    }
}