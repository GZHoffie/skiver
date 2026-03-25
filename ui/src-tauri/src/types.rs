use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct SketchParams {
    pub input_files: Vec<String>,
    pub output_path: String,
    pub k: u8,
    pub v: u8,
    pub c: Option<usize>,
    pub trim_front: usize,
    pub trim_back: usize,
    pub forward_only: bool,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct AnalyzeParams {
    pub kvmer_path: String,
    pub output_prefix: String,
    pub k: u8,
    pub v: u8,
    pub lower_bound: Option<u32>,
    pub forward_only: bool,
    pub use_all: bool,
    pub outlier_threshold: f64,
    pub ignore_largest_t: usize,
    pub ignore_smallest_t: usize,
    pub hazard_model: String,
    pub num_experiments: u32,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct LogEvent {
    pub stream: String,
    pub line: String,
    pub exit_code: Option<i32>,
}
