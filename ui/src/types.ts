export type AppPhase = "idle" | "sketching" | "analyzing" | "done" | "error";

export interface SketchParams {
  input_files: string[];
  output_path: string;
  k: number;
  v: number;
  c: number | null;
  trim_front: number;
  trim_back: number;
  forward_only: boolean;
}

export interface AnalyzeParams {
  input_files: string[];
  output_prefix: string;
  k: number;
  v: number;
  lower_bound: number | null;
  forward_only: boolean;
  use_all: boolean;
  outlier_threshold: number;
  ignore_largest_t: number;
  ignore_smallest_t: number;
  hazard_model: string;
  num_experiments: number;
  reference: string | null;
}

export interface LogLine {
  stream: "stdout" | "stderr" | "exit";
  line: string;
  exit_code?: number;
}

export interface CsvPaths {
  summaryErrorRate: string;
  kvmer: string;
  hazardRate: string;
  survivalRate: string;
  summaryErrorSpectrum: string;
  summaryErrorSpectrumDepT: string;
  summaryPhred: string | null;
  summaryGcContent: string;
  summaryReadPosition: string;
}

export interface WorkflowState {
  inputFiles: string[];
  outputDir: string;
  outputPrefix: string;
  sketchParams: SketchParams;
  analyzeParams: AnalyzeParams;
  phase: AppPhase;
  logs: LogLine[];
  csvPaths: CsvPaths | null;
  errorMessage: string | null;
}
