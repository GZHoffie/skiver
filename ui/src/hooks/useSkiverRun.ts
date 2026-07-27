import { useState, useCallback, useRef } from "react";
import { invoke } from "@tauri-apps/api/core";
import { listen, UnlistenFn } from "@tauri-apps/api/event";
import {
  AppPhase,
  AnalyzeParams,
  CsvPaths,
  LogLine,
  SketchParams,
} from "../types";

export function useSkiverRun() {
  const [phase, setPhase] = useState<AppPhase>("idle");
  const [logs, setLogs] = useState<LogLine[]>([]);
  const [csvPaths, setCsvPaths] = useState<CsvPaths | null>(null);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);
  const unlistenRef = useRef<UnlistenFn | null>(null);
  const lastSketchRef = useRef<SketchParams | null>(null);

  const appendLog = useCallback((entry: LogLine) => {
    setLogs((prev) => [...prev, entry]);
  }, []);

  function isSequenceInput(path: string): boolean {
    const lower = path.toLowerCase();
    return [
      ".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz",
      ".fq", ".fnq", ".fastq", ".fq.gz", ".fnq.gz", ".fastq.gz",
      ".bam",
    ].some((extension) => lower.endsWith(extension));
  }

  const run = useCallback(
    async (sketchParams: SketchParams, analyzeParams: AnalyzeParams) => {
      setLogs([]);
      setCsvPaths(null);
      setErrorMessage(null);

      const sequenceInputs = sketchParams.input_files.filter(isSequenceInput);
      const sketchInputs = sketchParams.input_files.filter(
        (path) => !isSequenceInput(path)
      );
      if (
        sketchInputs.length > 1 ||
        (sketchInputs.length === 1 && sequenceInputs.length > 0)
      ) {
        setErrorMessage(
          "Choose either exactly one pre-computed sketch file or one or more FASTA/FASTQ/BAM files."
        );
        setPhase("error");
        return;
      }

      const shouldRunSketch = sketchInputs.length === 0;
      const analyzeInputFiles = shouldRunSketch
        ? analyzeParams.input_files
        : sketchInputs;
      const effectiveAnalyzeParams: AnalyzeParams = {
        ...analyzeParams,
        input_files: analyzeInputFiles,
      };
      const sketchChanged =
        !lastSketchRef.current ||
        JSON.stringify(sketchParams) !== JSON.stringify(lastSketchRef.current);

      // Set up log listener
      if (unlistenRef.current) unlistenRef.current();
      unlistenRef.current = await listen<LogLine>("skiver-log", (event) => {
        appendLog(event.payload);
      });

      try {
        // --- Sketch ---
        if (!shouldRunSketch) {
          appendLog({ stream: "stdout", line: "=== skiver sketch (skipped — pre-computed sketch provided) ===" });
        } else if (sketchChanged) {
          setPhase("sketching");
          appendLog({ stream: "stdout", line: "=== skiver sketch ===" });
          await invoke("run_sketch", { params: sketchParams });
          lastSketchRef.current = sketchParams;
        } else {
          appendLog({ stream: "stdout", line: "=== skiver sketch (skipped — params unchanged) ===" });
        }

        // --- Analyze ---
        setPhase("analyzing");
        appendLog({ stream: "stdout", line: "=== skiver analyze ===" });
        await invoke("run_analyze", { params: effectiveAnalyzeParams });

        // Build CSV paths
        const prefix = effectiveAnalyzeParams.output_prefix;
        const paths: CsvPaths = {
          summaryErrorRate: `${prefix}.summary_error_rate.csv`,
          kvmer: `${prefix}.kvmer.csv`,
          hazardRate: `${prefix}.hazard_rate.csv`,
          survivalRate: `${prefix}.survival_rate.csv`,
          summaryErrorSpectrum: `${prefix}.summary_error_spectrum.csv`,
          summaryErrorSpectrumDepT: `${prefix}.summary_error_spectrum_dependence_on_t.csv`,
          summaryPhred: `${prefix}.summary_phred.csv`,
          summaryGcContent: `${prefix}.summary_gc_content.csv`,
          summaryReadPosition: `${prefix}.summary_read_position.csv`,
        };

        setCsvPaths(paths);
        setPhase("done");
      } catch (err) {
        setErrorMessage(String(err));
        setPhase("error");
      } finally {
        if (unlistenRef.current) {
          unlistenRef.current();
          unlistenRef.current = null;
        }
      }
    },
    [appendLog]
  );

  const reset = useCallback(() => {
    setPhase("idle");
    setLogs([]);
    setCsvPaths(null);
    setErrorMessage(null);
  }, []);

  return { phase, logs, csvPaths, errorMessage, run, reset };
}
