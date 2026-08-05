import { useState } from "react";
import { SketchParams, AnalyzeParams } from "./types";
import { useSkiverRun } from "./hooks/useSkiverRun";
import { FileSelector } from "./components/workflow/FileSelector";
import { ReferenceSelector } from "./components/workflow/ReferenceSelector";
import { OutputConfig } from "./components/workflow/OutputConfig";
import { ParameterPanel } from "./components/workflow/ParameterPanel";
import { LogConsole } from "./components/workflow/LogConsole";
import { HazardSurvivalPlot } from "./components/plots/HazardSurvivalPlot";
import { ErrorSpectrumHeatmap } from "./components/plots/ErrorSpectrumHeatmap";
import { SBS96Plot } from "./components/plots/SBS96Plot";
import { SpectrumDependenceOnT } from "./components/plots/SpectrumDependenceOnT";
import { QScoreCalibration } from "./components/plots/QScoreCalibration";
import { GCContentPlot } from "./components/plots/GCContentPlot";
import { ReadPositionPlot } from "./components/plots/ReadPositionPlot";
import { CoverageHistogram } from "./components/plots/CoverageHistogram";
import { SummaryBoxes } from "./components/SummaryBoxes";
import "./App.css";

const DEFAULT_SKETCH: SketchParams = {
  input_files: [],
  output_path: "",
  k: 17,
  v: 17,
  c: null,
  trim_front: 0,
  trim_back: 0,
  forward_only: false,
};

const DEFAULT_ANALYZE: AnalyzeParams = {
  input_files: [],
  output_prefix: "",
  k: 17,
  v: 17,
  lower_bound: null,
  forward_only: false,
  use_all: false,
  outlier_threshold: 1e-12,
  ignore_largest_t: 2,
  ignore_smallest_t: 2,
  hazard_model: "weibull",
  num_experiments: 100,
  reference: null,
};

type Tab = "run" | "results";

const PLOT_TABS = [
  "Hazard & Survival",
  "Error Spectrum",
  "SBS96",
  "Spectrum vs t",
  "Q-Score",
  "GC Content",
  "Read Position",
  "Coverage",
] as const;

export default function App() {
  const [tab, setTab] = useState<Tab>("run");
  const [plotTab, setPlotTab] = useState<(typeof PLOT_TABS)[number]>(PLOT_TABS[0]);
  const [inputFiles, setInputFiles] = useState<string[]>([]);
  const [outputDir, setOutputDir] = useState("");
  const [outputPrefix, setOutputPrefix] = useState("sample");
  const [sketchParams, setSketchParams] = useState<SketchParams>(DEFAULT_SKETCH);
  const [analyzeParams, setAnalyzeParams] = useState<AnalyzeParams>(DEFAULT_ANALYZE);

  const { phase, logs, csvPaths, errorMessage, run } = useSkiverRun();

  const busy = phase === "sketching" || phase === "analyzing";

  const [errors, setErrors] = useState({ inputFiles: false, outputDir: false, outputPrefix: false });

  function handleRun() {
    const e = {
      inputFiles: inputFiles.length === 0,
      outputDir: !outputDir,
      outputPrefix: !outputPrefix.trim(),
    };
    setErrors(e);
    if (e.inputFiles || e.outputDir || e.outputPrefix) return;

    const kvmerPath = `${outputDir}/${outputPrefix}.kvmer`;
    const fullPrefix = `${outputDir}/${outputPrefix}`;

    const sp: SketchParams = {
      ...sketchParams,
      input_files: inputFiles,
      output_path: kvmerPath,
    };
    const ap: AnalyzeParams = {
      ...analyzeParams,
      input_files: [kvmerPath],
      output_prefix: fullPrefix,
    };

    run(sp, ap).then(() => setTab("results"));
  }

  const phaseLabel: Record<string, string> = {
    idle: "",
    sketching: "Running skiver sketch",
    analyzing: "Running skiver analyze",
    done: "Done",
    error: "Error",
  };

  return (
    <div className="app">
      <header className="app-header">
        <h1>Skiver</h1>
        <nav className="tab-bar">
          <button
            className={tab === "run" ? "active" : ""}
            onClick={() => setTab("run")}
          >
            Run
          </button>
          <button
            className={tab === "results" ? "active" : ""}
            onClick={() => setTab("results")}
            disabled={!csvPaths}
          >
            Results
          </button>
        </nav>
        {busy && <span className="phase-badge busy">{phaseLabel[phase]}</span>}
        {phase === "done" && <span className="phase-badge done">Done</span>}
        {phase === "error" && <span className="phase-badge error">Error</span>}
      </header>

      {tab === "run" && (
        <main className="run-panel">
          <FileSelector
            files={inputFiles}
            onChange={(f) => { setInputFiles(f); if (f.length > 0) setErrors((e) => ({ ...e, inputFiles: false })); }}
            disabled={busy}
            error={errors.inputFiles}
          />
          <ReferenceSelector
            reference={analyzeParams.reference}
            onChange={(ref) => setAnalyzeParams((p) => ({ ...p, reference: ref, use_all: ref !== null ? true : p.use_all }))}
            disabled={busy}
          />
          <OutputConfig
            outputDir={outputDir}
            outputPrefix={outputPrefix}
            onDirChange={(d) => { setOutputDir(d); if (d) setErrors((e) => ({ ...e, outputDir: false })); }}
            onPrefixChange={(p) => { setOutputPrefix(p); if (p.trim()) setErrors((e) => ({ ...e, outputPrefix: false })); }}
            disabled={busy}
            dirError={errors.outputDir}
            prefixError={errors.outputPrefix}
          />
          <ParameterPanel
            sketch={sketchParams}
            analyze={analyzeParams}
            onSketchChange={setSketchParams}
            onAnalyzeChange={setAnalyzeParams}
            disabled={busy}
          />

          <div className="run-controls">
            {phase === "idle" || phase === "error" ? (
              <button className="run-btn" onClick={handleRun}>
                Run Skiver
              </button>
            ) : phase === "done" ? (
              <button className="run-btn" onClick={handleRun}>
                Run again
              </button>
            ) : (
              <button className="run-btn" disabled>
                {phaseLabel[phase]}
              </button>
            )}
          </div>

          {errorMessage && (
            <div className="error-banner">{errorMessage}</div>
          )}
          <LogConsole logs={logs} />
        </main>
      )}

      {tab === "results" && csvPaths && (
        <main className="results-panel">
          <SummaryBoxes errorRatePath={csvPaths.summaryErrorRate} />
          <nav className="plot-tab-bar">
            {PLOT_TABS.map((pt) => (
              <button
                key={pt}
                className={plotTab === pt ? "active" : ""}
                onClick={() => setPlotTab(pt)}
                disabled={pt === "Q-Score" && !csvPaths.summaryPhred}
                title={
                  pt === "Q-Score" && !csvPaths.summaryPhred
                    ? "No Q-score calibration data was produced"
                    : undefined
                }
              >
                {pt}
              </button>
            ))}
          </nav>

          <div className="plot-area">
            {plotTab === "Hazard & Survival" && (
              <HazardSurvivalPlot
                hazardPath={csvPaths.hazardRate}
                errorRatePath={csvPaths.summaryErrorRate}
                survivalPath={csvPaths.survivalRate}
              />
            )}
            {plotTab === "Error Spectrum" && (
              <ErrorSpectrumHeatmap
                spectrumPath={csvPaths.summaryErrorSpectrum}
                errorRatePath={csvPaths.summaryErrorRate}
              />
            )}
            {plotTab === "SBS96" && (
              <SBS96Plot spectrumPath={csvPaths.summaryErrorSpectrum} />
            )}
            {plotTab === "Spectrum vs t" && (
              <SpectrumDependenceOnT depTPath={csvPaths.summaryErrorSpectrumDepT} />
            )}
            {plotTab === "Q-Score" && (
              <QScoreCalibration phredPath={csvPaths.summaryPhred} />
            )}
            {plotTab === "GC Content" && (
              <GCContentPlot gcContentPath={csvPaths.summaryGcContent} />
            )}
            {plotTab === "Read Position" && (
              <ReadPositionPlot readPositionPath={csvPaths.summaryReadPosition} />
            )}
            {plotTab === "Coverage" && (
              <CoverageHistogram
                kvmerPath={csvPaths.kvmer}
                errorRatePath={csvPaths.summaryErrorRate}
              />
            )}
          </div>
        </main>
      )}
    </div>
  );
}
