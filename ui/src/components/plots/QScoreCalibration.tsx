import { useState } from "react";
import Plot from "react-plotly.js";
import { useCSVData } from "../../hooks/useCSVData";

interface PhredRow {
  qscore: number;
  empirical_qscore: number;
  num_observed: number;
  num_substitution: number;
  num_insertion: number;
  num_deletion: number;
  bayesian_substitution_rate: number;
  bayesian_insertion_rate: number;
  bayesian_deletion_rate: number;
  bayesian_error_rate: number;
}

interface Props {
  phredPath: string;
}

const MIN_COVERAGE = 100;

export function QScoreCalibration({ phredPath }: Props) {
  const { data: rows, loading } = useCSVData<PhredRow>(phredPath);
  const [logScale, setLogScale] = useState(true);

  if (loading) return <div className="plot-loading">Loading…</div>;
  if (!rows) return null;

  const filtered = rows.filter((r) => r.num_observed >= MIN_COVERAGE);
  if (filtered.length === 0)
    return <div className="plot-loading">No data above coverage threshold.</div>;

  const qVals   = filtered.map((r) => r.qscore);
  const subRate = filtered.map((r) => r.bayesian_substitution_rate);
  const insRate = filtered.map((r) => r.bayesian_insertion_rate);
  const delRate = filtered.map((r) => r.bayesian_deletion_rate);
  const errRate = filtered.map((r) => r.bayesian_error_rate);
  const counts  = filtered.map((r) => r.num_observed);

  const maxQ = Math.max(...qVals);
  const minQ = Math.min(...qVals);
  const qTheory = Array.from({ length: 300 }, (_, i) => minQ + (i * (maxQ - minQ)) / 299);
  const theoryRate = qTheory.map((q) => Math.pow(10, -q / 10));

  const traces: Plotly.Data[] = [
    // ── Subplot 1: substitution / insertion / deletion rates ──
    {
      type: "scatter",
      x: qVals,
      y: subRate,
      mode: "lines+markers",
      name: "Substitution",
      line: { color: "steelblue", width: 3 },
      marker: { size: 6, symbol: "circle" },
      xaxis: "x",
      yaxis: "y",
    },
    {
      type: "scatter",
      x: qVals,
      y: insRate,
      mode: "lines+markers",
      name: "Insertion",
      line: { color: "darkorange", width: 3 },
      marker: { size: 6, symbol: "square" },
      xaxis: "x",
      yaxis: "y",
    },
    {
      type: "scatter",
      x: qVals,
      y: delRate,
      mode: "lines+markers",
      name: "Deletion",
      line: { color: "forestgreen", width: 3 },
      marker: { size: 6, symbol: "triangle-up" },
      xaxis: "x",
      yaxis: "y",
    },
    // ── Subplot 2: theoretical vs. empirical (Bayesian) error rate ──
    {
      type: "scatter",
      x: qTheory,
      y: theoryRate,
      mode: "lines",
      name: "Theoretical (10^(-Q/10))",
      line: { color: "black", dash: "dash", width: 3 },
      xaxis: "x2",
      yaxis: "y2",
    },
    {
      type: "scatter",
      x: qVals,
      y: errRate,
      mode: "lines+markers",
      name: "Normalized error rate",
      line: { color: "slategray", width: 3 },
      marker: { size: 6 },
      xaxis: "x2",
      yaxis: "y2",
    },
    // ── Subplot 3: histogram of observed bases ──
    {
      type: "bar",
      x: qVals,
      y: counts,
      name: "# observed bases",
      marker: { color: "slategray", opacity: 0.6 },
      xaxis: "x3",
      yaxis: "y3",
    },
  ];

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    xaxis:  { title: "",                                                  domain: [0, 1], anchor: "y",  showticklabels: false },
    yaxis:  { title: { text: "Normalized rate",    font: { size: 14 } }, domain: [0.62, 1],    anchor: "x",  type: logScale ? "log" : "linear" },
    xaxis2: { title: "",                                                  domain: [0, 1], anchor: "y2", showticklabels: false, matches: "x" },
    yaxis2: { title: { text: "Error rate",         font: { size: 14 } }, domain: [0.24, 0.58], anchor: "x2", type: logScale ? "log" : "linear" },
    xaxis3: { title: { text: "Phred Q-score",      font: { size: 14 } }, domain: [0, 1], anchor: "y3", matches: "x" },
    yaxis3: { title: { text: "# observed bases",   font: { size: 14 } }, domain: [0, 0.20],   anchor: "x3" },
    annotations: [
      { text: "Error component rates by quality score", xref: "paper", yref: "paper", x: 0.5, y: 1.02, showarrow: false, xanchor: "center", font: { size: 14 } },
      { text: "Quality-score calibration",              xref: "paper", yref: "paper", x: 0.5, y: 0.60, showarrow: false, xanchor: "center", font: { size: 14 } },
    ],
    legend: { x: 1, xanchor: "right", y: 1, font: { size: 13 } },
    font: { family: "DejaVu Sans", size: 13 },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 40, l: 80, r: 20, b: 50 },
  };

  return (
    <div style={{ width: "100%" }}>
      <p style={{ fontSize: "13px", color: "#555", fontStyle: "italic", margin: "0 20px 8px 20px" }}>
        Rates are Bayesian-estimated given the Phred quality score, conditioned on the previous
        (<em>t</em>&#8209;1) bases agreeing with the consensus.
      </p>
      <div style={{ display: "flex", justifyContent: "flex-end", padding: "0 20px 8px 0" }}>
        <button
          onClick={() => setLogScale((v) => !v)}
          style={{ fontSize: "13px", padding: "3px 10px", cursor: "pointer" }}
        >
          {logScale ? "Switch to Linear scale" : "Switch to Log scale"}
        </button>
      </div>
      <Plot
        data={traces}
        layout={layout}
        useResizeHandler
        style={{ width: "100%", height: "1100px" }}
      />
    </div>
  );
}
