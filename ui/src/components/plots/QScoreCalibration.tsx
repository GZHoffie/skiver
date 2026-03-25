import { useState } from "react";
import Plot from "react-plotly.js";
import { useCSVData } from "../../hooks/useCSVData";

interface PhredRow {
  qscore: number;
  empirical_qscore: number;
  num_correct: number;
  num_error: number;
  error_rate: number;
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

  const filtered = rows.filter(
    (r) => r.num_correct + r.num_error >= MIN_COVERAGE
  );
  if (filtered.length === 0)
    return <div className="plot-loading">No data above coverage threshold.</div>;

  const qVals = filtered.map((r) => r.qscore);
  const empRate = filtered.map((r) => r.error_rate);
  const counts = filtered.map((r) => r.num_correct + r.num_error);

  const maxQ = Math.max(...qVals);
  const qTheory = Array.from({ length: 300 }, (_, i) => 1 + (i * (maxQ - 1)) / 299);
  const theoryRate = qTheory.map((q) => Math.pow(10, -q / 10));

  const traces: Plotly.Data[] = [
    {
      type: "scatter",
      x: qTheory,
      y: theoryRate,
      mode: "lines",
      name: "Theoretical (10^(-Q/10))",
      line: { color: "black", dash: "dash", width: 4 },
      xaxis: "x",
      yaxis: "y",
    },
    {
      type: "scatter",
      x: qVals,
      y: empRate,
      mode: "lines+markers",
      name: "Empirical error rate",
      line: { color: "slategray", width: 4 },
      marker: { size: 10 },
      xaxis: "x",
      yaxis: "y",
    },
    {
      type: "bar",
      x: qVals,
      y: counts,
      name: "Number of bases used for estimation",
      marker: { color: "slategray", opacity: 0.6 },
      xaxis: "x2",
      yaxis: "y2",
    },
  ];

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    xaxis:  { title: "",              domain: [0, 1],     anchor: "y",  showticklabels: false },
    yaxis:  {                         domain: [0.32, 1],  anchor: "x",  type: logScale ? "log" : "linear" },
    xaxis2: { title: { text: "Phred Q-score" }, domain: [0, 1],     anchor: "y2" },
    yaxis2: {                         domain: [0, 0.28],  anchor: "x2" },
    legend: { x: 1, xanchor: "right", y: 1 },
    font: { family: "DejaVu Sans" },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 30, l: 80, r: 20, b: 50 },
    annotations: [
      { text: "Error rate", xref: "paper", yref: "paper", x: -0.06, y: 0.66, showarrow: false, textangle: -90, font: { size: 14 } },
      { text: "Count",      xref: "paper", yref: "paper", x: -0.06, y: 0.14, showarrow: false, textangle: -90, font: { size: 14 } },
    ],
  };

  return (
    <div style={{ width: "100%" }}>
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
        style={{ width: "100%", height: "500px" }}
      />
    </div>
  );
}
