import { useState } from "react";
import Plot from "react-plotly.js";
import { useCSVData } from "../../hooks/useCSVData";

interface PhredRow {
  qscore: number;
  empirical_qscore: number;
  num_correct: number;
  num_error: number;
  per_base_error_rate?: number;
  per_base_error_rate_median?: number;
  "per_base_error_rate_5-95th_percentile"?: string;
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
  const empRate = filtered.map((r) => r.per_base_error_rate ?? r.per_base_error_rate_median ?? 0);
  const counts = filtered.map((r) => r.num_correct + r.num_error);
  const ci = filtered.map((r) => {
    const [lo, hi] = (r["per_base_error_rate_5-95th_percentile"] ?? "").split("~");
    return { lo: Number(lo), hi: Number(hi) };
  });
  const ciLower = ci.map((v) => v.lo);
  const ciUpper = ci.map((v) => v.hi);
  const hasCi = ci.every((v) => Number.isFinite(v.lo) && Number.isFinite(v.hi));

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
    ...(hasCi
      ? [
          {
            type: "scatter" as const,
            x: qVals,
            y: ciLower,
            mode: "lines" as const,
            line: { color: "transparent" },
            showlegend: false,
            hoverinfo: "skip" as const,
            xaxis: "x",
            yaxis: "y",
          },
          {
            type: "scatter" as const,
            x: qVals,
            y: ciUpper,
            mode: "lines" as const,
            fill: "tonexty" as const,
            fillcolor: "rgba(112,128,144,0.25)",
            line: { color: "transparent" },
            name: "5-95% percentile",
            xaxis: "x",
            yaxis: "y",
          },
        ]
      : []),
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
    yaxis:  { title: { text: "Error rate", font: { size: 14 } }, domain: [0.32, 1], anchor: "x", type: logScale ? "log" : "linear" },
    xaxis2: { title: { text: "Phred Q-score", font: { size: 14 } }, domain: [0, 1], anchor: "y2", matches: "x" },
    yaxis2: { title: { text: "Count", font: { size: 14 } }, domain: [0, 0.28], anchor: "x2" },
    legend: { x: 1, xanchor: "right", y: 1, font: { size: 14 } },
    font: { family: "DejaVu Sans", size: 13 },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 30, l: 80, r: 20, b: 50 },
  };

  return (
    <div style={{ width: "100%" }}>
      <p style={{ fontSize: "13px", color: "#555", fontStyle: "italic", margin: "0 20px 8px 20px" }}>
        The empirical error rate is calculated to be the probability that the base disagrees with the
        consensus, given that the previous (<em>t</em>&#8209;1) bases agree.
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
        style={{ width: "100%", height: "750px" }}
      />
    </div>
  );
}
