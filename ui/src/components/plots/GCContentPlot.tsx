import { useMemo, useState } from "react";
import Plot from "react-plotly.js";
import { useCSVData } from "../../hooks/useCSVData";

interface GcContentRow {
  gc_content_min: number;
  gc_content_max_exclusive: number;
  per_base_error_rate?: number;
  per_base_error_rate_median?: number;
  "per_base_error_rate_5-95th_percentile"?: string;
  num_correct: number;
  num_error: number;
}

interface Props {
  gcContentPath: string;
}

const MIN_BASE_OPTIONS = [0, 100, 500, 1000, 5000, 10000];

function errorRate(row: GcContentRow): number {
  return row.per_base_error_rate ?? row.per_base_error_rate_median ?? 0;
}

function parseCi(row: GcContentRow): { lo: number; hi: number } {
  const [lo, hi] = (row["per_base_error_rate_5-95th_percentile"] ?? "").split("~");
  return { lo: Number(lo), hi: Number(hi) };
}

export function GCContentPlot({ gcContentPath }: Props) {
  const { data: rows, loading } = useCSVData<GcContentRow>(gcContentPath);
  const [logScale, setLogScale] = useState(false);
  const [minBases, setMinBases] = useState(5000);

  const prepared = useMemo(() => {
    return (rows ?? [])
      .map((row) => {
        const gcMin = Number(row.gc_content_min);
        const gcMax = Number(row.gc_content_max_exclusive);
        const numCorrect = Number(row.num_correct);
        const numError = Number(row.num_error);
        return {
          ...row,
          gcMid: (gcMin + gcMax - 1) / 2,
          binWidth: gcMax - gcMin,
          numTotal: numCorrect + numError,
          errorRate: errorRate(row),
          ci: parseCi(row),
        };
      })
      .filter((row) => Number.isFinite(row.gcMid) && row.numTotal > 0)
      .sort((a, b) => a.gcMid - b.gcMid);
  }, [rows]);

  if (loading) return <div className="plot-loading">Loading…</div>;
  if (!rows) return null;
  if (prepared.length === 0) return <div className="plot-loading">No GC content rows available.</div>;

  const lineRows = prepared.filter((row) => row.numTotal >= minBases);
  const hasCi = lineRows.every((row) => Number.isFinite(row.ci.lo) && Number.isFinite(row.ci.hi));

  const xLine = lineRows.map((row) => row.gcMid);
  const yLine = lineRows.map((row) => row.errorRate);
  const ciLower = lineRows.map((row) => row.ci.lo);
  const ciUpper = lineRows.map((row) => row.ci.hi);
  const xHist = prepared.map((row) => row.gcMid);
  const counts = prepared.map((row) => row.numTotal);
  const binWidth = prepared[0]?.binWidth ?? 1;

  const traces: Plotly.Data[] = [
    ...(hasCi && lineRows.length > 0
      ? [
          {
            type: "scatter" as const,
            x: xLine,
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
            x: xLine,
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
      x: xLine,
      y: yLine,
      mode: "lines+markers",
      name: "Empirical error rate",
      line: { color: "slategray", width: 4 },
      marker: { size: 9 },
      xaxis: "x",
      yaxis: "y",
    },
    {
      type: "bar",
      x: xHist,
      y: counts,
      width: binWidth * 0.9,
      name: "Number of bases",
      marker: { color: "slategray", opacity: 0.6 },
      xaxis: "x2",
      yaxis: "y2",
    },
  ];

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    xaxis: {
      title: "",
      domain: [0, 1],
      anchor: "y",
      showticklabels: false,
      matches: "x2",
    },
    yaxis: {
      title: { text: logScale ? "Error rate (log scale)" : "Error rate", font: { size: 14 } },
      domain: [0.32, 1],
      anchor: "x",
      type: logScale ? "log" : "linear",
    },
    xaxis2: {
      title: { text: "GC content (%)", font: { size: 14 } },
      domain: [0, 1],
      anchor: "y2",
    },
    yaxis2: {
      title: { text: "# bases", font: { size: 14 } },
      domain: [0, 0.28],
      anchor: "x2",
    },
    legend: { x: 1, xanchor: "right", y: 1, font: { size: 14 } },
    font: { family: "DejaVu Sans", size: 13 },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 30, l: 80, r: 20, b: 50 },
  };

  return (
    <div style={{ width: "100%" }}>
      <div style={{ display: "flex", justifyContent: "flex-end", gap: "12px", padding: "0 20px 8px 0" }}>
        <label style={{ display: "flex", alignItems: "center", gap: "6px", fontSize: "13px" }}>
          Minimum bases
          <select
            value={minBases}
            onChange={(e) => setMinBases(Number(e.target.value))}
            style={{ fontSize: "13px" }}
          >
            {MIN_BASE_OPTIONS.map((value) => (
              <option key={value} value={value}>
                {value}
              </option>
            ))}
          </select>
        </label>
        <button
          onClick={() => setLogScale((v) => !v)}
          style={{ fontSize: "13px", padding: "3px 10px", cursor: "pointer" }}
        >
          {logScale ? "Switch to Linear scale" : "Switch to Log scale"}
        </button>
      </div>
      {lineRows.length === 0 && (
        <div className="plot-loading">No GC bins above the minimum-bases threshold.</div>
      )}
      <Plot
        data={traces}
        layout={layout}
        useResizeHandler
        style={{ width: "100%", height: "650px" }}
      />
    </div>
  );
}
