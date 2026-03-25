import { useState, useMemo } from "react";
import Plot from "react-plotly.js";
import { useCSVData } from "../../hooks/useCSVData";

interface ReadPosRow {
  index: number;
  from_start: boolean;
  num_correct: number;
  num_error: number;
  error_rate: number;
}

interface Props {
  readPositionPath: string;
}

function movingAverage(arr: number[], window: number): number[] {
  const half = Math.floor(window / 2);
  return arr.map((_, i) => {
    const lo = Math.max(0, i - half);
    const hi = Math.min(arr.length - 1, i + half);
    const slice = arr.slice(lo, hi + 1);
    return slice.reduce((a, b) => a + b, 0) / slice.length;
  });
}

export function ReadPositionPlot({ readPositionPath }: Props) {
  const { data: rows, loading } = useCSVData<ReadPosRow>(readPositionPath);

  const allStart = useMemo(
    () => (rows ?? []).filter((r) => r.from_start).sort((a, b) => a.index - b.index),
    [rows]
  );
  const allEnd = useMemo(
    () => (rows ?? []).filter((r) => !r.from_start).sort((a, b) => a.index - b.index),
    [rows]
  );

  const maxPos = useMemo(
    () => Math.max(
      allStart.length > 0 ? allStart[allStart.length - 1].index : 0,
      allEnd.length > 0 ? allEnd[allEnd.length - 1].index : 0,
      1
    ),
    [allStart, allEnd]
  );

  const [xLimit, setXLimit] = useState<number>(100);
  const effectiveLimit = xLimit;

  if (loading) return <div className="plot-loading">Loading…</div>;
  if (!rows) return null;

  const dfStart = allStart;
  const dfEnd = allEnd;

  const window = 5; // smoothing window for error rate

  function makeTraces(
    df: ReadPosRow[],
    xAxis: string,
    yAxis: string,
    yAxis2: string,
    showLegend: boolean
  ): Plotly.Data[] {
    const idx = df.map((r) => r.index);
    const rate = df.map((r) => r.error_rate);
    const counts = df.map((r) => r.num_correct + r.num_error);
    const smoothed = movingAverage(rate, window);

    return [
      {
        type: "scatter",
        x: idx,
        y: rate,
        mode: "lines",
        name: "Error rate",
        line: { color: "lightgray", width: 1 },
        xaxis: xAxis,
        yaxis: yAxis,
        showlegend: showLegend,
      },
      {
        type: "scatter",
        x: idx,
        y: smoothed,
        mode: "lines",
        name: "Error rate (smoothed)",
        line: { color: "indianred", width: 2 },
        xaxis: xAxis,
        yaxis: yAxis,
        showlegend: showLegend,
      },
      {
        type: "bar",
        x: idx,
        y: counts,
        name: "Count",
        marker: { color: "steelblue", opacity: 0.5 },
        xaxis: xAxis,
        yaxis: yAxis2,
        showlegend: false,
      },
    ];
  }

  const traces = [
    ...makeTraces(dfStart, "x", "y", "y2", true),
    ...makeTraces(dfEnd, "x3", "y3", "y4", false),
  ];
  // Fix bar traces to use the correct x-axes for the bottom subplots
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  (traces[2] as any).xaxis = "x2";
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  (traces[5] as any).xaxis = "x4";

  const xRange: [number, number] = [0, effectiveLimit];

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    // top-left: from-start error rate
    xaxis:  { title: "",         domain: [0, 0.45],    anchor: "y",  showticklabels: false, range: xRange },
    yaxis:  {                    domain: [0.35, 1],    anchor: "x"  },
    xaxis2: { title: { text: "Position" }, domain: [0, 0.45],    anchor: "y2", range: xRange },
    yaxis2: {                              domain: [0, 0.3],     anchor: "x2" },
    xaxis3: { title: "",                   domain: [0.55, 1],    anchor: "y3", showticklabels: false, range: xRange },
    yaxis3: {                              domain: [0.35, 1],    anchor: "x3" },
    xaxis4: { title: { text: "Position" }, domain: [0.55, 1],    anchor: "y4", range: xRange },
    yaxis4: {                    domain: [0, 0.3],     anchor: "x4" },
    annotations: [
      { text: "From start of read", xref: "paper", yref: "paper", x: 0.15,  y: 1.07, showarrow: false, font: { size: 15 } },
      { text: "From end of read",   xref: "paper", yref: "paper", x: 0.85,  y: 1.07, showarrow: false, font: { size: 15 } },
      { text: "Error rate",         xref: "paper", yref: "paper", x: -0.07, y: 0.75, showarrow: false, textangle: -90, font: { size: 12 } },
      { text: "# bases used for estimation", xref: "paper", yref: "paper", x: -0.07, y: -0.05, showarrow: false, textangle: -90, font: { size: 12 } },
    ],
    font: { family: "DejaVu Sans" },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 50, l: 80, r: 20, b: 50 },
  };

  return (
    <div style={{ width: "100%" }}>
      <div style={{ display: "flex", alignItems: "center", gap: "12px", padding: "0 80px 20px 80px" }}>
        <label style={{ whiteSpace: "nowrap", fontSize: "13px" }}>
          X limit: {effectiveLimit}
        </label>
        <input
          type="range"
          min={1}
          max={Math.min(maxPos, 1000)}
          value={effectiveLimit}
          onChange={(e) => setXLimit(Number(e.target.value))}
          style={{ flex: 1 }}
        />
      </div>
      <Plot
        data={traces}
        layout={layout}
        useResizeHandler
        style={{ width: "100%", height: "560px" }}
      />
    </div>
  );
}
