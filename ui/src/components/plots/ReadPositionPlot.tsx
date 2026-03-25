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

  if (loading) return <div className="plot-loading">Loading…</div>;
  if (!rows) return null;

  const dfStart = rows
    .filter((r) => r.from_start)
    .sort((a, b) => a.index - b.index)
    .slice(0, 100);
  const dfEnd = rows
    .filter((r) => !r.from_start)
    .sort((a, b) => a.index - b.index)
    .slice(0, 100);

  const window = Math.max(3, Math.floor(dfStart.length * 0.05));

  function makeTraces(
    df: ReadPosRow[],
    xAxis: string,
    yAxis: string,
    yAxis2: string,
    label: string
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
        name: `${label} (raw)`,
        line: { color: "lightgray", width: 1 },
        xaxis: xAxis,
        yaxis: yAxis,
        showlegend: true,
      },
      {
        type: "scatter",
        x: idx,
        y: smoothed,
        mode: "lines",
        name: `${label} (smoothed)`,
        line: { color: "indianred", width: 2 },
        xaxis: xAxis,
        yaxis: yAxis,
        showlegend: true,
      },
      {
        type: "bar",
        x: idx,
        y: counts,
        name: `${label} count`,
        marker: { color: "steelblue", opacity: 0.5 },
        xaxis: xAxis,
        yaxis: yAxis2,
        showlegend: false,
      },
    ];
  }

  const traces = [
    ...makeTraces(dfStart, "x", "y", "y2", "From start"),
    ...makeTraces(dfEnd, "x3", "y3", "y4", "From end"),
  ];
  // Fix bar traces to use the correct x-axes for the bottom subplots
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  (traces[2] as any).xaxis = "x2";
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  (traces[5] as any).xaxis = "x4";

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    // top-left: from-start error rate
    xaxis:  { title: "",         domain: [0, 0.45],    anchor: "y",  showticklabels: false },
    yaxis:  {                    domain: [0.35, 1],    anchor: "x"  },
    xaxis2: { title: "Position", domain: [0, 0.45],    anchor: "y2" },
    yaxis2: {                    domain: [0, 0.3],     anchor: "x2" },
    xaxis3: { title: "",         domain: [0.55, 1],    anchor: "y3", showticklabels: false },
    yaxis3: {                    domain: [0.35, 1],    anchor: "x3" },
    xaxis4: { title: "Position", domain: [0.55, 1],    anchor: "y4" },
    yaxis4: {                    domain: [0, 0.3],     anchor: "x4" },
    annotations: [
      { text: "From start of read", xref: "paper", yref: "paper", x: 0.22,  y: 1.04, showarrow: false, font: { size: 13 } },
      { text: "From end of read",   xref: "paper", yref: "paper", x: 0.78,  y: 1.04, showarrow: false, font: { size: 13 } },
      { text: "Error rate",         xref: "paper", yref: "paper", x: -0.07, y: 0.67, showarrow: false, textangle: -90, font: { size: 12 } },
      { text: "Count",              xref: "paper", yref: "paper", x: -0.07, y: 0.15, showarrow: false, textangle: -90, font: { size: 12 } },
    ],
    margin: { t: 50, l: 80, r: 20, b: 50 },
  };

  return (
    <Plot
      data={traces}
      layout={layout}
      useResizeHandler
      style={{ width: "100%", height: "560px" }}
    />
  );
}
