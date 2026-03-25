import Plot from "react-plotly.js";
import { useCSVData } from "../../hooks/useCSVData";

interface SpectrumRow {
  operation: string;
  prev_base: string;
  next_base: string;
  total: number;
  forward: number;
}

interface ErrorRateRow {
  per_base_error_rate: number;
}

interface Props {
  spectrumPath: string;
  errorRatePath: string;
}

const BASES = ["A", "C", "G", "T", "-"];

function buildMatrix(
  rows: SpectrumRow[],
  col: "total" | "forward"
): (number | null)[][] {
  const matrix: (number | null)[][] = BASES.map(() => BASES.map(() => null));
  for (const row of rows) {
    const parts = row.operation.split(">");
    if (parts.length !== 2) continue;
    const [from, to] = parts;
    const fi = BASES.indexOf(from);
    const ti = BASES.indexOf(to);
    if (fi === -1 || ti === -1 || fi === ti) continue;
    matrix[fi][ti] = (matrix[fi][ti] ?? 0) + row[col];
  }
  return matrix;
}

function buildBarData(rows: SpectrumRow[], col: "total" | "forward") {
  let sub = 0, ins = 0, del = 0;
  for (const row of rows) {
    if (row.operation.includes("->")) ins += row[col];
    else if (row.operation.includes(">-")) del += row[col];
    else sub += row[col];
  }
  return { sub, ins, del };
}

export function ErrorSpectrumHeatmap({ spectrumPath, errorRatePath }: Props) {
  const { data: spectrum, loading: sLoading } = useCSVData<SpectrumRow>(spectrumPath);
  const { data: rates, loading: rLoading } = useCSVData<ErrorRateRow>(errorRatePath);

  if (sLoading || rLoading) return <div className="plot-loading">Loading…</div>;
  if (!spectrum || !rates || rates.length === 0) return null;

  const scale = rates[0].per_base_error_rate;

  const matTotal = buildMatrix(spectrum, "total").map((row) =>
    row.map((v) => (v === null ? null : v * scale))
  );
  const matFwd = buildMatrix(spectrum, "forward").map((row) =>
    row.map((v) => (v === null ? null : v * scale))
  );
  const barTotal = buildBarData(spectrum, "total");
  const barFwd = buildBarData(spectrum, "forward");

  const heatmapConfig = (
    z: (number | null)[][],
    axis: string,
    title: string
  ): Plotly.Data => ({
    type: "heatmap",
    z,
    x: BASES,
    y: BASES,
    colorscale: "slategray",
    showscale: true,
    xaxis: `x${axis}`,
    yaxis: `y${axis}`,
    name: title,
  });

  const barConfig = (
    b: { sub: number; ins: number; del: number },
    axis: string,
    title: string
  ): Plotly.Data => ({
    type: "bar",
    x: ["Substitution", "Insertion", "Deletion"],
    y: [b.sub * scale, b.ins * scale, b.del * scale],
    marker: { color: ["slategray", "slategray", "slategray"] },
    xaxis: `x${axis}`,
    yaxis: `y${axis}`,
    name: title,
  });

  const traces: Plotly.Data[] = [
    heatmapConfig(matTotal, "", "Total — heatmap"),
    barConfig(barTotal, "2", "Total — bar"),
    heatmapConfig(matFwd, "3", "Forward — heatmap"),
    barConfig(barFwd, "4", "Forward — bar"),
  ];

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    // top-left: heatmap total
    xaxis:  { title: "To",         domain: [0, 0.45],    anchor: "y"  },
    yaxis:  {                      domain: [0.52, 1],    anchor: "x"  },
    xaxis2: { title: "Error type", domain: [0.55, 1],    anchor: "y2" },
    yaxis2: {                      domain: [0.52, 1],    anchor: "x2" },
    xaxis3: { title: "To",         domain: [0, 0.45],    anchor: "y3" },
    yaxis3: {                      domain: [0, 0.46],    anchor: "x3" },
    xaxis4: { title: "Error type", domain: [0.55, 1],    anchor: "y4" },
    yaxis4: {                      domain: [0, 0.46],    anchor: "x4" },
    annotations: [
      { text: "Both strands",   xref: "paper", yref: "paper", x: 0.5,   y: 1.04, showarrow: false, font: { size: 14 } },
      { text: "Forward strand", xref: "paper", yref: "paper", x: 0.5,   y: 0.49, showarrow: false, font: { size: 14 } },
      { text: "From",           xref: "paper", yref: "paper", x: -0.07, y: 0.76, showarrow: false, textangle: -90, font: { size: 12 } },
      { text: "Rate",           xref: "paper", yref: "paper", x: 0.49,  y: 0.76, showarrow: false, textangle: -90, font: { size: 12 } },
      { text: "From",           xref: "paper", yref: "paper", x: -0.07, y: 0.23, showarrow: false, textangle: -90, font: { size: 12 } },
      { text: "Rate",           xref: "paper", yref: "paper", x: 0.49,  y: 0.23, showarrow: false, textangle: -90, font: { size: 12 } },
    ],
    margin: { t: 60, l: 80, r: 20, b: 50 },
  };

  return (
    <Plot
      data={traces}
      layout={layout}
      useResizeHandler
      style={{ width: "100%", height: "700px" }}
    />
  );
}
