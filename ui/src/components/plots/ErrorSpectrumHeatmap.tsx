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

  // normalize matrix so that they sum to 100.
  const total = matrix.reduce((sum: number, row) => sum + row.reduce((rSum: number, v) => rSum + (v ?? 0), 0), 0);
  if (total > 0) {
    for (let i = 0; i < matrix.length; i++) {
      for (let j = 0; j < matrix[i].length; j++) {
        if (matrix[i][j] !== null) {
          matrix[i][j] = matrix[i][j]! / total * 100;
        }
      }
    }
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
  // normalize these counts so that they sum to 100.
  const total = sub + ins + del;
  if (total > 0) {
    sub = sub / total * 100;
    ins = ins / total * 100;
    del = del / total * 100;
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

  // Row domains with a generous gap between them
  const topY:    [number, number] = [0.58, 1.0];
  const bottomY: [number, number] = [0,    0.40];
  // Heatmap column; colorbar sits in the gap; bar chart fills the rest
  const heatX: [number, number] = [0, 0.34];
  const barX:  [number, number] = [0.54, 1];
  const colorbarX = 0.36;

  const COLORSCALE: [number, string][] = [[0, "whitesmoke"], [1, "slategray"]];

  const heatmapConfig = (
    z: (number | null)[][],
    axis: string,
    title: string,
    colorbarY: number,
    colorbarLen: number,
  ): Plotly.Data => ({
    type: "heatmap",
    z,
    x: BASES,
    y: BASES,
    colorscale: COLORSCALE,
    zmin: 0,
    xgap: 6,
    ygap: 6,
    showscale: true,
    colorbar: {
      x: colorbarX,
      xanchor: "left",
      y: colorbarY,
      yanchor: "middle",
      len: colorbarLen,
      thickness: 12,
      tickfont: { size: 13 },
      outlinewidth: 0,
    },
    text: z.map((row) => row.map((v) => (v === null ? "" : v.toFixed(3)))) as unknown as string[],
    texttemplate: "%{text}",
    textfont: { size: 13 },
    xaxis: `x${axis}`,
    yaxis: `y${axis}`,
    name: title,
  });

  const barConfig = (
    b: { sub: number; ins: number; del: number },
    axis: string,
    title: string
  ): Plotly.Data => {
    const yVals = [b.sub * scale, b.ins * scale, b.del * scale];
    return {
      type: "bar",
      x: ["Substitution", "Insertion", "Deletion"],
      y: yVals,
      text: yVals.map((v) => v.toFixed(3)) as unknown as string[],
      textposition: "outside" as const,
      textfont: { size: 12 },
      width: 0.5,
      marker: { color: "slategray", opacity: 0.6 },
      xaxis: `x${axis}`,
      yaxis: `y${axis}`,
      name: title,
      showlegend: false,
    };
  };

  const traces: Plotly.Data[] = [
    heatmapConfig(matTotal, "",  "Both strands — heatmap",  0.78,  0.44),
    barConfig(barTotal,     "2", "Both strands — bar"),
    heatmapConfig(matFwd,   "3", "Forward strand — heatmap",  0.215,  0.43),
    barConfig(barFwd,       "4", "Forward strand — bar"),
  ];

  const gridOff = { showgrid: false, zeroline: false };

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    xaxis:  { title: { text: "Observed base", font: { size: 14 } }, domain: heatX, anchor: "y",  ...gridOff },
    yaxis:  { title: { text: "Original base", font: { size: 14 } }, domain: topY,    anchor: "x",  autorange: "reversed", ...gridOff },
    xaxis2: { title: { text: "Error type", font: { size: 14 } },   domain: barX,    anchor: "y2" },
    yaxis2: { title: { text: "Error rate (%)", font: { size: 14 } }, domain: topY,   anchor: "x2" },
    xaxis3: { title: { text: "Observed base", font: { size: 14 } }, domain: heatX,  anchor: "y3", ...gridOff },
    yaxis3: { title: { text: "Original base", font: { size: 14 } }, domain: bottomY, anchor: "x3", autorange: "reversed", ...gridOff },
    xaxis4: { title: { text: "Error type", font: { size: 14 } },   domain: barX,    anchor: "y4" },
    yaxis4: { title: { text: "Error rate (%)", font: { size: 14 } }, domain: bottomY, anchor: "x4" },
    annotations: [
      { text: "Estimation of error rates (%) using both forward and reverse strands of the reads",   xref: "paper", yref: "paper", x: 0.5,   y: 1.04, showarrow: false, xanchor: "center", font: { size: 16 } },
      { text: "Estimation of error rates (%) using only forward strand of the reads", xref: "paper", yref: "paper", x: 0.5,   y: 0.45, showarrow: false, xanchor: "center", font: { size: 16 } },
    ],
    font: { family: "DejaVu Sans", size: 13 },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 60, l: 80, r: 20, b: 50 },
    showlegend: false,
  };

  return (
    <Plot
      data={traces}
      layout={layout}
      useResizeHandler
      style={{ width: "100%", height: "850px" }}
    />
  );
}
