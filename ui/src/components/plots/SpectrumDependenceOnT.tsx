import Plot from "react-plotly.js";
import { useCSVData } from "../../hooks/useCSVData";

// The CSV has columns: operation, prev_base, next_base, total, freq_at_t1, freq_at_t2, ...
interface DepTRow {
  operation: string;
  [key: string]: string | number;
}

interface Props {
  depTPath: string;
}

function classifyOp(op: string): "Substitution" | "Insertion" | "Deletion" {
  if (op.startsWith("->")) return "Insertion";
  if (op.endsWith(">-")) return "Deletion";
  return "Substitution";
}

export function SpectrumDependenceOnT({ depTPath }: Props) {
  const { data: rows, loading } = useCSVData<DepTRow>(depTPath);

  if (loading) return <div className="plot-loading">Loading…</div>;
  if (!rows || rows.length === 0) return null;

  // Find position columns
  const firstRow = rows[0];
  const posCols = Object.keys(firstRow).filter((k) => k.startsWith("freq_at_t"));
  const positions = posCols.map((c) => parseInt(c.replace("freq_at_t", ""), 10));

  // Aggregate by type
  const types = ["Substitution", "Insertion", "Deletion"] as const;
  const typeColors: Record<string, string> = {
    Substitution: "steelblue",
    Insertion: "indianred",
    Deletion: "seagreen",
  };

  const typeAggregated: Record<string, number[]> = {
    Substitution: new Array(positions.length).fill(0),
    Insertion: new Array(positions.length).fill(0),
    Deletion: new Array(positions.length).fill(0),
  };

  for (const row of rows) {
    const type = classifyOp(row.operation as string);
    for (let i = 0; i < posCols.length; i++) {
      typeAggregated[type][i] += (row[posCols[i]] as number) ?? 0;
    }
  }

  // Normalize to proportions per position
  const totals = positions.map((_, i) =>
    types.reduce((s, t) => s + typeAggregated[t][i], 0)
  );
  const typeProportions: Record<string, number[]> = {};
  for (const type of types) {
    typeProportions[type] = typeAggregated[type].map((v, i) =>
      totals[i] > 0 ? v / totals[i] : 0
    );
  }

  const leftTraces: Plotly.Data[] = types.map((type) => ({
    type: "scatter",
    x: positions,
    y: typeProportions[type],
    mode: "lines",
    name: type,
    line: { color: typeColors[type], width: 2 },
    xaxis: "x",
    yaxis: "y",
    legend: "legend",
  }));

  // Right panel: individual operations, normalised by total across all ops at each t
  const uniqueOps = [...new Set(rows.map((r) => r.operation as string))];

  const colTotals = positions.map((_, i) =>
    rows.reduce((s, r) => s + ((r[posCols[i]] as number) ?? 0), 0)
  );

  // Colour palettes: substitutions → steelblue shades,
  //                 insertions    → indianred shades,
  //                 deletions     → seagreen shades
  const SUB_BLUES  = ["#c6dff5","#93c4ea","#5ea8de","#3b8ccc",
                      "#1e70b8","#1058a0","#084288","#042e6e",
                      "#021c50","#01112e","#000a1c","#00040e"];
  const INS_REDS   = ["#fcc4c4","#f08080","#cd5c5c","#9e2222"];
  const DEL_GREENS = ["#b8e8cc","#70c898","#2e8b57","#125e35"];

  const subOps = uniqueOps.filter((op) => !op.startsWith("->") && !op.endsWith(">-"));
  const insOps = uniqueOps.filter((op) => op.startsWith("->"));
  const delOps = uniqueOps.filter((op) => op.endsWith(">-"));

  const opColorMap: Record<string, string> = {};
  subOps.forEach((op, i) => { opColorMap[op] = SUB_BLUES[i % SUB_BLUES.length]; });
  insOps.forEach((op, i) => { opColorMap[op] = INS_REDS[i % INS_REDS.length]; });
  delOps.forEach((op, i) => { opColorMap[op] = DEL_GREENS[i % DEL_GREENS.length]; });

  const rightTraces: Plotly.Data[] = uniqueOps.map((op) => {
    const opRows = rows.filter((r) => r.operation === op);
    const vals = positions.map((_, i) => {
      const raw = opRows.reduce((s, r) => s + ((r[posCols[i]] as number) ?? 0), 0);
      return colTotals[i] > 0 ? raw / colTotals[i] : 0;
    });
    return {
      type: "scatter",
      x: positions,
      y: vals,
      mode: "lines",
      name: op,
      line: { color: opColorMap[op] ?? "#888", width: 1.5 },
      xaxis: "x2",
      yaxis: "y2",
      legend: "legend2",
    };
  });

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    xaxis:  { title: { text: "<i>t</i>", font: { size: 12 } }, domain: [0, 0.45], anchor: "y"  },
    yaxis:  {                                domain: [0, 1],    anchor: "x"  },
    xaxis2: { title: { text: "<i>t</i>", font: { size: 12 } }, domain: [0.55, 1], anchor: "y2" },
    yaxis2: {                                domain: [0, 1],    anchor: "x2" },
    legend: {
      x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", orientation: "h",
    },
    legend2: {
      x: 1.02, y: 1.0, xanchor: "left", yanchor: "top",
    },
    annotations: [
      { text: "Error type composition", xref: "paper", yref: "paper", x: 0.22, y: 1.18, showarrow: false, xanchor: "center", font: { size: 15 } },
      { text: "Individual operations",  xref: "paper", yref: "paper", x: 0.78, y: 1.18, showarrow: false, xanchor: "center", font: { size: 15 } },
      { text: "Proportion",             xref: "paper", yref: "paper", x: -0.07, y: 0.5, showarrow: false, textangle: -90, font: { size: 12 } },
      { text: "Proportion",             xref: "paper", yref: "paper", x: 0.49,  y: 0.5, showarrow: false, textangle: -90, font: { size: 12 } },
    ],
    font: { family: "DejaVu Sans" },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 90, l: 80, r: 160, b: 50 },
  };

  return (
    <Plot
      data={[...leftTraces, ...rightTraces]}
      layout={layout}
      useResizeHandler
      style={{ width: "100%", height: "460px" }}
    />
  );
}
