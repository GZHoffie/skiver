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
  }));

  // Right panel: individual operations
  const uniqueOps = [...new Set(rows.map((r) => r.operation as string))];
  const opColors = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"];

  const rightTraces: Plotly.Data[] = uniqueOps.map((op, idx) => {
    const opRows = rows.filter((r) => r.operation === op);
    const vals = positions.map((_, i) =>
      opRows.reduce((s, r) => s + ((r[posCols[i]] as number) ?? 0), 0)
    );
    return {
      type: "scatter",
      x: positions,
      y: vals,
      mode: "lines",
      name: op,
      line: { color: opColors[idx % opColors.length], width: 1 },
      xaxis: "x2",
      yaxis: "y2",
    };
  });

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    xaxis:  { title: "Position t", domain: [0, 0.45],  anchor: "y"  },
    yaxis:  {                      domain: [0, 1],     anchor: "x"  },
    xaxis2: { title: "Position t", domain: [0.55, 1],  anchor: "y2" },
    yaxis2: {                      domain: [0, 1],     anchor: "x2" },
    annotations: [
      { text: "Error type composition", xref: "paper", yref: "paper", x: 0.22,  y: 1.05, showarrow: false, font: { size: 13 } },
      { text: "Individual operations",  xref: "paper", yref: "paper", x: 0.78,  y: 1.05, showarrow: false, font: { size: 13 } },
      { text: "Proportion",             xref: "paper", yref: "paper", x: -0.07, y: 0.5,  showarrow: false, textangle: -90, font: { size: 12 } },
      { text: "Frequency",              xref: "paper", yref: "paper", x: 0.49,  y: 0.5,  showarrow: false, textangle: -90, font: { size: 12 } },
    ],
    margin: { t: 50, l: 80, r: 20, b: 50 },
  };

  return (
    <Plot
      data={[...leftTraces, ...rightTraces]}
      layout={layout}
      useResizeHandler
      style={{ width: "100%", height: "420px" }}
    />
  );
}
