import Plot from "react-plotly.js";
import { useCSVData } from "../../hooks/useCSVData";

interface SpectrumRow {
  operation: string;
  prev_base: string;
  next_base: string;
  total: number;
}

interface Props {
  spectrumPath: string;
}

const DNA = ["A", "C", "G", "T"];
const CANONICAL_MUTS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"];
const MUT_COLORS: Record<string, string> = {
  "C>A": "#03bcee",
  "C>G": "#010101",
  "C>T": "#e32926",
  "T>A": "#cbcacb",
  "T>C": "#a1ce63",
  "T>G": "#ecc6c5",
};

export function SBS96Plot({ spectrumPath }: Props) {
  const { data: spectrum, loading } = useCSVData<SpectrumRow>(spectrumPath);

  if (loading) return <div className="plot-loading">Loading…</div>;
  if (!spectrum) return null;

  // Filter to substitutions only
  const subs = spectrum.filter(
    (r) => !r.operation.includes("->") && !r.operation.includes(">-")
  );

  // Detect if we have all 12 substitution types or only canonical (C>*, T>*)
  const hasNonCanonical = subs.some(
    (r) => !CANONICAL_MUTS.includes(r.operation)
  );
  const mutations = hasNonCanonical
    ? [...new Set(subs.map((r) => r.operation))].sort()
    : CANONICAL_MUTS;

  // Build lookup map
  const lookup = new Map<string, number>();
  for (const row of subs) {
    const key = `${row.prev_base}|${row.operation}|${row.next_base}`;
    lookup.set(key, (lookup.get(key) ?? 0) + row.total);
  }

  const traces: Plotly.Data[] = mutations.map((mut) => {
    const xLabels: string[] = [];
    const yValues: number[] = [];
    for (const prev of DNA) {
      for (const next of DNA) {
        xLabels.push(`${prev}[${mut}]${next}`);
        yValues.push(lookup.get(`${prev}|${mut}|${next}`) ?? 0);
      }
    }
    return {
      type: "bar",
      name: mut,
      x: xLabels,
      y: yValues,
      marker: { color: MUT_COLORS[mut] ?? "gray" },
    };
  });

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    barmode: "group",
    xaxis: {
      tickangle: -90,
      tickfont: { size: 9 },
      title: { text: "Trinucleotide context", y: 0.65 },
    },
    yaxis: {},
    legend: { orientation: "h", y: 1.1 },
    font: { family: "DejaVu Sans" },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 60, b: 120, l: 80, r: 20 },
    annotations: [
      { text: "Count", xref: "paper", yref: "paper", x: -0.06, y: 0.5, showarrow: false, textangle: -90, font: { size: 12 } },
    ],
  };

  return (
    <Plot
      data={traces}
      layout={layout}
      useResizeHandler
      style={{ width: "100%", height: "480px" }}
    />
  );
}
