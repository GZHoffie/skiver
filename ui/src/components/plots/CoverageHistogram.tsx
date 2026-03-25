import { useEffect, useState } from "react";
import Plot from "react-plotly.js";
import { invoke } from "@tauri-apps/api/core";
import Papa from "papaparse";
import { useCSVData } from "../../hooks/useCSVData";

interface ErrorRateRow {
  lambda: number;
  beta: number;
}

interface KvmerRow {
  total_count: number;
  passes_filter: string | boolean;
  key: string;
}

interface Props {
  kvmerPath: string;
  errorRatePath: string;
}

const PERCENTILE = 0.999;

function percentileThreshold(arr: number[], p: number): number {
  const sorted = [...arr].sort((a, b) => a - b);
  return sorted[Math.floor(sorted.length * p)] ?? Infinity;
}

export function CoverageHistogram({ kvmerPath, errorRatePath }: Props) {
  const { data: rates, loading: rLoading } = useCSVData<ErrorRateRow>(errorRatePath);
  const [allCov, setAllCov] = useState<number[] | null>(null);
  const [passingCov, setPassingCov] = useState<number[] | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    if (!kvmerPath || !rates || rates.length === 0) return;

    const { lambda, beta } = rates[0];

    // Estimate k from first few characters of any key
    // We'll determine correction factor after parsing first row
    let correctionFactor = 1;
    let kDetermined = false;

    setLoading(true);
    const all: number[] = [];
    const passing: number[] = [];

    invoke<string>("read_csv_file", { path: kvmerPath })
      .then((csv) => {
        Papa.parse<KvmerRow>(csv, {
          header: true,
          dynamicTyping: true,
          skipEmptyLines: true,
          step(result) {
            const row = result.data;
            if (!kDetermined && row.key) {
              const k = (row.key as string).length;
              correctionFactor = Math.exp(lambda * Math.pow(k, beta));
              kDetermined = true;
            }
            const corrected = (row.total_count as number) * correctionFactor;
            all.push(corrected);
            const passes =
              row.passes_filter === true ||
              row.passes_filter === "true" ||
              row.passes_filter === "True";
            if (passes) passing.push(corrected);
          },
          complete() {
            setAllCov(all);
            setPassingCov(passing);
            setLoading(false);
          },
          error() {
            setLoading(false);
          },
        });
      })
      .catch(() => setLoading(false));
  }, [kvmerPath, rates]);

  if (rLoading || loading) return <div className="plot-loading">Loading…</div>;
  if (!allCov) return null;

  const cap = percentileThreshold(allCov, PERCENTILE);
  const allFiltered = allCov.filter((v) => v <= cap);
  const passingFiltered = (passingCov ?? []).filter((v) => v <= cap);

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const traces: any[] = [
    {
      type: "histogram",
      x: allFiltered,
      name: "All keys",
      opacity: 0.5,
      marker: { color: "slategray" },
      nbinsx: 100,
    },
    {
      type: "histogram",
      x: passingFiltered,
      name: "Passing filter",
      opacity: 0.7,
      marker: { color: "steelblue" },
      nbinsx: 100,
    },
  ];

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    barmode: "overlay",
    xaxis: { title: { text: "Estimated true coverage" } },
    yaxis: { type: "log" },
    font: { family: "DejaVu Sans" },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 30, l: 80, r: 20, b: 50 },
    annotations: [
      { text: "Count", xref: "paper", yref: "paper", x: -0.06, y: 0.5, showarrow: false, textangle: -90, font: { size: 14 } },
    ],
  };

  return (
    <Plot
      data={traces}
      layout={layout}
      useResizeHandler
      style={{ width: "100%", height: "420px" }}
    />
  );
}
