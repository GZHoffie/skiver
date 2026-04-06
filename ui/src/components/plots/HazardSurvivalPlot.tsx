import Plot from "react-plotly.js";
import { useCSVData } from "../../hooks/useCSVData";

interface HazardRow {
  t: number;
  hazard_ratio: number;
  "5th_percentile": number;
  "95th_percentile": number;
}

interface ErrorRateRow {
  lambda: number;
  beta: number;
}

interface Props {
  hazardPath: string;
  errorRatePath: string;
}

export function HazardSurvivalPlot({ hazardPath, errorRatePath }: Props) {
  const { data: hazard, loading: hLoading } = useCSVData<HazardRow>(hazardPath);
  const { data: rates, loading: rLoading } = useCSVData<ErrorRateRow>(errorRatePath);

  if (hLoading || rLoading) return <div className="plot-loading">Loading…</div>;
  if (!hazard || !rates || rates.length === 0) return null;

  const { lambda, beta } = rates[0];
  const t = hazard.map((r) => r.t);
  const hr = hazard.map((r) => r.hazard_ratio);
  const lo = hazard.map((r) => r["5th_percentile"]);
  const hi = hazard.map((r) => r["95th_percentile"]);

  const fittedHazard = t.map(
    (ti) => 1 - Math.exp(-lambda * (Math.pow(ti, beta) - Math.pow(ti - 1, beta)))
  );

  // Survival curve
  const tSurv = Array.from({ length: 100 }, (_, i) => i + 1);
  const survival = tSurv.map((ti) => Math.exp(-lambda * Math.pow(ti, beta)));

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const hazardTraces: any[] = [
    {
      type: "scatter", x: t, y: lo, mode: "lines",
      line: { color: "transparent" }, showlegend: false,
      name: "5th pct", xaxis: "x", yaxis: "y", legend: "legend",
    },
    {
      type: "scatter", x: t, y: hi, mode: "lines",
      fill: "tonexty", fillcolor: "rgba(112,128,144,0.25)",
      line: { color: "transparent" }, name: "5–95% percentile",
      xaxis: "x", yaxis: "y", legend: "legend",
    },
    {
      type: "scatter", x: t, y: hr, mode: "lines",
      line: { color: "slategray", width: 4 }, name: "Estimated hazard rate",
      xaxis: "x", yaxis: "y", legend: "legend",
    },
    {
      type: "scatter", x: t, y: fittedHazard, mode: "lines",
      line: { color: "steelblue", width: 4, dash: "dash" }, name: "Fitted hazard rate",
      xaxis: "x", yaxis: "y", legend: "legend",
    },
  ];

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const survivalTraces: any[] = [
    {
      type: "scatter", x: tSurv, y: survival, mode: "lines",
      line: { color: "steelblue", width: 4, dash: "dash" }, name: "Fitted survival rate",
      xaxis: "x2", yaxis: "y2", legend: "legend2",
    },
  ];

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const layout: any = {
    xaxis:  { title: { text: "<i>t</i>", font: { size: 14 } }, domain: [0, 0.44],  anchor: "y"  },
    yaxis:  { title: { text: "Hazard rate <i>h</i>(<i>t</i>)", font: { size: 14 } }, domain: [0, 0.75], anchor: "x"  },
    xaxis2: { title: { text: "<i>t</i>", font: { size: 14 } }, domain: [0.56, 1],  anchor: "y2" },
    yaxis2: { title: { text: "Survival rate <i>S</i>(<i>t</i>)", font: { size: 14 } }, domain: [0, 0.75], anchor: "x2" },
    legend: {
      x: 0, y: 0.75, xanchor: "left", yanchor: "bottom",
      font: { size: 14 },
    },
    legend2: {
      x: 0.56, y: 0.75, xanchor: "left", yanchor: "bottom",
      orientation: "h", font: { size: 14 },
    },
    font: { family: "DejaVu Sans", size: 13 },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 30, l: 90, r: 30, b: 70 },
    annotations: [
      // subplot titles
      { text: "Hazard rate",   xref: "paper", yref: "paper", x: 0.22, y: 1, showarrow: false, xanchor: "center", font: { size: 16 } },
      { text: "Probability of next base being erroneous given previous <i>t</i> bases are correct",
        xref: "paper", yref: "paper", x: 0.22, y: 0.95, showarrow: false, xanchor: "center",
        font: { size: 12, color: "#9099b0" } },
      { text: "Survival rate", xref: "paper", yref: "paper", x: 0.78, y: 1, showarrow: false, xanchor: "center", font: { size: 16 } },
      { text: "Probability that a <i>t</i>-mer is free of sequencing errors",
        xref: "paper", yref: "paper", x: 0.78, y: 0.95, showarrow: false, xanchor: "center",
        font: { size: 12, color: "#9099b0" } },
    ],
  };

  return (
    <Plot
      data={[...hazardTraces, ...survivalTraces]}
      layout={layout}
      useResizeHandler
      style={{ width: "100%", height: "550px" }}
    />
  );
}
