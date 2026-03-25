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
    xaxis:  { title: { text: "<i>t</i>" }, domain: [0, 0.44],  anchor: "y"  },
    yaxis:  {                       domain: [0, 1],     anchor: "x"  },
    xaxis2: { title: { text: "<i>t</i>" }, domain: [0.56, 1],  anchor: "y2" },
    yaxis2: {                       domain: [0, 1],     anchor: "x2" },
    legend: {
      x: 0, y: 1.02, xanchor: "left", yanchor: "bottom",
      orientation: "h",
    },
    legend2: {
      x: 0.56, y: 1.02, xanchor: "left", yanchor: "bottom",
      orientation: "h",
    },
    font: { family: "DejaVu Sans" },
    hoverlabel: { font: { family: "DejaVu Sans" } },
    margin: { t: 110, l: 80, r: 20, b: 50 },
    annotations: [
      // subplot titles
      { text: "Hazard rate",   xref: "paper", yref: "paper", x: 0.22, y: 1.26, showarrow: false, xanchor: "center", font: { size: 15 } },
      { text: "Error rate at each position",
        xref: "paper", yref: "paper", x: 0.22, y: 1.19, showarrow: false, xanchor: "center",
        font: { size: 11, color: "#9099b0" } },
      { text: "Survival rate", xref: "paper", yref: "paper", x: 0.78, y: 1.26, showarrow: false, xanchor: "center", font: { size: 15 } },
      { text: "Probability that a <i>t</i>-mer is free of sequencing errors",
        xref: "paper", yref: "paper", x: 0.78, y: 1.19, showarrow: false, xanchor: "center",
        font: { size: 11, color: "#9099b0" } },
      // y-axis titles (rotated annotations)
      { text: "Hazard rate <i>h</i>(<i>t</i>)", xref: "paper", yref: "paper", x: -0.07, y: 0.5, showarrow: false, textangle: -90, font: { size: 12 } },
      { text: "Survival rate <i>S</i>(<i>t</i>)",    xref: "paper", yref: "paper", x: 0.52,  y: 0.5, showarrow: false, textangle: -90, font: { size: 12 } },
    ],
  };

  return (
    <Plot
      data={[...hazardTraces, ...survivalTraces]}
      layout={layout}
      useResizeHandler
      style={{ width: "100%", height: "480px" }}
    />
  );
}
