import React from "react";
import { useCSVData } from "../hooks/useCSVData";

interface ErrorRateRow {
  per_base_error_rate: number;
  "per_base_error_rate_5-95th_percentile": string;
  effective_error_rate: number;
  "effective_error_rate_5-95th_percentile": string;
  lambda: number;
  "lambda_5-95th_percentile": string;
  beta: number;
  "beta_5-95th_percentile": string;
  key_median_coverage: number;
  "key_coverage_5-95th_percentile": string;
  true_median_coverage: number;
  "true_coverage_5-95th_percentile": string;
}

interface Props {
  errorRatePath: string;
}

function fmtPct(v: number, digits = 4): string {
  if (v === undefined || v === null || isNaN(v)) return "—";
  return (v * 100).toPrecision(digits) + "%";
}

function fmtCoverage(v: number): string {
  if (v === undefined || v === null || isNaN(v)) return "—";
  return Math.round(v) + "x";
}

function transformCi(ci: string, transform: (v: number) => string): string {
  if (!ci) return ci;
  const [lo, hi] = ci.split("~");
  const loN = parseFloat(lo);
  const hiN = parseFloat(hi);
  if (isNaN(loN) || isNaN(hiN)) return ci;
  return `${transform(loN)}~${transform(hiN)}`;
}

function Box({
  label,
  value,
  ci,
  description,
}: {
  label: string;
  value: string;
  ci: string;
  description?: React.ReactNode;
}) {
  const [lo, hi] = ci ? ci.split("~") : ["", ""];
  return (
    <div className="summary-box">
      <div className="summary-box-label">
        {label}
        {description && (
          <span className="info-icon">
            i
            <span className="info-tooltip">{description}</span>
          </span>
        )}
      </div>
      <div className="summary-box-value">{value}</div>
      {lo && hi && (
        <div className="summary-box-ci">
          [{lo} – {hi}]
        </div>
      )}
    </div>
  );
}

export function SummaryBoxes({ errorRatePath }: Props) {
  const { data, loading } = useCSVData<ErrorRateRow>(errorRatePath);

  if (loading) return <div className="summary-boxes-loading">Loading summary…</div>;
  if (!data || data.length === 0) return null;

  const r = data[0];

  return (
    <div className="summary-boxes">
      <Box
        label="Per-base error rate"
        value={fmtPct(r.per_base_error_rate)}
        ci={transformCi(r["per_base_error_rate_5-95th_percentile"], (v) => fmtPct(v))}
        description={<>Probability that a random base in the read is erroneous. Note: This value is calculated by extrapolating from the hazard rate and can be sensitive to the parameters.</>}
      />
      <Box
        label="Effective error rate"
        value={fmtPct(r.effective_error_rate)}
        ci={transformCi(r["effective_error_rate_5-95th_percentile"], (v) => fmtPct(v))}
        description={<>Similar to the gap-compressed error rate. This value estimates the probability of a random base is erroneous after excluding the regions where errors are clustered together.</>}
      />
      <Box
        label="Median coverage"
        value={fmtCoverage(r.true_median_coverage)}
        ci={transformCi(r["true_coverage_5-95th_percentile"], fmtCoverage)}
        description={<>Median coverage of the keys in the (<i>k</i>, <i>v</i>)-mer sketch, corrected by the error rate estimation.</>}
      />
    </div>
  );
}
