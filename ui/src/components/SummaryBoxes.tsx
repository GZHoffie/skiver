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

function fmt(v: number, digits = 4): string {
  if (v === undefined || v === null || isNaN(v)) return "—";
  return v.toPrecision(digits);
}

function Box({
  label,
  value,
  ci,
}: {
  label: string;
  value: string;
  ci: string;
}) {
  const [lo, hi] = ci ? ci.split("~") : ["", ""];
  return (
    <div className="summary-box">
      <div className="summary-box-label">{label}</div>
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
        value={fmt(r.per_base_error_rate)}
        ci={r["per_base_error_rate_5-95th_percentile"]}
      />
      <Box
        label="Effective error rate"
        value={fmt(r.effective_error_rate)}
        ci={r["effective_error_rate_5-95th_percentile"]}
      />
      <Box
        label="λ (lambda)"
        value={fmt(r.lambda)}
        ci={r["lambda_5-95th_percentile"]}
      />
      <Box
        label="β (beta)"
        value={fmt(r.beta)}
        ci={r["beta_5-95th_percentile"]}
      />
      <Box
        label="Median k-mer coverage"
        value={fmt(r.key_median_coverage, 5)}
        ci={r["key_coverage_5-95th_percentile"]}
      />
      <Box
        label="True median coverage"
        value={fmt(r.true_median_coverage, 5)}
        ci={r["true_coverage_5-95th_percentile"]}
      />
    </div>
  );
}
