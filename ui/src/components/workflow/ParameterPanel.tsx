import { ReactNode } from "react";
import { SketchParams, AnalyzeParams } from "../../types";

interface Props {
  sketch: SketchParams;
  analyze: AnalyzeParams;
  onSketchChange: (p: SketchParams) => void;
  onAnalyzeChange: (p: AnalyzeParams) => void;
  disabled: boolean;
}

function InfoIcon({ description }: { description: ReactNode }) {
  return (
    <span className="info-icon">
      i
      <span className="info-tooltip">{description}</span>
    </span>
  );
}

function FieldLabel({ label, description }: { label: ReactNode; description?: ReactNode }) {
  return (
    <label>
      {label}
      {description && <InfoIcon description={description} />}
    </label>
  );
}

function NumField({
  label,
  description,
  value,
  onChange,
  min,
  max,
  step,
  disabled,
}: {
  label: ReactNode;
  description?: ReactNode;
  value: number;
  onChange: (v: number) => void;
  min?: number;
  max?: number;
  step?: number;
  disabled: boolean;
}) {
  return (
    <div className="field-row">
      <FieldLabel label={label} description={description} />
      <input
        type="number"
        value={value}
        min={min}
        max={max}
        step={step ?? 1}
        disabled={disabled}
        onChange={(e) => onChange(Number(e.target.value))}
      />
    </div>
  );
}

function NumFieldOpt({
  label,
  description,
  value,
  onChange,
  min,
  step,
  placeholder,
  disabled,
}: {
  label: ReactNode;
  description?: ReactNode;
  value: number | null;
  onChange: (v: number | null) => void;
  min?: number;
  step?: number;
  placeholder?: string;
  disabled: boolean;
}) {
  return (
    <div className="field-row">
      <FieldLabel label={label} description={description} />
      <input
        type="number"
        value={value ?? ""}
        placeholder={placeholder ?? "auto"}
        min={min}
        step={step ?? 1}
        disabled={disabled}
        onChange={(e) => onChange(e.target.value === "" ? null : Number(e.target.value))}
      />
    </div>
  );
}

function CheckField({
  label,
  description,
  value,
  onChange,
  disabled,
}: {
  label: ReactNode;
  description?: string;
  value: boolean;
  onChange: (v: boolean) => void;
  disabled: boolean;
}) {
  return (
    <div className="field-row">
      <FieldLabel label={label} description={description} />
      <input
        type="checkbox"
        checked={value}
        disabled={disabled}
        onChange={(e) => onChange(e.target.checked)}
      />
    </div>
  );
}

export function ParameterPanel({
  sketch,
  analyze,
  onSketchChange,
  onAnalyzeChange,
  disabled,
}: Props) {
  function setSketch<K extends keyof SketchParams>(key: K, val: SketchParams[K]) {
    onSketchChange({ ...sketch, [key]: val });
  }
  function setAnalyze<K extends keyof AnalyzeParams>(key: K, val: AnalyzeParams[K]) {
    onAnalyzeChange({ ...analyze, [key]: val });
  }

  return (
    <div className="section">
      <h3>Parameters</h3>
      <div className="param-grid">
        <div>
          <h4>Sketch / Analyze</h4>
          <NumField
            label={<><i>k</i> (key length)</>}
            description={<>Length of the key in the (<i>k</i>, <i>v</i>)-mer sketch. Must be ≤ 31.</>}
            value={sketch.k}
            onChange={(v) => { setSketch("k", v); setAnalyze("k", v); }}
            min={1} max={31} disabled={disabled}
          />
          <NumField
            label={<><i>v</i> (value length)</>}
            description={<>Length of the value in the (<i>k</i>, <i>v</i>)-mer sketch. Must be ≤ 31.</>}
            value={sketch.v}
            onChange={(v) => { setSketch("v", v); setAnalyze("v", v); }}
            min={1} max={31} disabled={disabled}
          />
          <NumField
            label="Trim at the front of the reads (bp)"
            description="Number of bases to discard from the start of each read before sketching. Useful to remove low-quality leading bases."
            value={sketch.trim_front}
            onChange={(v) => setSketch("trim_front", v)}
            min={0} disabled={disabled}
          />
          <NumField
            label="Trim at the back of the reads (bp)"
            description="Number of bases to discard from the end of each read before sketching. Useful to remove low-quality trailing bases."
            value={sketch.trim_back}
            onChange={(v) => setSketch("trim_back", v)}
            min={0} disabled={disabled}
          />
          <NumFieldOpt
            label={<>Subsampling parameter <i>c</i></>}
            description={<>Retain only 1-in-<i>c</i> keys to reduce memory usage. If not set, this is chosen automatically from the estimated input size.</>}
            value={sketch.c}
            onChange={(v) => setSketch("c", v)}
            min={1} placeholder="auto" disabled={disabled}
          />
          <CheckField
            label="Use only the forward strand"
            description="Sketch only the forward strand of each read, not their reverse complements."
            value={sketch.forward_only}
            onChange={(v) => { setSketch("forward_only", v); setAnalyze("forward_only", v); }}
            disabled={disabled}
          />
        </div>
        <div>
          <h4>Analyze</h4>
          <NumField
            label={<>Ignore smallest <i>t</i></>}
            description={<>Exclude the smallest <i>t</i> values from hazard model fitting, as the smallest <i>t</i> values may be unreliable in error spectrum analysis.</>}
            value={analyze.ignore_smallest_t}
            onChange={(v) => setAnalyze("ignore_smallest_t", v)}
            min={0} disabled={disabled}
          />
          <NumField
            label={<>Ignore largest <i>t</i></>}
            description={<>Exclude the largest <i>t</i> values from hazard model fitting, as the largest <i>t</i> values may be unreliable in error spectrum analysis.</>}
            value={analyze.ignore_largest_t}
            onChange={(v) => setAnalyze("ignore_largest_t", v)}
            min={0} disabled={disabled}
          />
          <NumField
            label="Outlier threshold"
            description={<>P-value threshold for identifying and removing outlier keys before fitting the error model.</>}
            value={analyze.outlier_threshold}
            onChange={(v) => setAnalyze("outlier_threshold", v)}
            min={0} step={1e-12} disabled={disabled}
          />
          <NumFieldOpt
            label="Lower bound for key coverage"
            description={<>Minimum coverage a key must have to be included in the analysis. Leave blank to use the Rust default: 0 with a reference, 10 otherwise.</>}
            value={analyze.lower_bound}
            onChange={(v) => setAnalyze("lower_bound", v)}
            min={0} placeholder="auto" disabled={disabled}
          />
          <div className="field-row">
            <FieldLabel
              label="Hazard model"
              description={<>Statistical model used to fit the hazard rate curve. "Weibull" allows a non-constant hazard; "Constant" assumes uniform error rate across positions.</>}
            />
            <select
              value={analyze.hazard_model}
              disabled={disabled}
              onChange={(e) => setAnalyze("hazard_model", e.target.value)}
            >
              <option value="weibull">Weibull</option>
              <option value="constant">Constant</option>
            </select>
          </div>
          <CheckField
            label="Use all (no outlier filter)"
            description="Disable outlier filtering and use all keys for model fitting. Recommended to set to true if a reference is provided."
            value={analyze.use_all}
            onChange={(v) => setAnalyze("use_all", v)}
            disabled={disabled}
          />
        </div>
      </div>
    </div>
  );
}
