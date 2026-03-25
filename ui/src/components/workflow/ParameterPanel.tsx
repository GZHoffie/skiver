import { SketchParams, AnalyzeParams } from "../../types";

interface Props {
  sketch: SketchParams;
  analyze: AnalyzeParams;
  onSketchChange: (p: SketchParams) => void;
  onAnalyzeChange: (p: AnalyzeParams) => void;
  disabled: boolean;
}

function NumField({
  label,
  value,
  onChange,
  min,
  max,
  step,
  disabled,
}: {
  label: string;
  value: number;
  onChange: (v: number) => void;
  min?: number;
  max?: number;
  step?: number;
  disabled: boolean;
}) {
  return (
    <div className="field-row">
      <label>{label}</label>
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

function CheckField({
  label,
  value,
  onChange,
  disabled,
}: {
  label: string;
  value: boolean;
  onChange: (v: boolean) => void;
  disabled: boolean;
}) {
  return (
    <div className="field-row">
      <label>{label}</label>
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
          <NumField label="k (key length)" value={sketch.k} onChange={(v) => { setSketch("k", v); setAnalyze("k", v); }} min={1} max={31} disabled={disabled} />
          <NumField label="v (value length)" value={sketch.v} onChange={(v) => { setSketch("v", v); setAnalyze("v", v); }} min={1} max={31} disabled={disabled} />
          <NumField label="Trim front (bp)" value={sketch.trim_front} onChange={(v) => setSketch("trim_front", v)} min={0} disabled={disabled} />
          <NumField label="Trim back (bp)" value={sketch.trim_back} onChange={(v) => setSketch("trim_back", v)} min={0} disabled={disabled} />
          <CheckField label="Forward only" value={sketch.forward_only} onChange={(v) => { setSketch("forward_only", v); setAnalyze("forward_only", v); }} disabled={disabled} />
        </div>
        <div>
          <h4>Analyze</h4>
          <NumField label="Ignore smallest t" value={analyze.ignore_smallest_t} onChange={(v) => setAnalyze("ignore_smallest_t", v)} min={0} disabled={disabled} />
          <NumField label="Ignore largest t" value={analyze.ignore_largest_t} onChange={(v) => setAnalyze("ignore_largest_t", v)} min={0} disabled={disabled} />
          <NumField label="Bootstrap experiments" value={analyze.num_experiments} onChange={(v) => setAnalyze("num_experiments", v)} min={1} disabled={disabled} />
          <div className="field-row">
            <label>Hazard model</label>
            <select
              value={analyze.hazard_model}
              disabled={disabled}
              onChange={(e) => setAnalyze("hazard_model", e.target.value)}
            >
              <option value="weibull">Weibull</option>
              <option value="constant">Constant</option>
            </select>
          </div>
          <CheckField label="Use all (no outlier filter)" value={analyze.use_all} onChange={(v) => setAnalyze("use_all", v)} disabled={disabled} />
        </div>
      </div>
    </div>
  );
}
