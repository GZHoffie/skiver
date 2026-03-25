import { invoke } from "@tauri-apps/api/core";

interface Props {
  outputDir: string;
  outputPrefix: string;
  onDirChange: (d: string) => void;
  onPrefixChange: (p: string) => void;
  disabled: boolean;
  dirError?: boolean;
  prefixError?: boolean;
}

export function OutputConfig({
  outputDir,
  outputPrefix,
  onDirChange,
  onPrefixChange,
  disabled,
  dirError,
  prefixError,
}: Props) {
  async function pick() {
    const dir = await invoke<string>("pick_output_dir");
    if (dir) onDirChange(dir);
  }

  return (
    <div className="section">
      <h3>Output</h3>
      <div className="field-row">
        <label>Directory</label>
        <button onClick={pick} disabled={disabled}>
          Browse…
        </button>
        <span className={`path-display${dirError ? " error" : ""}`}>{outputDir || "Not selected"}</span>
      </div>
      <div className="field-row">
        <label>Prefix</label>
        <input
          type="text"
          value={outputPrefix}
          onChange={(e) => onPrefixChange(e.target.value)}
          disabled={disabled}
          placeholder="e.g. sample1"
          className={prefixError ? "error" : undefined}
        />
      </div>
    </div>
  );
}
