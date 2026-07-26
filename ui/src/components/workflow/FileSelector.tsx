import { useState } from "react";
import { invoke } from "@tauri-apps/api/core";

interface Props {
  files: string[];
  onChange: (files: string[]) => void;
  disabled: boolean;
  error?: boolean;
}

export function FileSelector({ files, onChange, disabled, error }: Props) {
  const [pattern, setPattern] = useState("");

  async function pick() {
    const picked = await invoke<string[]>("pick_input_files");
    if (picked && picked.length > 0) {
      onChange([...files, ...picked.filter((f) => !files.includes(f))]);
    }
  }

  function addPattern() {
    const p = pattern.trim();
    if (!p || files.includes(p)) return;
    onChange([...files, p]);
    setPattern("");
  }

  function remove(idx: number) {
    onChange(files.filter((_, i) => i !== idx));
  }

  return (
    <div className={`section${error ? " error" : ""}`}>
      <h3>Input Files</h3>
      <p className="hint">
        One pre-computed sketch (any extension), or FASTA/FASTQ/BAM files and glob patterns
      </p>
      <div style={{ display: "flex", gap: "8px" }}>
        <input
          type="text"
          value={pattern}
          placeholder="Path or glob pattern"
          onChange={(e) => setPattern(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && addPattern()}
          disabled={disabled}
          style={{ flex: 1 }}
        />
        <button onClick={addPattern} disabled={disabled || !pattern.trim()}>Add</button>
        <button onClick={pick} disabled={disabled}>Browse…</button>
      </div>
      {files.length > 0 && (
        <ul className="file-list">
          {files.map((f, i) => (
            <li key={f}>
              <span className="file-name" title={f}>
                {f.split("/").pop()}
              </span>
              <button
                className="remove-btn"
                onClick={() => remove(i)}
                disabled={disabled}
              >
                ✕
              </button>
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}
