import { invoke } from "@tauri-apps/api/core";

interface Props {
  files: string[];
  onChange: (files: string[]) => void;
  disabled: boolean;
}

export function FileSelector({ files, onChange, disabled }: Props) {
  async function pick() {
    const picked = await invoke<string[]>("pick_input_files");
    if (picked && picked.length > 0) {
      onChange([...files, ...picked.filter((f) => !files.includes(f))]);
    }
  }

  function remove(idx: number) {
    onChange(files.filter((_, i) => i !== idx));
  }

  return (
    <div className="section">
      <h3>Input Files</h3>
      <p className="hint">FASTA/FASTQ files or pre-computed .kvmer sketches</p>
      <button onClick={pick} disabled={disabled}>
        + Add files…
      </button>
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
