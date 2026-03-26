import { invoke } from "@tauri-apps/api/core";

interface Props {
  reference: string | null;
  onChange: (ref: string | null) => void;
  disabled: boolean;
}

export function ReferenceSelector({ reference, onChange, disabled }: Props) {
  async function pick() {
    try {
      const picked = await invoke<string>("pick_reference_file");
      if (picked) onChange(picked);
    } catch {
      // user cancelled
    }
  }

  return (
    <div className="section">
      <h3>Reference File (optional)</h3>
      <p className="hint">FASTA/FASTA.gz reference genome. Supports glob patterns.</p>
      <div style={{ display: "flex", gap: "8px" }}>
        <input
          type="text"
          value={reference ?? ""}
          placeholder="Path or glob pattern…"
          onChange={(e) => onChange(e.target.value || null)}
          disabled={disabled}
          style={{ flex: 1 }}
        />
        <button onClick={pick} disabled={disabled}>Browse…</button>
        {reference && (
          <button onClick={() => onChange(null)} disabled={disabled}>✕</button>
        )}
      </div>
    </div>
  );
}
