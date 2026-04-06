import { useState, useEffect } from "react";
import { invoke } from "@tauri-apps/api/core";
import Papa from "papaparse";

export function useCSVData<T>(path: string | null): {
  data: T[] | null;
  loading: boolean;
  error: string | null;
} {
  const [data, setData] = useState<T[] | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!path) return;
    setLoading(true);
    setError(null);
    setData(null);

    invoke<string>("read_csv_file", { path })
      .then((csv) => {
        const result = Papa.parse<T>(csv, {
          header: true,
          dynamicTyping: true,
          skipEmptyLines: true,
        });
        setData(result.data);
      })
      .catch((e) => setError(String(e)))
      .finally(() => setLoading(false));
  }, [path]);

  return { data, loading, error };
}
