import { useEffect, useRef } from "react";
import { LogLine } from "../../types";

interface Props {
  logs: LogLine[];
}

export function LogConsole({ logs }: Props) {
  const bottomRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    bottomRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [logs]);

  if (logs.length === 0) return null;

  return (
    <div className="log-console">
      {logs.map((log, i) => (
        <div key={i} className={`log-line log-${log.stream}${log.stream === "stdout" && log.line.startsWith("$ ") ? " log-cmd" : ""}`}>
          {log.stream === "exit"
            ? `[exit code ${log.exit_code}]`
            : log.line}
        </div>
      ))}
      <div ref={bottomRef} />
    </div>
  );
}
