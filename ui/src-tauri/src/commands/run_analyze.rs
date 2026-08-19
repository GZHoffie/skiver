use tauri::Emitter;
use tauri_plugin_shell::process::CommandEvent;
use tauri_plugin_shell::ShellExt;

use std::path::Path;

use crate::types::{AnalyzeParams, AnalyzeResult, LogEvent};

fn csv_has_data_rows(path: &Path) -> bool {
    std::fs::read_to_string(path)
        .map(|csv| csv.lines().skip(1).any(|line| !line.trim().is_empty()))
        .unwrap_or(false)
}

#[tauri::command]
pub async fn run_analyze(
    app: tauri::AppHandle,
    params: AnalyzeParams,
) -> Result<AnalyzeResult, String> {
    let mut args: Vec<String> = vec!["analyze".into()];

    for f in &params.input_files {
        args.push(f.clone());
    }

    args.extend(["-k".into(), params.k.to_string()]);
    args.extend(["-v".into(), params.v.to_string()]);
    args.extend(["-o".into(), params.output_prefix.clone()]);
    args.extend([
        "--ignore-largest-t".into(),
        params.ignore_largest_t.to_string(),
    ]);
    args.extend([
        "--ignore-smallest-t".into(),
        params.ignore_smallest_t.to_string(),
    ]);
    args.extend(["--threads".into(), params.threads.to_string()]);
    args.extend(["--hazard-model".into(), params.hazard_model.clone()]);
    args.extend([
        "--num-experiments".into(),
        params.num_experiments.to_string(),
    ]);
    args.extend(["-e".into(), params.outlier_threshold.to_string()]);
    if let Some(lb) = params.lower_bound {
        args.extend(["-l".into(), lb.to_string()]);
    }
    if let Some(ref r) = params.reference {
        args.extend(["-r".into(), r.clone()]);
    }
    if params.forward_only {
        args.push("--forward-only".into());
    }
    if params.use_all {
        args.push("--use-all".into());
    }

    app.emit(
        "skiver-log",
        LogEvent {
            stream: "stdout".into(),
            line: format!("$ skiver {}", args.join(" ")),
            exit_code: None,
        },
    )
    .ok();

    let (mut rx, _child) = app
        .shell()
        .sidecar("skiver")
        .map_err(|e| e.to_string())?
        .args(&args)
        .spawn()
        .map_err(|e| e.to_string())?;

    while let Some(event) = rx.recv().await {
        let log = match event {
            CommandEvent::Stdout(line) => LogEvent {
                stream: "stdout".into(),
                line: String::from_utf8_lossy(&line).into_owned(),
                exit_code: None,
            },
            CommandEvent::Stderr(line) => LogEvent {
                stream: "stderr".into(),
                line: String::from_utf8_lossy(&line).into_owned(),
                exit_code: None,
            },
            CommandEvent::Terminated(payload) => {
                let code = payload.code;
                app.emit(
                    "skiver-log",
                    LogEvent {
                        stream: "exit".into(),
                        line: String::new(),
                        exit_code: code,
                    },
                )
                .ok();
                if code != Some(0) {
                    return Err(format!("skiver analyze exited with code {:?}", code));
                }
                let phred_path = format!("{}.summary_phred.csv", params.output_prefix);
                return Ok(AnalyzeResult {
                    has_qscore_data: csv_has_data_rows(Path::new(&phred_path)),
                });
            }
            _ => continue,
        };
        app.emit("skiver-log", log).ok();
    }
    Err("skiver analyze ended without reporting an exit status".into())
}

#[cfg(test)]
mod tests {
    use super::csv_has_data_rows;
    use std::fs;

    #[test]
    fn detects_optional_csv_data_rows() {
        let dir = std::env::temp_dir();
        let empty_path = dir.join(format!("skiver-empty-phred-{}.csv", std::process::id()));
        let data_path = dir.join(format!("skiver-data-phred-{}.csv", std::process::id()));

        fs::write(&empty_path, "qscore,num_correct,num_error\n").unwrap();
        fs::write(&data_path, "qscore,num_correct,num_error\n20,100,1\n").unwrap();

        assert!(!csv_has_data_rows(&empty_path));
        assert!(csv_has_data_rows(&data_path));

        fs::remove_file(empty_path).unwrap();
        fs::remove_file(data_path).unwrap();
    }
}
