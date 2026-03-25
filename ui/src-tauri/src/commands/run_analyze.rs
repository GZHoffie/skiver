use tauri::Emitter;
use tauri_plugin_shell::ShellExt;
use tauri_plugin_shell::process::CommandEvent;

use crate::types::{AnalyzeParams, LogEvent};

#[tauri::command]
pub async fn run_analyze(app: tauri::AppHandle, params: AnalyzeParams) -> Result<(), String> {
    let mut args: Vec<String> = vec!["analyze".into(), params.kvmer_path.clone()];

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
    args.extend(["--hazard-model".into(), params.hazard_model.clone()]);
    args.extend([
        "--num-experiments".into(),
        params.num_experiments.to_string(),
    ]);
    args.extend(["-e".into(), params.outlier_threshold.to_string()]);
    if let Some(lb) = params.lower_bound {
        args.extend(["-l".into(), lb.to_string()]);
    }
    if params.forward_only {
        args.push("--forward-only".into());
    }
    if params.use_all {
        args.push("--use-all".into());
    }

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
                return Ok(());
            }
            _ => continue,
        };
        app.emit("skiver-log", log).ok();
    }
    Ok(())
}
