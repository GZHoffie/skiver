use tauri::Emitter;
use tauri_plugin_shell::ShellExt;
use tauri_plugin_shell::process::CommandEvent;

use crate::types::{LogEvent, SketchParams};

#[tauri::command]
pub async fn run_sketch(app: tauri::AppHandle, params: SketchParams) -> Result<(), String> {
    let mut args: Vec<String> = vec!["sketch".into()];

    for f in &params.input_files {
        args.push(f.clone());
    }
    args.extend(["-k".into(), params.k.to_string()]);
    args.extend(["-v".into(), params.v.to_string()]);
    args.extend(["-o".into(), params.output_path.clone()]);
    args.extend(["-f".into(), params.trim_front.to_string()]);
    args.extend(["-b".into(), params.trim_back.to_string()]);
    if let Some(c) = params.c {
        args.extend(["-c".into(), c.to_string()]);
    }
    if params.forward_only {
        args.push("--forward-only".into());
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
                    return Err(format!("skiver sketch exited with code {:?}", code));
                }
                return Ok(());
            }
            _ => continue,
        };
        app.emit("skiver-log", log).ok();
    }
    Ok(())
}
