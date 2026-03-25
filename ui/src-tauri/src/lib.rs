mod commands;
mod types;

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .invoke_handler(tauri::generate_handler![
            commands::file_ops::pick_input_files,
            commands::file_ops::pick_reference_file,
            commands::file_ops::pick_output_dir,
            commands::file_ops::read_csv_file,
            commands::file_ops::file_exists,
            commands::run_sketch::run_sketch,
            commands::run_analyze::run_analyze,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
