use tauri_plugin_dialog::DialogExt;

#[tauri::command]
pub async fn pick_input_files(app: tauri::AppHandle) -> Result<Vec<String>, String> {
    let paths = app
        .dialog()
        .file()
        .blocking_pick_files()
        .ok_or("No files selected")?;
    Ok(paths
        .into_iter()
        .filter_map(|p| p.into_path().ok())
        .map(|p| p.to_string_lossy().into_owned())
        .collect())
}

#[tauri::command]
pub async fn pick_output_dir(app: tauri::AppHandle) -> Result<String, String> {
    let path = app
        .dialog()
        .file()
        .blocking_pick_folder()
        .ok_or("No directory selected")?;
    path.into_path()
        .map(|p| p.to_string_lossy().into_owned())
        .map_err(|e| e.to_string())
}

#[tauri::command]
pub async fn pick_reference_file(app: tauri::AppHandle) -> Result<String, String> {
    let path = app
        .dialog()
        .file()
        .add_filter("Reference genomes", &["fasta", "fa", "fna", "gz"])
        .blocking_pick_file()
        .ok_or("No file selected")?;
    path.into_path()
        .map(|p| p.to_string_lossy().into_owned())
        .map_err(|e| e.to_string())
}

#[tauri::command]
pub async fn read_csv_file(path: String) -> Result<String, String> {
    std::fs::read_to_string(&path).map_err(|e| e.to_string())
}

#[tauri::command]
pub async fn file_exists(path: String) -> bool {
    std::path::Path::new(&path).exists()
}
