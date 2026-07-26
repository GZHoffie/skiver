//! Generate the plots produced by `scripts/plot_all.py` as vector PDFs.

/**
 * NOTE: This module is created with the help of ChatGPT and verified by the author.
 */

use kuva::prelude::*;
use log::info;
use std::error::Error;
use std::fs;
use std::path::Path;

const SLATE: &str = "slategray";
const RED: &str = "indianred";

type Result<T> = std::result::Result<T, Box<dyn Error>>;

/// Plot settings used by `analyze`.
pub struct PlotOptions {
    pub normalize: bool,
    pub log_scale: bool,
    pub t_min: usize,
    pub t_max: usize,
    pub num_bases: usize,
    pub min_coverage: f64,
    pub min_bases: f64,
}

impl Default for PlotOptions {
    fn default() -> Self {
        Self {
            normalize: false,
            log_scale: false,
            t_min: 1,
            t_max: 100,
            num_bases: 100,
            min_coverage: 5000.0,
            min_bases: 5000.0,
        }
    }
}

#[derive(Clone)]
struct Table {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Table {
    fn read(path: &str) -> Result<Self> {
        let mut reader = csv::ReaderBuilder::new().flexible(true).from_path(path)?;
        let headers: Vec<String> = reader.headers()?.iter().map(str::to_owned).collect();
        let mut rows: Vec<Vec<String>> = reader
            .records()
            .collect::<std::result::Result<Vec<_>, _>>()?
            .into_iter()
            .map(|r| r.iter().map(str::to_owned).collect())
            .collect();
        // Some historical GC summaries omit the trailing empty field.
        for row in &mut rows {
            row.resize(headers.len(), String::new());
        }
        Ok(Self { headers, rows })
    }

    fn index(&self, name: &str) -> Result<usize> {
        self.headers
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| format!("missing CSV column {name:?}").into())
    }

    fn strings(&self, name: &str) -> Result<Vec<String>> {
        let i = self.index(name)?;
        Ok(self.rows.iter().map(|r| r[i].clone()).collect())
    }

    fn nums(&self, name: &str) -> Result<Vec<f64>> {
        let i = self.index(name)?;
        self.rows
            .iter()
            .enumerate()
            .map(|(row, r)| {
                r[i].parse::<f64>().map_err(|e| {
                    format!("invalid number in column {name:?}, row {}: {e}", row + 2).into()
                })
            })
            .collect()
    }

    fn num(&self, name: &str) -> Result<f64> {
        self.nums(name)?
            .first()
            .copied()
            .ok_or_else(|| format!("CSV has no rows for column {name:?}").into())
    }
}

fn output(prefix: &str, kind: &str) -> String {
    format!("{prefix}.plot_{kind}.pdf")
}

fn save_scene(path: &str, scene: kuva::render::render::Scene) -> Result<()> {
    let bytes = PdfBackend::new().render_scene(&scene)?;
    fs::write(path, bytes)?;
    Ok(())
}

fn save_one(path: &str, plots: Vec<Plot>, layout: Layout) -> Result<()> {
    save_scene(path, kuva::render::render::render_multiple(plots, layout))
}

fn line(data: impl IntoIterator<Item = (f64, f64)>, color: &str, label: &str) -> Plot {
    LinePlot::new()
        .with_data(data)
        .with_color(color)
        .with_stroke_width(2.5)
        .with_legend(label)
        .into()
}

fn band(x: Vec<f64>, lo: Vec<f64>, hi: Vec<f64>, color: &str, label: &str) -> Plot {
    BandPlot::new(x, lo, hi)
        .with_color(color)
        .with_opacity(0.25)
        .with_legend(label)
        .into()
}

fn layout(plots: &[Plot], title: &str, x: &str, y: &str, log_y: bool) -> Layout {
    let l = Layout::auto_from_plots(plots)
        .with_title(title)
        .with_x_label(x)
        .with_y_label(y)
        .with_show_grid(false)
        .with_legend_position(LegendPosition::InsideTopRight)
        .with_legend_box(false);
    if log_y { l.with_log_y() } else { l }
}

fn bin_edges(centers: &[f64]) -> Vec<f64> {
    match centers {
        [] => Vec::new(),
        [x] => vec![x - 0.5, x + 0.5],
        _ => {
            let mut edges = Vec::with_capacity(centers.len() + 1);
            edges.push(centers[0] - (centers[1] - centers[0]) / 2.0);
            edges.extend(centers.windows(2).map(|w| (w[0] + w[1]) / 2.0));
            let n = centers.len();
            edges.push(centers[n - 1] + (centers[n - 1] - centers[n - 2]) / 2.0);
            edges
        }
    }
}

fn scale_counts(counts: Vec<f64>, label: &str) -> (Vec<f64>, String) {
    let max_count = counts.iter().copied().reduce(f64::max).unwrap_or(0.0);
    let exponent = if max_count >= 1000.0 {
        ((max_count.log10().floor() as i32) / 3) * 3
    } else {
        0
    };
    if exponent == 0 {
        (counts, label.to_owned())
    } else {
        let scale = 10f64.powi(exponent);
        (
            counts.into_iter().map(|count| count / scale).collect(),
            format!("{label} (×10^{exponent})"),
        )
    }
}

fn split_ci(s: &str) -> Result<(f64, f64)> {
    let (lo, hi) = s
        .split_once('~')
        .ok_or_else(|| format!("invalid percentile interval {s:?}"))?;
    Ok((lo.parse()?, hi.parse()?))
}

fn plot_spectrum(prefix: &str, normalize: bool) -> Result<()> {
    let spectrum = Table::read(&format!("{prefix}.summary_error_spectrum.csv"))?;
    let operations = spectrum.strings("operation")?;
    let target = if normalize {
        1.0
    } else {
        Table::read(&format!("{prefix}.summary_error_rate.csv"))?.num("per_base_error_rate")?
    };
    let suffix = if normalize { " (normalized)" } else { "" };
    let mut panels = Vec::new();
    let mut layouts = Vec::new();

    for (column, strand) in [("total", "both strands"), ("forward", "forward strand")] {
        let counts = spectrum.nums(column)?;
        let total: f64 = counts.iter().sum();
        let scale = if total > 0.0 {
            target * 100.0 / total
        } else {
            0.0
        };
        let mut matrix = vec![vec![0.0; 5]; 5];
        let mut types = [0.0; 3];
        for (op, count) in operations.iter().zip(counts) {
            let (from, to) = op.split_once('>').unwrap_or(("-", "-"));
            let bi = |b: &str| {
                ["A", "C", "G", "T", "-"]
                    .iter()
                    .position(|x| *x == b)
                    .unwrap_or(4)
            };
            matrix[bi(from)][bi(to)] += count * scale;
            if op.starts_with("->") {
                types[0] += count * scale;
            } else if op.ends_with(">-") {
                types[1] += count * scale;
            } else {
                types[2] += count * scale;
            }
        }
        // Hide the no-error diagonal without changing the off-diagonal scale.
        for (i, row) in matrix.iter_mut().enumerate() {
            row[i] = 0.0;
        }
        let spectrum_map = ColorMap::Custom(std::sync::Arc::new(|t: f64| {
            let t = t.clamp(0.0, 1.0);
            let channel = |light: u8, dark: u8| {
                (light as f64 + (dark as f64 - light as f64) * t).round() as u8
            };
            format!(
                "rgb({},{},{})",
                channel(255, 112),
                channel(255, 128),
                channel(255, 144)
            )
        }));
        let heat: Plot = Heatmap::new()
            .with_data(matrix.clone())
            .with_color_map(spectrum_map)
            .into();
        let hp = vec![heat];
        let mut hl = Layout::auto_from_plots(&hp)
            .with_title(format!("Error spectrum{suffix} — {strand}"))
            .with_x_label("Observed base")
            .with_y_label("Original base")
            .with_x_categories(["A", "C", "G", "T", "-"].map(str::to_owned).to_vec())
            .with_y_categories(["A", "C", "G", "T", "-"].map(str::to_owned).to_vec());
        for (row, values) in matrix.iter().enumerate() {
            for (column, value) in values.iter().enumerate() {
                if row != column {
                    hl = hl.with_annotation(
                        TextAnnotation::new(
                            format!("{value:.3}"),
                            column as f64 + 1.0,
                            // TextAnnotation places un-arrowed text six pixels above its
                            // anchor. Shift the data-space anchor down slightly so the
                            // visible label is centered in the heatmap cell.
                            row as f64 + 0.93,
                        )
                        .with_font_size(11),
                    );
                }
            }
        }
        panels.push(hp);
        layouts.push(hl);

        let bars: Plot = BarPlot::new()
            .with_bars([
                ("Insertion", types[0]),
                ("Deletion", types[1]),
                ("Substitution", types[2]),
            ])
            .with_color(SLATE)
            .into();
        let bp = vec![bars];
        let max_value = types.iter().copied().reduce(f64::max).unwrap_or(0.0);
        let mut bl = layout(
            &bp,
            &format!("Error type distribution{suffix} — {strand}"),
            "Error type",
            if normalize {
                "Proportion (%)"
            } else {
                "Error rate (%)"
            },
            false,
        );
        for (i, value) in types.iter().enumerate() {
            bl = bl.with_annotation(
                TextAnnotation::new(format!("{value:.3}"), i as f64 + 1.0, value * 1.04)
                    .with_font_size(12),
            );
        }
        if max_value > 0.0 {
            bl.y_range.1 = max_value * 1.15;
            bl.data_y_range = Some((0.0, max_value * 1.15));
        }
        panels.push(bp);
        layouts.push(bl);
    }
    // Reorder from [heat0, bar0, heat1, bar1] into the intended 2×2 grid.
    let scene = Figure::new(2, 2)
        .with_plots(panels)
        .with_layouts(layouts)
        .with_cell_size(520.0, 400.0)
        .render();
    save_scene(&output(prefix, "spectrum"), scene)
}

fn percentile(mut values: Vec<f64>, p: f64) -> f64 {
    values.sort_by(f64::total_cmp);
    if values.is_empty() {
        return f64::NAN;
    }
    let at = p * (values.len() - 1) as f64;
    let lo = at.floor() as usize;
    let hi = at.ceil() as usize;
    values[lo] + (values[hi] - values[lo]) * (at - lo as f64)
}

fn histogram_counts(values: &[f64], min: f64, max: f64, bins: usize) -> Vec<f64> {
    let width = (max - min) / bins as f64;
    let mut counts = vec![0.0; bins];
    for &value in values {
        let index = if value == max {
            bins - 1
        } else {
            ((value - min) / width).floor() as usize
        };
        if index < bins {
            counts[index] += 1.0;
        }
    }
    counts
}

fn plot_coverage(prefix: &str) -> Result<()> {
    let report = Table::read(&format!("{prefix}.summary_error_rate.csv"))?;
    let kv = Table::read(&format!("{prefix}.kvmer.csv"))?;
    let lambda = report.num("lambda")?;
    let beta = report.num("beta")?;
    let key_len = kv
        .strings("key")?
        .first()
        .map(String::len)
        .ok_or("empty kvmer CSV")?;
    let survival = (-lambda * (key_len as f64).powf(beta)).exp();
    let all: Vec<f64> = kv
        .nums("total_count")?
        .into_iter()
        .map(|v| v / survival)
        .collect();
    let pass = kv.strings("passes_filter")?;
    let passing: Vec<f64> = all
        .iter()
        .zip(pass)
        .filter(|(_, p)| p == "true" || p == "True" || p == "1")
        .map(|(v, _)| *v)
        .collect();
    let threshold = percentile(all.clone(), 0.9999);
    let all: Vec<_> = all.into_iter().filter(|v| *v <= threshold).collect();
    let passing: Vec<_> = passing.into_iter().filter(|v| *v <= threshold).collect();
    let min = all.iter().copied().reduce(f64::min).unwrap();
    let max = all.iter().copied().reduce(f64::max).unwrap();
    let bins = 100;
    let edges: Vec<f64> = (0..=bins)
        .map(|i| min + (max - min) * i as f64 / bins as f64)
        .collect();
    let log_counts = |values: &[f64]| {
        histogram_counts(values, min, max, bins)
            .into_iter()
            .map(|count| if count > 0.0 { count.log10() } else { 0.0 })
            .collect::<Vec<f64>>()
    };
    let all_counts = log_counts(&all);
    let passing_counts = log_counts(&passing);
    let max_log_count = all_counts
        .iter()
        .chain(&passing_counts)
        .copied()
        .reduce(f64::max)
        .unwrap();
    let plots: Vec<Plot> = vec![
        Histogram::from_bins(edges.clone(), all_counts)
            .with_color("#70809080")
            .with_legend("All keys")
            .into(),
        Histogram::from_bins(edges, passing_counts)
            .with_color("#4682b4b0")
            .with_legend("Passing filter")
            .into(),
    ];
    let mut l = layout(
        &plots,
        "Estimated true coverage",
        "Coverage",
        "Count",
        false,
    )
    .with_ticks(10)
    .with_x_tick_format(TickFormat::Integer)
    .with_legend_entries(vec![
        LegendEntry {
            label: "All keys".into(),
            color: "#70809080".into(),
            shape: LegendShape::Rect,
            dasharray: None,
        },
        LegendEntry {
            label: "Passing filter".into(),
            color: "#4682b4b0".into(),
            shape: LegendShape::Rect,
            dasharray: None,
        },
    ])
    .with_y_tick_format(TickFormat::Custom(std::sync::Arc::new(|exponent| {
        let count = 10f64.powf(exponent);
        if count >= 1.0 {
            format!("{count:.0}")
        } else {
            format!("{count:.2}")
        }
    })));
    l.x_range.0 = 0.0;
    l.data_x_range = Some((0.0, max));
    l.y_range = (0.0, max_log_count * 1.05);
    l.data_y_range = Some(l.y_range);
    save_one(&output(prefix, "coverage"), plots, l)
}

fn plot_hazard(prefix: &str, t_min: usize, t_max: usize, log_y: bool) -> Result<()> {
    let report = Table::read(&format!("{prefix}.summary_error_rate.csv"))?;
    let hr = Table::read(&format!("{prefix}.hazard_rate.csv"))?;
    let lambda = report.num("lambda")?;
    let beta = report.num("beta")?;
    let t = hr.nums("t")?;
    let empirical = hr.nums("hazard_ratio")?;
    let lo = hr.nums("5th_percentile")?;
    let hi = hr.nums("95th_percentile")?;
    let fitted: Vec<f64> = t
        .iter()
        .map(|v| 1.0 - (-lambda * (v.powf(beta) - (v - 1.0).powf(beta))).exp())
        .collect();
    let hp = vec![
        band(t.clone(), lo, hi, SLATE, "5%-95% percentile"),
        line(
            t.iter().copied().zip(empirical),
            SLATE,
            "Estimated hazard rate",
        ),
        LinePlot::new()
            .with_data(t.iter().copied().zip(fitted))
            .with_color(RED)
            .with_stroke_width(2.5)
            .with_line_style(LineStyle::Dashed)
            .with_legend("Fitted hazard rate")
            .into(),
    ];
    let hl = layout(&hp, "Estimated hazard rate", "t", "h(t)", log_y);
    let survival: Vec<(f64, f64)> = (t_min..=t_max)
        .map(|v| (v as f64, (-lambda * (v as f64).powf(beta)).exp()))
        .collect();
    let sp = vec![
        LinePlot::new()
            .with_data(survival)
            .with_color(RED)
            .with_stroke_width(2.5)
            .with_line_style(LineStyle::Dashed)
            .with_legend("Estimated survival rate")
            .into(),
    ];
    let sl = layout(&sp, "Estimated survival rate", "t", "S(t)", log_y);
    let scene = Figure::new(1, 2)
        .with_plots(vec![hp, sp])
        .with_layouts(vec![hl, sl])
        .with_cell_size(500.0, 380.0)
        .render();
    save_scene(&output(prefix, "hazard_survival"), scene)
}

fn filtered_rate_data(
    df: &Table,
    x_col: &str,
    min_count: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    let x = df.nums(x_col)?;
    let correct = df.nums("num_correct")?;
    let error = df.nums("num_error")?;
    let rate_col = if df.headers.iter().any(|h| h == "per_base_error_rate") {
        "per_base_error_rate"
    } else {
        "per_base_error_rate_median"
    };
    let rate = df.nums(rate_col)?;
    let ci = df.strings("per_base_error_rate_5-95th_percentile")?;
    let mut out = (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for i in 0..x.len() {
        let count = correct[i] + error[i];
        if count >= min_count {
            let (lo, hi) = split_ci(&ci[i])?;
            out.0.push(x[i]);
            out.1.push(rate[i]);
            out.2.push(lo);
            out.3.push(hi);
            out.4.push(count);
        }
    }
    Ok(out)
}

fn plot_qscore(prefix: &str, log_y: bool, min_count: f64) -> Result<()> {
    let df = Table::read(&format!("{prefix}.summary_phred.csv"))?;
    let (x, rate, lo, hi, counts) = filtered_rate_data(&df, "qscore", min_count)?;
    let max_q = x
        .iter()
        .copied()
        .reduce(f64::max)
        .ok_or("no qscore rows pass --min-coverage")?;
    let theory = (0..300).map(|i| {
        let q = 1.0 + (max_q - 1.0) * i as f64 / 299.0;
        (q, 10f64.powf(-q / 10.0))
    });
    let main = vec![
        band(x.clone(), lo, hi, SLATE, "5%-95% percentile"),
        line(x.iter().copied().zip(rate), SLATE, "Empirical error rate"),
        LinePlot::new()
            .with_data(theory)
            .with_color("black")
            .with_line_style(LineStyle::Dashed)
            .with_stroke_width(2.5)
            .with_legend("Theoretical (10^-Q/10)")
            .into(),
    ];
    let mut ml = layout(
        &main,
        "Quality-score calibration",
        "Phred quality score (Q)",
        "Error rate",
        log_y,
    );
    let edges = bin_edges(&x);
    let (counts, count_label) = scale_counts(counts, "# bases used for estimation");
    let hist = vec![
        Histogram::from_bins(edges.clone(), counts)
            .with_color(SLATE)
            .into(),
    ];
    let mut hl = layout(&hist, "", "Phred quality score (Q)", &count_label, false)
        .with_x_tick_step(2.0)
        .with_x_tick_format(TickFormat::Integer);
    let shared_range = (edges[0], *edges.last().unwrap());
    ml.x_range = shared_range;
    hl.x_range = shared_range;
    ml.data_x_range = Some(shared_range);
    hl.data_x_range = Some(shared_range);
    let scene = Figure::new(2, 1)
        .with_plots(vec![main, hist])
        .with_layouts(vec![ml, hl])
        .with_shared_x(0)
        .with_spacing(5.0)
        .with_row_height(0, 500.0)
        .with_row_height(1, 220.0)
        .with_figure_size(800.0, 700.0)
        .render();
    save_scene(&output(prefix, "qscore_calibration"), scene)
}

fn plot_gc(prefix: &str, log_y: bool, min_count: f64) -> Result<()> {
    let mut df = Table::read(&format!("{prefix}.summary_gc_content.csv"))?;
    // Older files have num_correct/num_error shifted one column to the right.
    if df.headers.iter().any(|h| h == "beta") {
        let e = df.index("num_error")?;
        if df.rows.iter().all(|r| r[e].is_empty()) {
            let b = df.index("beta")?;
            let c = df.index("num_correct")?;
            for row in &mut df.rows {
                row[e] = row[c].clone();
                row[c] = row[b].clone();
            }
        }
    }
    let min = df.nums("gc_content_min")?;
    let max = df.nums("gc_content_max_exclusive")?;
    let mid: Vec<f64> = min
        .iter()
        .zip(&max)
        .map(|(a, b)| (a + b - 1.0) / 2.0)
        .collect();
    let mut with_mid = df.clone();
    with_mid.headers.push("gc_mid".into());
    for (row, value) in with_mid.rows.iter_mut().zip(&mid) {
        row.push(value.to_string());
    }
    let (x, rate, lo, hi, _) = filtered_rate_data(&with_mid, "gc_mid", min_count)?;
    let correct = df.nums("num_correct")?;
    let error = df.nums("num_error")?;
    let main = vec![
        band(x.clone(), lo, hi, SLATE, "5%-95% percentile"),
        line(x.iter().copied().zip(rate), SLATE, "Empirical error rate"),
    ];
    let mut ml = layout(
        &main,
        "Error rate vs. GC content",
        "GC content (%)",
        "Error rate",
        log_y,
    );
    let counts: Vec<f64> = correct
        .iter()
        .zip(error)
        .filter_map(|(c, e)| {
            let total = c + e;
            (total >= min_count).then_some(total)
        })
        .collect();
    let (counts, count_label) = scale_counts(counts, "# bases");
    let edges = bin_edges(&x);
    let hist = vec![Histogram::from_bins(edges, counts).with_color(SLATE).into()];
    let mut hl = layout(&hist, "", "GC content (%)", &count_label, false)
        .with_x_tick_step(5.0)
        .with_x_tick_format(TickFormat::Integer);
    let bin_width = max[0] - min[0];
    let shared_range = (
        x.iter().copied().reduce(f64::min).unwrap() - bin_width / 2.0,
        x.iter().copied().reduce(f64::max).unwrap() + bin_width / 2.0,
    );
    ml.x_range = shared_range;
    hl.x_range = shared_range;
    ml.data_x_range = Some(shared_range);
    hl.data_x_range = Some(shared_range);
    let scene = Figure::new(2, 1)
        .with_plots(vec![main, hist])
        .with_layouts(vec![ml, hl])
        .with_shared_x(0)
        .with_spacing(5.0)
        .with_row_height(0, 500.0)
        .with_row_height(1, 220.0)
        .with_figure_size(800.0, 700.0)
        .render();
    save_scene(&output(prefix, "gc_content"), scene)
}

fn wilson(errors: f64, total: f64) -> (f64, f64) {
    if total <= 0.0 {
        return (0.0, 0.0);
    }
    let z = 1.6448536269514722_f64;
    let p = errors / total;
    let d = 1.0 + z * z / total;
    let center = (p + z * z / (2.0 * total)) / d;
    let half = z * (p * (1.0 - p) / total + z * z / (4.0 * total * total)).sqrt() / d;
    ((center - half).max(0.0), (center + half).min(1.0))
}

fn moving_average(y: &[f64], window: usize) -> Vec<f64> {
    (0..y.len())
        .map(|i| {
            let radius = window / 2;
            let lo = i.saturating_sub(radius);
            let hi = (i + radius + 1).min(y.len());
            y[lo..hi].iter().sum::<f64>() / (hi - lo) as f64
        })
        .collect()
}

fn plot_read_position(prefix: &str, num_bases: usize) -> Result<()> {
    let df = Table::read(&format!("{prefix}.summary_read_position.csv"))?;
    let index = df.nums("index")?;
    let from_start = df.strings("from_start")?;
    let correct = df.nums("num_correct")?;
    let errors = df.nums("num_error")?;
    let rate = df.nums("error_rate")?;
    let mut panels = Vec::new();
    let mut layouts = Vec::new();
    for (want_start, name) in [(true, "from start"), (false, "from end")] {
        let ids: Vec<usize> = (0..index.len())
            .filter(|i| {
                let flag =
                    from_start[*i] == "true" || from_start[*i] == "True" || from_start[*i] == "1";
                flag == want_start && index[*i] <= num_bases as f64
            })
            .collect();
        let x: Vec<_> = ids.iter().map(|i| index[*i]).collect();
        let y: Vec<_> = ids.iter().map(|i| rate[*i]).collect();
        let counts: Vec<_> = ids.iter().map(|i| correct[*i] + errors[*i]).collect();
        let intervals: Vec<_> = ids
            .iter()
            .map(|i| wilson(errors[*i], correct[*i] + errors[*i]))
            .collect();
        let mut plots = vec![
            band(
                x.clone(),
                intervals.iter().map(|v| v.0).collect(),
                intervals.iter().map(|v| v.1).collect(),
                SLATE,
                "5%-95% Wilson interval",
            ),
            line(
                x.iter().copied().zip(y.iter().copied()),
                SLATE,
                "Error rate",
            ),
        ];
        let window = (num_bases / 10).max(3);
        if y.len() >= window {
            plots.push(line(
                x.iter().copied().zip(moving_average(&y, window)),
                RED,
                &format!("Smoothed (window={window})"),
            ));
        }
        let l = layout(
            &plots,
            &format!("Error rate ({name})"),
            "Position in read",
            "Error rate",
            false,
        );
        panels.push(plots);
        layouts.push(l);
        let bp = vec![
            Histogram::from_bins(bin_edges(&x), counts)
                .with_color(SLATE)
                .into(),
        ];
        let bl = layout(
            &bp,
            "",
            "Position in read",
            "# bases used for estimation",
            false,
        )
        .with_x_tick_step(10.0)
        .with_x_tick_format(TickFormat::Integer);
        panels.push(bp);
        layouts.push(bl);
    }
    // Current order is start-line, start-count, end-line, end-count; place lines on top.
    let start_line = panels.remove(0);
    let start_bar = panels.remove(0);
    let end_line = panels.remove(0);
    let end_bar = panels.remove(0);
    let panels = vec![start_line, end_line, start_bar, end_bar];
    let start_line = layouts.remove(0);
    let start_bar = layouts.remove(0);
    let end_line = layouts.remove(0);
    let end_bar = layouts.remove(0);
    let layouts = vec![start_line, end_line, start_bar, end_bar];
    let scene = Figure::new(2, 2)
        .with_plots(panels)
        .with_layouts(layouts)
        .with_shared_x(0)
        .with_shared_x(1)
        .with_shared_y(0)
        .with_spacing(5.0)
        .with_row_height(0, 480.0)
        .with_row_height(1, 230.0)
        .with_figure_size(1200.0, 800.0)
        .render();
    save_scene(&output(prefix, "read_position"), scene)
}

fn run_if_present<F>(required: &[String], f: F) -> bool
where
    F: FnOnce() -> Result<()>,
{
    if required.iter().any(|p| !Path::new(p).is_file()) {
        return false;
    }
    f().is_ok()
}

pub fn generate(prefix: &str, options: &PlotOptions) {
    let file = |suffix: &str| format!("{prefix}.{suffix}");
    let report = file("summary_error_rate.csv");
    let spectrum = file("summary_error_spectrum.csv");
    let results = [
        run_if_present(&[spectrum.clone(), report.clone()], || {
            plot_spectrum(prefix, options.normalize)
        }),
        run_if_present(&[file("kvmer.csv"), report.clone()], || {
            plot_coverage(prefix)
        }),
        run_if_present(&[file("hazard_rate.csv"), report], || {
            plot_hazard(prefix, options.t_min, options.t_max, options.log_scale)
        }),
        run_if_present(&[file("summary_phred.csv")], || {
            plot_qscore(prefix, options.log_scale, options.min_coverage)
        }),
        run_if_present(&[file("summary_gc_content.csv")], || {
            plot_gc(prefix, options.log_scale, options.min_bases)
        }),
        run_if_present(&[file("summary_read_position.csv")], || {
            plot_read_position(prefix, options.num_bases)
        }),
    ];
    let written = results.iter().filter(|&&ok| ok).count();
    info!(
        "Plotting complete: wrote {} plot(s), skipped {}.",
        written,
        results.len() - written
    );
}