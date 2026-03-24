//! Runtime diagnostics configuration shared by Scan1/Scan2/validator code.
//!
//! The hot paths still expose the historical `CIRI_TRACE_*` / `CIRI_PROFILE_*`
//! environment-variable behavior for ad hoc debugging, but the preferred entry
//! point is now the CLI. This module keeps that wiring in one place so the
//! algorithm code can query "should I trace/profile?" without knowing whether
//! the request came from an environment variable or an explicit command-line
//! flag.

use anyhow::Result;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::{Mutex, OnceLock};

/// Process-wide runtime diagnostics options initialized from the CLI.
struct RuntimeConfig {
    debug_reads: HashSet<String>,
    debug_writer: Option<Mutex<BufWriter<File>>>,
    perf_writer: Option<Mutex<BufWriter<File>>>,
}

static RUNTIME_CONFIG: OnceLock<RuntimeConfig> = OnceLock::new();

thread_local! {
    /// Tracks whether the current thread is inside a traced `is_bsj_hg2` call.
    ///
    /// This keeps `--debug` focused on the requested read IDs instead of turning
    /// on verbose validator branch logging for every candidate in the process.
    static TRACE_HG2_ACTIVE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Initializes the optional debug/perf outputs for the current process.
///
/// This should be called once near startup. Later calls are ignored so tests
/// and library users do not need to coordinate teardown.
pub fn init_runtime(
    debug_reads: Option<&str>,
    debug_log_path: Option<&str>,
    perf_path: Option<&str>,
) -> Result<()> {
    let reads = debug_reads
        .map(|raw| {
            raw.split(',')
                .filter_map(|token| {
                    let t = token.trim();
                    (!t.is_empty()).then(|| t.to_string())
                })
                .collect::<HashSet<_>>()
        })
        .unwrap_or_default();
    let debug_writer = match (debug_log_path, reads.is_empty()) {
        (Some(path), false) => Some(Mutex::new(BufWriter::new(File::create(path)?))),
        _ => None,
    };
    let perf_writer = match perf_path {
        Some(path) => Some(Mutex::new(BufWriter::new(File::create(path)?))),
        None => None,
    };
    let _ = RUNTIME_CONFIG.set(RuntimeConfig {
        debug_reads: reads,
        debug_writer,
        perf_writer,
    });
    Ok(())
}

/// Returns whether the given read ID is selected for targeted tracing.
#[inline]
pub fn should_trace_read(read_id: &str) -> bool {
    if let Some(config) = RUNTIME_CONFIG.get() {
        if config.debug_reads.contains(read_id) {
            return true;
        }
    }
    if let Ok(raw) = std::env::var("CIRI_TRACE_READS") {
        for token in raw.split(',') {
            let t = token.trim();
            if !t.is_empty() && t == read_id {
                return true;
            }
        }
    }
    false
}

/// Returns whether any CLI-provided debug reads are active.
#[inline]
pub fn cli_debug_enabled() -> bool {
    RUNTIME_CONFIG
        .get()
        .is_some_and(|config| !config.debug_reads.is_empty())
}

/// Executes a closure while enabling `is_bsj_hg2` branch tracing on this thread.
#[inline]
pub fn with_trace_hg2_scope<T>(enabled: bool, f: impl FnOnce() -> T) -> T {
    let prev = TRACE_HG2_ACTIVE.with(|flag| {
        let prev = flag.get();
        flag.set(enabled);
        prev
    });
    let out = f();
    TRACE_HG2_ACTIVE.with(|flag| flag.set(prev));
    out
}

/// Returns whether `is_bsj_hg2` branch-level trace output is active.
#[inline]
pub fn trace_hg2_enabled() -> bool {
    TRACE_HG2_ACTIVE.with(|flag| flag.get())
        || matches!(std::env::var("CIRI_TRACE_HG2"), Ok(v) if !v.is_empty() && v != "0")
}

/// Returns whether Scan1 profiling is enabled.
#[inline]
pub fn scan1_profile_enabled() -> bool {
    RUNTIME_CONFIG
        .get()
        .is_some_and(|config| config.perf_writer.is_some())
        || matches!(std::env::var("CIRI_PROFILE_SCAN1"), Ok(v) if !v.is_empty() && v != "0")
}

/// Returns whether Scan2 profiling is enabled.
#[inline]
pub fn scan2_profile_enabled() -> bool {
    RUNTIME_CONFIG
        .get()
        .is_some_and(|config| config.perf_writer.is_some())
        || matches!(std::env::var("CIRI_PROFILE_SCAN2"), Ok(v) if !v.is_empty() && v != "0")
}

/// Emits one debug trace line either to the configured file or to stderr.
pub fn emit_debug_line(line: &str) {
    if let Some(config) = RUNTIME_CONFIG.get() {
        if let Some(writer) = &config.debug_writer {
            let mut writer = writer.lock().unwrap();
            let _ = writeln!(writer, "{}", line);
            let _ = writer.flush();
            return;
        }
    }
    eprintln!("{}", line);
}

/// Emits one profiling line either to the configured file or to stderr.
pub fn emit_perf_line(line: &str) {
    if let Some(config) = RUNTIME_CONFIG.get() {
        if let Some(writer) = &config.perf_writer {
            let mut writer = writer.lock().unwrap();
            let _ = writeln!(writer, "{}", line);
            let _ = writer.flush();
            return;
        }
    }
    eprintln!("{}", line);
}
