//! Thin binary wrapper for the user-facing `ciri` command.

fn main() -> anyhow::Result<()> {
    ciri_toolkit::cli::ciri::main()
}
