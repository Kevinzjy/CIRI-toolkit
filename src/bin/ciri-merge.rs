//! Thin binary wrapper for the multi-sample `ciri-merge` command.

fn main() -> anyhow::Result<()> {
    ciri_toolkit::cli::merge::main()
}
