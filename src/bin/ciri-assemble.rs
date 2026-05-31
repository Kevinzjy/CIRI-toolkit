//! Thin binary wrapper for the multi-sample `ciri-assemble` command.

fn main() -> anyhow::Result<()> {
    ciri_toolkit::cli::assemble::main()
}
