//! Thin binary wrapper for the development-only `ciri-simulator` command.

fn main() -> anyhow::Result<()> {
    ciri_toolkit::cli::simulator::main()
}
