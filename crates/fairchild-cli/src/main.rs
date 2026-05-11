use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use clap::{Parser, ValueEnum};

use fairchild_core::{ac_analysis, dc_op_nr, freq_decade, freq_linear, freq_oct, tran_nr_tr, DeviceRegistry};
use fairchild_parser::{parse_spice, AcVariation, Analysis};

#[derive(Parser)]
#[command(name = "fairchild", about = "SPICE-compatible analog circuit simulator")]
struct Cli {
    /// Input SPICE netlist file
    #[arg(short, long)]
    file: PathBuf,

    /// Output format
    #[arg(long, default_value = "csv")]
    format: Format,

    /// Output file (default: stdout)
    #[arg(short, long)]
    output: Option<PathBuf>,
}

#[derive(Clone, ValueEnum)]
enum Format {
    Csv,
    Nutmeg,
}

fn main() {
    let cli = Cli::parse();

    let src = fs::read_to_string(&cli.file).unwrap_or_else(|e| {
        eprintln!("error: cannot read {:?}: {e}", cli.file);
        std::process::exit(1);
    });

    let netlist = parse_spice(&src).unwrap_or_else(|e| {
        eprintln!("error: parse failed: {e}");
        std::process::exit(1);
    });

    let title = netlist.title.clone();

    let writer: Box<dyn Write> = match &cli.output {
        Some(path) => {
            let f = fs::File::create(path).unwrap_or_else(|e| {
                eprintln!("error: cannot create {:?}: {e}", path);
                std::process::exit(1);
            });
            Box::new(BufWriter::new(f))
        }
        None => Box::new(BufWriter::new(io::stdout())),
    };
    let mut w = writer;

    let mut ran_something = false;

    for analysis in &netlist.analyses {
        match analysis {
            Analysis::Op => {
                let result = dc_op_nr(&netlist).unwrap_or_else(|e| {
                    eprintln!("error: DC op failed: {e}");
                    std::process::exit(1);
                });
                match cli.format {
                    Format::Csv => result.write_csv(&mut w),
                    Format::Nutmeg => result.write_nutmeg(&mut w, &title),
                }
                .unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                ran_something = true;
            }
            Analysis::Tran { step, stop } => {
                let result = tran_nr_tr(&netlist, *step, *stop).unwrap_or_else(|e| {
                    eprintln!("error: tran failed: {e}");
                    std::process::exit(1);
                });
                match cli.format {
                    Format::Csv => result.write_csv(&mut w),
                    Format::Nutmeg => result.write_nutmeg(&mut w, &title),
                }
                .unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                ran_something = true;
            }
            Analysis::Ac { variation, points, fstart, fstop } => {
                let freqs = match variation {
                    AcVariation::Dec => freq_decade(*fstart, *fstop, *points),
                    AcVariation::Oct => freq_oct(*fstart, *fstop, *points),
                    AcVariation::Lin => freq_linear(*fstart, *fstop, *points),
                };
                let registry = DeviceRegistry::default();
                let result = ac_analysis(&netlist, &freqs, None, &registry).unwrap_or_else(|e| {
                    eprintln!("error: AC analysis failed: {e}");
                    std::process::exit(1);
                });
                match cli.format {
                    Format::Csv => result.write_csv(&mut w),
                    Format::Nutmeg => result.write_nutmeg(&mut w, &title),
                }
                .unwrap_or_else(|e| eprintln!("warning: write error: {e}"));
                ran_something = true;
            }
        }
    }

    if !ran_something {
        eprintln!("warning: no analyses found in netlist (add .op, .tran, or .ac)");
    }
}
