use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};
use rsomics_vcf_expr::EvalContext;

use crate::newgt::NewGt;
use crate::record::rewrite_record_into;
use crate::target::Target;

pub struct SetGtConfig {
    pub target: Target,
    pub new_gt: NewGt,
    /// Required when `target == Target::Query`; `None` otherwise.
    pub filter: Option<EvalContext>,
}

pub struct SetGtStats {
    pub total: u64,
    pub changed: u64,
}

/// Stream the VCF, rewriting GT fields. Single-threaded, constant RSS.
pub fn set_gt(input: &Path, output: &mut dyn Write, cfg: &SetGtConfig) -> Result<SetGtStats> {
    let file = std::fs::File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;

    let is_gz = {
        use std::io::Read as _;
        let mut buf = [0u8; 2];
        let mut peek = BufReader::new(&file);
        let n = peek.read(&mut buf).map_err(RsomicsError::Io)?;
        n >= 2 && buf[0] == 0x1f && buf[1] == 0x8b
    };

    let file2 = std::fs::File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;

    let mut writer = BufWriter::new(output);

    if is_gz {
        let decoder = flate2::read::MultiGzDecoder::new(file2);
        stream_lines(BufReader::new(decoder), &mut writer, cfg)
    } else {
        stream_lines(BufReader::new(file2), &mut writer, cfg)
    }
}

fn stream_lines<R: Read, W: Write>(
    reader: BufReader<R>,
    writer: &mut W,
    cfg: &SetGtConfig,
) -> Result<SetGtStats> {
    let mut total = 0u64;
    let mut changed = 0u64;
    let mut out_buf = String::with_capacity(4096);

    for raw in reader.lines() {
        let line = raw.map_err(RsomicsError::Io)?;

        if line.starts_with('#') {
            writer
                .write_all(line.as_bytes())
                .map_err(RsomicsError::Io)?;
            writer.write_all(b"\n").map_err(RsomicsError::Io)?;
            continue;
        }
        if line.is_empty() {
            continue;
        }

        total += 1;
        out_buf.clear();

        let sample_pass: Option<Vec<bool>> = if cfg.target == Target::Query {
            let ctx = cfg.filter.as_ref().ok_or_else(|| {
                RsomicsError::InvalidInput("-t q requires -i or -e expression".into())
            })?;
            let n_samples = count_samples(&line);
            let result = ctx
                .eval_line(&line, n_samples)
                .map_err(|e| RsomicsError::InvalidInput(format!("filter expression error: {e}")))?;
            Some(result.pass)
        } else {
            None
        };

        rewrite_record_into(
            &line,
            cfg.target,
            &cfg.new_gt,
            &mut out_buf,
            sample_pass.as_deref(),
        );

        if out_buf != line {
            changed += 1;
        }

        writer
            .write_all(out_buf.as_bytes())
            .map_err(RsomicsError::Io)?;
        writer.write_all(b"\n").map_err(RsomicsError::Io)?;
    }

    Ok(SetGtStats { total, changed })
}

/// Columns after FORMAT (col 8); VCF has 9 fixed cols (0-8).
fn count_samples(line: &str) -> usize {
    let mut tab_count = 0usize;
    for &b in line.as_bytes() {
        if b == b'\t' {
            tab_count += 1;
        }
    }
    tab_count.saturating_sub(8)
}
