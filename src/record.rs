use crate::newgt::{NewGt, rewrite_gt_into};
use crate::target::{Target, matches_target};

/// Rewrite the FORMAT/GT column of one VCF data line.
///
/// For `Target::Query`, `sample_pass[i]` must be `true` iff sample `i` should be rewritten.
/// For other targets, `sample_pass` is ignored.
pub fn rewrite_record_into(
    line: &str,
    target: Target,
    new_gt: &NewGt,
    out: &mut String,
    sample_pass: Option<&[bool]>,
) {
    let bytes = line.as_bytes();

    let mut tab_positions = [0usize; 10];
    let mut n_tabs = 0usize;
    for (i, &b) in bytes.iter().enumerate() {
        if b == b'\t' {
            tab_positions[n_tabs] = i;
            n_tabs += 1;
            if n_tabs == 9 {
                break;
            }
        }
    }

    if n_tabs < 8 {
        out.push_str(line);
        return;
    }

    let format_start = tab_positions[7] + 1;
    let format_end = if n_tabs >= 8 {
        tab_positions[8]
    } else {
        bytes.len()
    };
    let format_col = &line[format_start..format_end];

    let gt_idx = format_col
        .split(':')
        .position(|k| k == "GT")
        .unwrap_or(usize::MAX);

    if gt_idx == usize::MAX {
        out.push_str(line);
        return;
    }

    out.push_str(&line[..format_end]);

    let mut cursor = format_end;
    let mut sample_idx = 0usize;
    while cursor < line.len() {
        let sample_start = cursor + 1;
        let sample_end = line[sample_start..]
            .find('\t')
            .map_or(line.len(), |p| sample_start + p);
        let sample_col = &line[sample_start..sample_end];
        out.push('\t');

        let mut field_i = 0usize;
        let mut gt_field_start = 0usize;
        let mut gt_field_end = sample_col.len();
        let sbytes = sample_col.as_bytes();
        let mut prev = 0usize;
        let mut found = false;

        for (j, &b) in sbytes.iter().enumerate() {
            if b == b':' {
                if field_i == gt_idx {
                    gt_field_start = prev;
                    gt_field_end = j;
                    found = true;
                    break;
                }
                field_i += 1;
                prev = j + 1;
            }
        }
        if !found && field_i == gt_idx {
            gt_field_start = prev;
            gt_field_end = sbytes
                .iter()
                .skip(prev)
                .position(|&b| b == b':')
                .map_or(sbytes.len(), |p| p + prev);
        }

        let gt_field = &sample_col[gt_field_start..gt_field_end];

        let should_rewrite = match target {
            Target::Query => sample_pass
                .and_then(|sp| sp.get(sample_idx))
                .copied()
                .unwrap_or(false),
            other => matches_target(gt_field, other),
        };

        if should_rewrite {
            out.push_str(&sample_col[..gt_field_start]);
            rewrite_gt_into(gt_field, new_gt, out);
            out.push_str(&sample_col[gt_field_end..]);
        } else {
            out.push_str(sample_col);
        }

        cursor = sample_end;
        sample_idx += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_line(gt1: &str, gt2: &str) -> String {
        format!("chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t{gt1}\t{gt2}")
    }

    #[test]
    fn rewrite_all_missing() {
        let line = make_line("./.", "0/1");
        let mut out = String::new();
        rewrite_record_into(&line, Target::All, &NewGt::Ref, &mut out, None);
        assert!(out.contains("\t0/0\t0/0"));
    }

    #[test]
    fn rewrite_any_missing_only() {
        let line = make_line("./.", "0/1");
        let mut out = String::new();
        rewrite_record_into(&line, Target::AnyMissing, &NewGt::Ref, &mut out, None);
        assert!(out.contains("\t0/0\t0/1"));
    }

    #[test]
    fn no_gt_field_passthrough() {
        let line = "chr1\t100\t.\tA\tT\t.\t.\t.\tDP\t30\t20";
        let mut out = String::new();
        rewrite_record_into(line, Target::All, &NewGt::Missing, &mut out, None);
        assert_eq!(out, line);
    }

    #[test]
    fn too_few_tabs_passthrough() {
        let line = "chr1\t100\t.\tA\tT\t.\t.\t.";
        let mut out = String::new();
        rewrite_record_into(line, Target::All, &NewGt::Missing, &mut out, None);
        assert_eq!(out, line);
    }

    #[test]
    fn query_target_uses_sample_pass() {
        let line = make_line("0/1", "0/1");
        let mut out = String::new();
        let pass = [true, false];
        rewrite_record_into(&line, Target::Query, &NewGt::Missing, &mut out, Some(&pass));
        assert!(out.contains("\t./.\t0/1"));
    }
}
