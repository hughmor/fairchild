//! Which output column a probe name selects.
//!
//! `--probe` is honoured in two mechanisms, because the analyses differ in
//! kind: a bounded result (`.ac`, `.dc`, `.op`) is filtered after it is
//! rendered, and a transient is filtered before it is, since rendering the
//! whole run is the thing streaming exists to avoid. Two mechanisms are fine.
//! Two *rules* would not be — a user whose `V(OUT)` selected a column in one
//! analysis and not in another would have no way to tell which. So the rule
//! lives here, once, and both mechanisms ask it.

/// Does `probe` name `column`?
///
/// Matching is case-insensitive. `V(out)` also selects an AC sweep's
/// `mag_V(out)` / `phase_deg_V(out)` pair, since those are the columns that
/// signal produces there.
pub fn matches(probe: &str, column: &str) -> bool {
    let p = probe.to_lowercase();
    let c = column.to_lowercase();
    c == p
        || c.strip_prefix("mag_").is_some_and(|h| h == p)
        || c.strip_prefix("phase_deg_").is_some_and(|h| h == p)
}

/// The probes in `probes` that name no column of `columns`.
///
/// A probe that matches nothing is an error rather than a silent drop (#72):
/// a result that looks complete and is not becomes a `KeyError` a long way
/// downstream from its cause.
pub fn unmatched(probes: &[String], columns: &[String]) -> Vec<String> {
    probes
        .iter()
        .filter(|p| !columns.iter().any(|c| matches(p, c)))
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn case_and_ac_pair() {
        assert!(matches("V(OUT)", "v(out)"));
        assert!(matches("v(out)", "mag_V(out)"));
        assert!(matches("v(out)", "phase_deg_V(out)"));
        assert!(!matches("v(out)", "v(out2)"));
        assert!(!matches("v(out)", "i(out)"));
    }

    #[test]
    fn unmatched_names_come_back() {
        let cols = vec!["time".to_string(), "V(out)".to_string()];
        let probes = vec!["V(out)".to_string(), "V(nope)".to_string()];
        assert_eq!(unmatched(&probes, &cols), vec!["V(nope)".to_string()]);
    }
}
