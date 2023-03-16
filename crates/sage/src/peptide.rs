use std::{collections::HashMap};

use itertools::Itertools;

use crate::{
    enzyme::{Digest, Position},
    mass::{Mass, Residue, H2O, VALID_AA},
};

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct Peptide {
    pub decoy: bool,
    pub sequence: Vec<Residue>,
    /// Modification on peptide C-terminus
    pub nterm: Option<f32>,
    /// Modification on peptide C-terminus
    pub cterm: Option<f32>,
    /// Monoisotopic mass, inclusive of N/C-terminal mods
    pub monoisotopic: f32,
    /// Number of missed cleavages for this sequence
    pub missed_cleavages: u8,
    /// Where is this peptide located in the protein?
    pub position: Position,
}

impl Peptide {
    fn set_nterm_mod(&mut self, m: f32) {
        if self.nterm.is_none() {
            self.nterm = Some(m);
            self.monoisotopic += m;
        }
    }

    fn set_cterm_mod(&mut self, m: f32) {
        if self.cterm.is_none() {
            self.cterm = Some(m);
            self.monoisotopic += m;
        }
    }

    fn set_nterm_protein_mod(&mut self, m: f32) {
        if self.position == Position::Full || self.position == Position::Nterm {
            self.set_nterm_mod(m)
        }
    }

    fn set_cterm_protein_mod(&mut self, m: f32) {
        if self.position == Position::Full || self.position == Position::Cterm {
            self.set_cterm_mod(m)
        }
    }

    pub fn label(&self) -> i32 {
        match self.decoy {
            true => -1,
            false => 1,
        }
    }

    /// Apply a static modification to a peptide in-place
    pub fn static_mod(&mut self, position_modifier: char, residue: char, mass: f32) {
        //
        if position_modifier == '^' {
            //If the residue specifier is generic "*" add the modification to the site
            if residue == '*'{
                self.set_nterm_mod(mass)
            //If the residue specifier is specific to an AA, only add the modification
            //To the position if the specified AA is present at the position
            } else {
                match *self.sequence.first().unwrap() == Residue::Just(residue) {
                    true => self.set_nterm_mod(mass),
                    false => {}
                }
            }

        } else if position_modifier == '$' {
            //If the residue specifier is generic "*" add the modification to the site
            if residue == '*'{
                self.set_cterm_mod(mass)
            //If the residue specifier is specific to an AA, only add the modification
            //To the position if the specified AA is present at the position
            } else {
                match *self.sequence.last().unwrap() == Residue::Just(residue) {
                    true => self.set_cterm_mod(mass),
                    false => {}
                }
            }

        } else if position_modifier == '[' {

            if residue == '*'{
                self.set_nterm_protein_mod(mass)
            //If the residue specifier is specific to an AA, only add the modification
            //To the position if the specified AA is present at the position
            } else {
                match *self.sequence.first().unwrap() == Residue::Just(residue) {
                    true => self.set_nterm_protein_mod(mass),
                    false => {}
                }
            }

        } else if position_modifier == ']' {
            if residue == '*'{
                self.set_cterm_protein_mod(mass)
            //If the residue specifier is specific to an AA, only add the modification
            //To the position if the specified AA is present at the position
            } else {
                match *self.sequence.last().unwrap() == Residue::Just(residue) {
                    true => self.set_cterm_protein_mod(mass),
                    false => {}
                }
            }

        } else {
            for resi in self.sequence.iter_mut() {
                // Don't overwrite an already modified amino acid!
                match resi {
                    Residue::Just(c) if *c == residue => {
                        self.monoisotopic += mass;
                        *resi = Residue::Mod(residue, mass);
                    }
                    _ => {}
                }
            }
        }
    }

    /// Apply all variable mods in `sites` to self
    fn apply_variable_mods(&mut self, sites: &[&(Site, f32)]) {
        for (site, mass) in sites {
            match site {
                Site::PeptideN => self.set_nterm_mod(*mass),
                Site::PeptideC => self.set_cterm_mod(*mass),
                Site::ProteinN => self.set_nterm_protein_mod(*mass),
                Site::ProteinC => self.set_cterm_protein_mod(*mass),
                Site::Sequence(index) => {
                    if let Residue::Just(c) = self.sequence[*index as usize] {
                        self.sequence[*index as usize] = Residue::Mod(c, *mass);
                        self.monoisotopic += mass;
                    }
                }
            }
        }
    }

    fn modification_sites(&self, residue: char, mass: f32) -> ModificationSites {
        ModificationSites {
            peptide: self,
            index: 0,
            residue,
            mass,
        }
    }

    /// Apply variable modifications, then static modifications to a peptide
    pub fn apply(
        mut self,
        variable_mods: &HashMap<String, f32>,//&[(String, f32)],
        static_mods: &HashMap<String, f32>,
        combinations: usize,
    ) -> Vec<Peptide> {

        //If only static mods
        if variable_mods.is_empty() {
            for (resi, mass) in static_mods {
                self.static_mod(resi.chars().nth(0).unwrap(), resi.chars().nth_back(0).unwrap(), *mass);
            }
            vec![self]
        } else {
            // Create list of all possible variable modifications
            let mods = variable_mods
                .iter()
                .fold(vec![], |mut acc, (residue, mass)| {
                    //acc.extend(self.modification_sites(residue.chars().nth_back(0).unwrap(), *mass));
                    acc.extend(self.modification_sites(residue.chars().nth_back(0).unwrap(), *mass));
                    acc
                });

            let mut modified = Vec::new();
            modified.push(self.clone());

            for n in 1..=combinations {
                for combination in mods.iter().combinations(n).filter(no_duplicates) {
                    let mut peptide = self.clone();
                    peptide.apply_variable_mods(&combination);
                    modified.push(peptide);
                }
            }

            // Apply static mods to all peptides
            for peptide in modified.iter_mut() {
                for (residue, &mass) in static_mods {
                    peptide.static_mod(residue.chars().nth(0).unwrap(), residue.chars().nth_back(0).unwrap(), mass);
                }
            }

            modified
        }
    }

    /// If `self` is a decoy peptide, un-reverse it
    pub fn pseudo_forward(&self) -> Option<Peptide> {
        if self.decoy {
            let mut fwd = self.clone();
            if fwd.sequence.len() > 2 {
                let n = fwd.sequence.len().saturating_sub(1);
                fwd.sequence[1..n].reverse();
            }
            return Some(fwd);
        }
        None
    }

    pub fn ambiguous(&self, other: &Peptide) -> bool {
        self.sequence.len() == other.sequence.len()
            && self.monoisotopic == other.monoisotopic
            && self
                .sequence
                .iter()
                .zip(other.sequence.iter())
                .all(|(l, r)| l.monoisotopic() == r.monoisotopic())
    }
}

fn no_duplicates(combination: &Vec<&(Site, f32)>) -> bool {
    let mut n = 0;
    let mut c = 0;
    for (site, _) in combination {
        match site {
            Site::PeptideN => n += 1,
            Site::PeptideC => c += 1,
            Site::ProteinN => n += 1,
            Site::ProteinC => c += 1,
            _ => {}
        }
    }

    n <= 1 && c <= 1
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Site {
    PeptideN,
    PeptideC,
    ProteinN,
    ProteinC,
    Sequence(u32),
}

struct ModificationSites<'a> {
    peptide: &'a Peptide,
    index: usize,
    residue: char,
    mass: f32,
}

impl<'a> Iterator for ModificationSites<'a> {
    type Item = (Site, f32);

    fn next(&mut self) -> Option<Self::Item> {
        match (self.residue, self.index, self.peptide.position) {
            ('^', 0, _) => {
                self.index = self.peptide.sequence.len();
                return Some((Site::PeptideN, self.mass));
            }
            ('[', 0, Position::Nterm) | ('[', 0, Position::Full) => {
                self.index = self.peptide.sequence.len();
                return Some((Site::ProteinN, self.mass));
            }
            ('$', 0, _) => {
                self.index = self.peptide.sequence.len();
                return Some((Site::PeptideC, self.mass));
            }
            (']', 0, Position::Cterm) | (']', 0, Position::Full) => {
                self.index = self.peptide.sequence.len();
                return Some((Site::ProteinC, self.mass));
            }
            _ => {}
        };
        while self.index < self.peptide.sequence.len() {
            let idx = self.index;
            self.index += 1;
            match self.peptide.sequence[idx] {
                Residue::Just(r) if r == self.residue => {
                    return Some((Site::Sequence(idx as u32), self.mass))
                }
                _ => continue,
            }
        }
        None
    }
}

impl TryFrom<&Digest> for Peptide {
    type Error = char;

    fn try_from(value: &Digest) -> Result<Self, Self::Error> {
        let mut sequence = Vec::with_capacity(value.sequence.len());
        let mut monoisotopic = H2O;

        for c in value.sequence.chars() {
            if !VALID_AA.contains(&c) {
                return Err(c);
            }
            monoisotopic += c.monoisotopic();
            sequence.push(Residue::Just(c));
        }

        Ok(Peptide {
            decoy: value.decoy,
            position: value.position,
            sequence,
            monoisotopic,
            nterm: None,
            cterm: None,
            missed_cleavages: value.missed_cleavages,
        })
    }
}

impl std::fmt::Display for Peptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(m) = self.nterm {
            if m.is_sign_positive() {
                write!(f, "[+{}]-", m)?;
            } else {
                write!(f, "[{}]-", m)?;
            }
        }
        f.write_str(
            &self
                .sequence
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(""),
        )?;
        if let Some(m) = self.cterm {
            if m.is_sign_positive() {
                write!(f, "-[+{}]", m)?;
            } else {
                write!(f, "-[{}]", m)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::enzyme::{Enzyme, EnzymeParameters};

    use super::*;

    fn var_mod_sequence(peptide: &Peptide, mods: &HashMap<String, f32>, combo: usize) -> Vec<String> {
        let static_mods = HashMap::default();
        peptide
            .clone()
            .apply(&mods, &static_mods, combo)
            .into_iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
    }

    #[test]
    fn test_variable_mods() {
        let variable_mods = [("*M".to_string(), 16.0f32), ("*C".to_string(), 57.)].into_iter().collect();
        let peptide = Peptide::try_from(&Digest {
            sequence: "GCMGCMG".into(),
            ..Default::default()
        })
        .unwrap();//.sort();//.apply(|x| x.to_owned()[..]);

        let expected = vec![
            "GCMGCMG",
            "GCM[+16]GCMG",
            "GCMGCM[+16]G",
            "GC[+57]MGCMG",
            "GCMGC[+57]MG",
            "GCM[+16]GCM[+16]G",
            "GC[+57]M[+16]GCMG",
            "GCM[+16]GC[+57]MG",
            "GC[+57]MGCM[+16]G",
            "GCMGC[+57]M[+16]G",
            "GC[+57]MGC[+57]MG",
        ].sort();//sort();//.apply(|x| x.to_owned()[..]);

        let peptides = var_mod_sequence(&peptide, &variable_mods, 2).sort();
        assert_eq!(peptides, expected);
    }

    #[test]
    fn test_variable_mods_no_effeect() {
        let variable_mods = [("*M".to_string(), 16.0f32), ("*C".to_string(), 57.)].into_iter().collect();
        let peptide = Peptide::try_from(&Digest {
            sequence: "AAAAAAAA".into(),
            ..Default::default()
        })
        .unwrap();

        let expected = vec!["AAAAAAAA"];
        let peptides = var_mod_sequence(&peptide, &variable_mods, 2);
        assert_eq!(peptides, expected);
    }

    /// Check that picked-peptide approach will match forward and reverse peptides
    #[test]
    fn test_psuedo_forward() {
        let trypsin = crate::enzyme::EnzymeParameters {
            missed_cleavages: 0,
            min_len: 3,
            max_len: 30,
            enyzme: Enzyme::new("KR", Some('P')),
        };

        let fwd = "MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGN";
        for digest in trypsin.digest(fwd) {
            let mut fwd = Peptide::try_from(&digest).unwrap();
            let mut rev = Peptide::try_from(&digest.reverse()).unwrap();

            assert_eq!(fwd.decoy, false);
            assert_eq!(rev.decoy, true);
            assert!(
                fwd.sequence.len() < 4 || fwd.sequence != rev.sequence,
                "{} {}",
                fwd,
                rev
            );
            assert_eq!(fwd.pseudo_forward(), None);
            assert_eq!(rev.pseudo_forward().unwrap().to_string(), fwd.to_string());

            fwd.static_mod('*', 'E', 15.0);
            rev.static_mod('*','E', 15.0);
            assert_eq!(rev.pseudo_forward().unwrap().to_string(), fwd.to_string());
        }
    }

    #[test]
    fn apply_mods() {
        let peptide = Peptide::try_from(&Digest {
            sequence: "AACAACAAK".into(),
            ..Default::default()
        })
        .unwrap();

        let mut expected = vec![
            "AAC[+57]AAC[+57]AAK",
            "AAC[+30]AAC[+57]AAK",
            "AAC[+57]AAC[+30]AAK",
            "AAC[+30]AAC[+30]AAK",
            "AAC[+57]AAC[+57]AAK[+5]",
            "AAC[+30]AAC[+57]AAK[+5]",
            "AAC[+57]AAC[+30]AAK[+5]",
            "AAC[+30]AAC[+30]AAK[+5]",
        ];//.sort();

        let mut static_mods = HashMap::new();
        static_mods.insert("*C".to_string(), 57.0);

        let variable_mods = [("*C".to_string(), 30.0),
                                                   ("$K".to_string(), 5.0)].into_iter().collect();

        let mut peptides = peptide
            .apply(&variable_mods, &static_mods, 3)
            .into_iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>();//.sort();
        assert_eq!(peptides.sort(), expected.sort());
    }

}
