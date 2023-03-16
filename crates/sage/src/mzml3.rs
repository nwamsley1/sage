//CC=gcc pyenv install 3.11.0
//CXX=/usr/local/bin/g++-12 cargo run --release modmzml.json
//pyenv activate venv_rust_pyo3  
use crate::mzml::Spectrum;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyList, PyLong, PyFloat, PyUnicode};

use crate::mass::Tolerance;
use crate::spectrum::Precursor;

use std::fs::File;
use std::collections::HashMap;
use std::fs::File;
use std::collections::HashMap;
use arrow2::array::{Array, PrimitiveArray, ListArray, Utf8Array};
use arrow2::chunk::Chunk;
use arrow2::datatypes::Schema;
use arrow2::error::Result;
use arrow2::io::ipc::{read::{FileReader, read_file_metadata}};

//use std::collections::HashMap;
use std::{time::{//Duration, 
          Instant
        }};

//#[derive(Default, Debug, Clone)]
//pub struct Spectrum {
//    pub ms_level: u8,
//    pub id: String,
//    // pub scan_id: Option<usize>,
//    pub precursors: Vec<Precursor>,
//    /// Profile or Centroided data
//    pub representation: Representation,
//    /// Scan start time
//    pub scan_start_time: f32,
//   /// Ion injection time
//  pub ion_injection_time: f32,
    /// Total ion current
//    pub total_ion_current: f32,
    /// M/z array
//    pub mz: Vec<f32>,
    /// Intensity array
//    pub intensity: Vec<f32>,
//}

//#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
//pub enum Representation {
//    Profile,
//    Centroid,
//}

//impl Default for Representation {
//    fn default() -> Self {
//        Self::Profile
//    }
//}

#[derive(Debug)]
pub enum RawFileError {
    Malformed,
    IOError(std::io::Error),
}

impl std::error::Error for RawFileError {}

impl From<std::io::Error> for RawFileError {
    fn from(residual: std::io::Error) -> Self {
        Self::IOError(residual)
    }
}

impl From<std::str::Utf8Error> for RawFileError {
    fn from(_: std::str::Utf8Error) -> Self {
        Self::Malformed
    }
}

impl std::fmt::Display for RawFileError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RawFileError::Malformed => f.write_str("MzMLError: malformed cvParam"),
            RawFileError::IOError(s) => write!(f, "MzMLError: IO error {}", s),
        }
    }
}

#[derive(Default)]
pub struct RawFileReader {
    ms_level: Option<u8>,
}

impl RawFileReader {

    /// Create a new [`RawFileReader`] with a minimum MS level filter
    ///
    /// # Example
    ///
    /// A minimum level of 2 will not parse or return MS1 scans
    pub fn with_level_filter(ms_level: u8) -> Self {
        Self {
            ms_level: Some(ms_level),
        }
    }

    //#[derive(Debug)]
    //pub enum RawFileError {
    //    Malformed
   // }
   // impl std::error::Error for RawFileError {}

    //impl From<std::str::Utf8Error> for RawFileError {
    //    fn from(_: std::str::Utf8Error) -> Self {
    //        Self::Malformed
    //    }
   // }
    pub fn parse(filename: &String) -> Result<Vec<Spectrum>, RawFileError>{
        println!("{:?}", filename);
        let (schema, scans) = read_chunks(&path)?;
        scans = scans;
        let TIC = scans[key_value["TIC"]].as_any().downcast_ref::<PrimitiveArray<f32>>().unwrap();
        let scanType = scans[key_value["scanType"]].as_any().downcast_ref::<Utf8Array<i32>>().unwrap();
        let msOrder = scans[key_value["msOrder"]].as_any().downcast_ref::<PrimitiveArray<i32>>().unwrap();
        let retentionTime = scans[key_value["retentionTime"]].as_any().downcast_ref::<PrimitiveArray<f32>>().unwrap();
        let scanNumber = scans[key_value["scanNumber"]].as_any().downcast_ref::<PrimitiveArray<i32>>().unwrap();
        let scanMasses = scans[key_value["masses"]].as_any().downcast_ref::<ListArray<i32>>().unwrap().as_any().downcast_ref::<ListArray<i32>>().unwrap();
        let scanIntensities = scans[key_value["intensities"]].as_any().downcast_ref::<ListArray<i32>>().unwrap().as_any().downcast_ref::<ListArray<i32>>().unwrap();
        let precursorMz = scans[key_value["precursorMZ"]].as_any().downcast_ref::<PrimitiveArray<f32>>().unwrap();
        let TIC = scans[key_value["TIC"]].as_any().downcast_ref::<PrimitiveArray<f32>>().unwrap();
        
        for scan_id in 0..scans.len(){
            let mut spectrum: Spectrum = Spectrum::default();
            //Apply any filters here
            if scanType.value(scan_id).contains("ITMS"){
                            //println!("skipping");
                return spectrum
            } 
            let mut precursor = Precursor::default();
            spectrum.ms_level = msOrder.value(scan_id);
            spectrum.mz = scanMasses.value(scan_id);
            spectrum.intensity = scanIntensities.value(scan_id);
            spectrum.scan_start_time = retentionTime.value(scan_id);
            spectrum.total_ion_current = TIC.value(scan_id);
            spectrum.id = scan_id;
            if spectrum.ms_level>1 {

                precursor.mz  = precursorMz.value(scan_id);
                let isolation_width: f32 = 1.0; //Need to edit my arrow files to include this. 
                //println!("{:?}",precursor.mz - isolation_width/2.0);
                //if precursor.mz != 0.0 {
                //    precursor.isolation_window = match (precursor.mz - isolation_width/2.0, precursor.mz + isolation_width/2.0) {
                //        (Some(lo), Some(hi)) => Some(Tolerance::Da(-lo, hi)),
                //        _ => None,
                //    };
                //    spectrum.precursors.push(precursor);
                //    precursor = Precursor::default();
                // }
                precursor.isolation_window = Some(Tolerance::Da(-(precursor.mz - isolation_width/2.0), 
                                                            precursor.mz + isolation_width/2.0)
                                                );

                spectrum.precursors.push(precursor);
                precursor = Precursor::default();
            }
            spectra.push(spectrum);
        }
        //Open a handle to the arrow file. Needs to be in the same directory as main.rs
        fn read_chunks(filename: &str) -> Result<(HashMap<String, usize>, Vec<Chunk<Box<dyn Array>>>)> {
            println!("Hello, world!");
            let mut reader = File::open(&path).unwrap();
            let metadata = read_file_metadata(&mut reader).unwrap();
        
            let schema = metadata.schema.clone();
            
            let mut reader = FileReader::new(reader, metadata, None, None);
        
            let chunks = reader.collect::<Result<Vec<_>>>().unwrap();
            let names = schema.fields.iter().map(|f| &f.name).collect::<Vec<_>>();
            let mut key_value = HashMap::new();
            for name in names.iter().enumerate() {
               // println!("{}", (name.0, name.1));
               key_value.insert(name.1.to_string(), name.0);
            }

            Ok((key_value, chunks))
        }

        match true {
            false => Err(RawFileError::Malformed),
            true => Ok(spectra),
        }
}
}