//CC=gcc pyenv install 3.11.0
//CXX=/usr/local/bin/g++-12 cargo run --release modmzml.json
//pyenv activate venv_rust_pyo3  
use crate::mzml::Spectrum;
use crate::mass::Tolerance;
use crate::spectrum::Precursor;

use std::fs::File;
use std::collections::HashMap;
use arrow2::array::{Array, PrimitiveArray, ListArray, Utf8Array, Float32Array};
use arrow2::chunk::Chunk;
use arrow2::datatypes::Schema;
use arrow2::error::Result;
use arrow2::io::ipc::{read::{FileReader, read_file_metadata}};
use arrow2::error::Error;
//use std::collections::HashMap;
use std::{time::{//Duration, 
          Instant
        }};


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
    //    } parse<B: AsyncBufRead + Unpin>(&self, b: B) -> Result<Vec<Spectrum>, MzMLError>
   // }
    pub fn parse(filename: &String) -> Vec<Spectrum>{
        println!("{:?}", filename);
        let mut spectra: Vec<Spectrum> = Vec::new();
        let (schema, scans) = read_chunks(&filename).unwrap();
        //scans = scans[0];
        println!("{}",filename);
        let TIC = scans[0][schema["TIC"]].as_any().downcast_ref::<PrimitiveArray<f32>>().unwrap();
        let scanType = scans[0][schema["scanType"]].as_any().downcast_ref::<Utf8Array<i32>>().unwrap();
        let msOrder = scans[0][schema["msOrder"]].as_any().downcast_ref::<PrimitiveArray<i32>>().unwrap();
        let retentionTime = scans[0][schema["retentionTime"]].as_any().downcast_ref::<PrimitiveArray<f32>>().unwrap();
        let scanNumber = scans[0][schema["scanNumber"]].as_any().downcast_ref::<PrimitiveArray<i32>>().unwrap();
        let scanMasses = scans[0][schema["masses"]].as_any().downcast_ref::<ListArray<i32>>().unwrap().as_any().downcast_ref::<ListArray<i32>>().unwrap();
        let scanIntensities = scans[0][schema["intensities"]].as_any().downcast_ref::<ListArray<i32>>().unwrap().as_any().downcast_ref::<ListArray<i32>>().unwrap();
        let precursorMz = scans[0][schema["precursorMZ"]].as_any().downcast_ref::<PrimitiveArray<f32>>().unwrap();
        let TIC = scans[0][schema["TIC"]].as_any().downcast_ref::<PrimitiveArray<f32>>().unwrap();
        
        for scan_id in 0..scans[0].len(){
            let mut spectrum: Spectrum = Spectrum::default();
            //Apply any filters here
            if scanType.value(scan_id).contains("ITMS"){
                            //println!("skipping");
                continue;
            } 
            let mut precursor = Precursor::default();
            spectrum.ms_level = msOrder.value(scan_id) as u8;
            spectrum.mz = scanMasses.value(scan_id).as_any().downcast_ref::<Float32Array>().unwrap().clone().values().to_vec();
            spectrum.intensity = scanIntensities.value(scan_id).as_any().downcast_ref::<Float32Array>().unwrap().clone().values().to_vec();
            spectrum.scan_start_time = retentionTime.value(scan_id);
            spectrum.total_ion_current = TIC.value(scan_id);
            spectrum.id = scanNumber.value(scan_id).to_string();
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
            } else{
                spectrum.precursors.push(precursor);
            }
            spectra.push(spectrum);
        }
        //Open a handle to the arrow file. Needs to be in the same directory as main.rs
        fn read_chunks(filename: &str) -> Result<(HashMap<String, usize>, Vec<Chunk<Box<dyn Array>>>)> {
            println!("Hello, world!");
            let mut reader = File::open(&filename).unwrap();
            let metadata = read_file_metadata(&mut reader).unwrap();
        
            let schema = metadata.schema.clone();
            
            let reader = FileReader::new(reader, metadata, None, None);
        
            let chunks = reader.collect::<Result<Vec<_>>>().unwrap();
            let names = schema.fields.iter().map(|f| &f.name).collect::<Vec<_>>();
            let mut key_value = HashMap::new();
            for name in names.iter().enumerate() {
               // println!("{}", (name.0, name.1));
               key_value.insert(name.1.to_string(), name.0);
            }

            Ok((key_value, chunks))
        }

        //match true {
            //false => Err(Error::Overflow),//Err(RawFileError::Malformed),
            //false => Err(arrow2::error::Error::OutOfSpec),
        //    true => Ok(spectra),
       // }
       //println!("{:?}", spectra[0]);
       println!("number of spectra {}", spectra.len());
       return spectra
}
}