//CC=gcc pyenv install 3.11.0
//CXX=/usr/local/bin/g++-12 cargo run --release modmzml.json
//pyenv activate venv_rust_pyo3  
use crate::mzml::Spectrum;
use crate::mzml::Representation;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyList, PyLong, PyFloat, PyUnicode};
use std::iter::Iterator;

use crate::mass::Tolerance;
use crate::spectrum::Precursor;

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

    pub fn parse() -> Result<Vec<Spectrum>, RawFileError>{
        pyo3::prepare_freethreaded_python();
        //env::set_var("RUST_BACKTRACE", "full");
        let spectra: Vec<Spectrum> = Python::with_gil(|py| {

            //Load python module that uses fisher_py to load the 
            //Thermo.CommonCore Dll's into python
            println!("yikes!");
            //let _raw_handle: &PyAny = get_raw_handle("HELA_uPAC_200cm_20221211_04.raw", py);
            
            //et _raw_handle: &PyAny = get_raw_handle("MA4365_FFPE_HPVpos_08_071522.raw", py);
            let _raw_handle: &PyAny = get_raw_handle("MA4427_16_WT_082722.raw", py);
            println!("mistake was not on line 88");
            //Get first and last scan numbers
            let first_scan_number: u32 = _raw_handle.getattr("run_header_ex")
            .unwrap().getattr("first_spectrum").unwrap().downcast::<PyLong>().unwrap().extract().unwrap();
            let last_scan_number: u32 = _raw_handle.getattr("run_header_ex")
            .unwrap().getattr("last_spectrum").unwrap().downcast::<PyLong>().unwrap().extract().unwrap();
            //let first_scan_number: u32 = 29990;
            //let last_scan_number: u32 = 30000;

            let now = Instant::now();
            let mut count = first_scan_number;
            let mut spectra = Vec::new();
            //Load all spectra into a Vec<Spectrum>

            loop {
                let spectrum: Spectrum = get_centroid_stream(count, _raw_handle);
                //println!("{:?}",spectrum.0.ion_injection_time);
                //println!("Hellow");
                spectra.push(spectrum);
                //if spectrum.1 > 0{}
                if count == last_scan_number {
                    break;
                }
                count += 1;
            }
            let new_now = Instant::now();
            println!("{:?}", new_now.saturating_duration_since(now));
            println!("{:?}", now.saturating_duration_since(new_now));
            //println!("{:?}", spectra[29000].0.get_heap_size());
            //println!("{:?}", spectra.get_heap_size());
            //check type of output
            fn print_type_of<T>(_: &T) {
                            println!("{}", std::any::type_name::<T>())
            }
            spectra
        });

    fn get_centroid_stream<'a>(scan_identifier: u32, _raw_handle : &'a PyAny) -> Spectrum {
        //Default spectrum and precursor
        let mut spectrum: Spectrum = Spectrum::default();
        let mut precursor = Precursor::default();

        //scan header string. Useful for filtering scans 
        let scan_filter: String = _raw_handle.getattr("get_scan_event_string_for_scan_number").unwrap().call1((scan_identifier,)).unwrap().downcast::<PyUnicode>().unwrap().extract().unwrap();
    
        //Apply any filters here
        if scan_filter.contains("ITMS"){
            //println!("skipping");
           return spectrum
        } 

        //Data need to parse the spectrum comes from "scan_event", "scan_stats" and "centroid_stream"
        //Can apply functions that filter based on "scan event"
        let scan_event: &PyAny = _raw_handle.getattr("get_scan_event_for_scan_number").unwrap().call1((scan_identifier,)).unwrap();
        let scan_stats: &PyAny = _raw_handle.getattr("get_scan_stats_for_scan_number")
        .unwrap().call1((scan_identifier,)).unwrap(); 
        let centroid_stream: &PyAny = _raw_handle.getattr("get_centroid_stream").unwrap().call1((scan_identifier, false)).unwrap();

        //println!("{:?}", scan_filter);
        spectrum.ms_level = scan_event.getattr("ms_order").unwrap()
        .getattr("value").unwrap().downcast::<PyLong>().unwrap().extract().unwrap();


        // Only get the spectrum if it is larger than zero. Otherwise return the default empty spectrum
        let spectrum_size: u32 = centroid_stream.getattr("length")
                                                .unwrap().downcast::<PyLong>().unwrap().extract().unwrap();  

        if spectrum_size > 0 {

            spectrum.mz = centroid_stream.getattr("masses")
            .unwrap().downcast::<PyList>().unwrap().extract().unwrap();

            spectrum.intensity = centroid_stream.getattr("intensities")
            .unwrap().downcast::<PyList>().unwrap().extract().unwrap();

            spectrum.scan_start_time = scan_stats.getattr("start_time").unwrap().downcast::<PyFloat>().unwrap().extract().unwrap();
            spectrum.total_ion_current = scan_stats.getattr("tic").unwrap().downcast::<PyFloat>().unwrap().extract().unwrap();
            spectrum.id = scan_stats.getattr("scan_number").unwrap().downcast::<PyLong>().unwrap().extract::<u32>().unwrap().to_string();
            //.extract().unwrap();
            //spectrum.ion_injection_time = scan_stats.get("Ion Injection Time (ms):").as_ref().unwrap().parse::<f32>().unwrap();
            //println!("{:?}", spectrum.ion_injection_time);
            //println!("{:?}", spectrum.ms_level);

            //Make precursor
            if spectrum.ms_level>1 {

                precursor.mz  = scan_event.getattr("get_mass").unwrap().call1((0, )).unwrap().downcast::<PyFloat>().unwrap().extract().unwrap();
                let isolation_width: f32 = scan_event.getattr("get_isolation_width").unwrap().call1((0, )).unwrap().downcast::<PyFloat>().unwrap().extract().unwrap();
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

        //precursor iosolation width
        //scan_event.getattr("get_isolation_width").unwrap().call1((0, )).unwrap().downcast::<PyFloat>().unwrap().extract().unwrap()
        
            spectrum.scan_start_time = _raw_handle.getattr("retention_time_from_scan_number").unwrap().call1((scan_identifier,)).unwrap()
            .downcast::<PyFloat>().unwrap().extract().unwrap();//::<PyLong>().unwrap();//::<PyLong>().unwrap().extract().unwrap();

            //spectrum.scan_start_time = scan_event.getattr( "")
            return spectrum
            } else {
            //println!("spectrum size was zero");
            return spectrum
        }
    }

    //Open a handle to the raw file. Needs to be in the same directory as main.rs
    //fn get_raw_handle<'a>(filename: &'a str, py: Python<'a>) -> Result<&'a PyAny, PyErr> {
    fn get_raw_handle<'a>(filename: &'a str, py: Python<'a>) -> &'a PyAny {

        //hardcode location of load_raw. This module uses fisher_py
        //to load the Thermo.CommonCore DLLs and open a handle to the raw file
        fn get_load_raw_module(py: Python) -> &PyModule {

            let code = std::fs::read_to_string("/Users/n.t.wamsley/Projects/SAGE_TESTING/sage/crates/sage/src/load_raw.py").unwrap();
            PyModule::from_code(py, &code,"load_raw","load_raw").expect("fail")

        }

        //load handle to the raw file
        get_load_raw_module(py).getattr("load_raw").unwrap()
                               .call1((filename,)).unwrap()
    }
    match true {
        false => Err(RawFileError::Malformed),
        true => Ok(spectra),
    }
}
}