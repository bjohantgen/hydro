// This is the main file of the hydro simulation.
//
// Author: Brayden JoHantgen
// Last Update: 9/16/2026

////////////
// Imports
////////////
#![allow(unused_imports)]
use std::time::Instant;
use std::f64::consts::PI;

pub mod math_func;
pub mod sim_func;
pub mod test_conds;

///////////////
// Simulation
///////////////
fn main() {
    let sod = test_conds::sod_shock();
    let m_dot = math_func::bondi_accretion_rate(1.0, sod.adiabatic_index, sod.p.0, sod.rho.0);
    println!("{}", m_dot);
}