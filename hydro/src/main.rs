// This is the main file of the hydro simulation.
//
// Author: Brayden JoHantgen
// Last Update: 9/13/2026

////////////
// Imports
////////////
#![allow(unused_imports)]
use std::time::Instant;
use std::f64::consts::PI;

pub mod math_func;
pub mod sim_func;

///////////////
// Simulation
///////////////
fn main() {

// Initial Conditions
    let icons = sim_func::InitConds{
                        adiabatic_index: 1.4,
                        discontinuity: 0.5,
                        p: (1.0, 0.125),
                        rho: (1.0, 0.1),
                        vx1: (0.0, 0.0), 
                        vx2: (0.0, 0.0),
                        vx3: (0.0, 0.0),
                        bx1: (0.0, 0.0),
                        bx2: (0.0, 0.0),
                        bx3: (0.0, 0.0),

    };

    let m_dot = math_func::bondi_accretion_rate(1.0, icons.adiabatic_index, icons.p.0, icons.rho.0);
    println!("{}", m_dot);

}