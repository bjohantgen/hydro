// This is the main file of the hydro simulation.
//
// Author: Brayden JoHantgen
// Last Update: 11/8/2026

////////////
// Imports
////////////
#![allow(unused_imports)]
use std::time::Instant;
use std::f64::consts::PI;

pub mod math_func;
pub mod sim_func;
pub mod test_conds;

/*
Look for Computation LRD papers, make a list, and add to Lit Review
Radiative Transfer
*/

///////////////
// Simulation
///////////////
fn main() {

    // Initial Conditions
    let init_setup = test_conds::bondi_accretion_state_5();

    let drive = sim_func::Driver{
                cfl: 0.3,
                tfinal: 5.381,
                checkpoint: 0.0125,
                x1_cell: 16,//512,
                x1_range: (0.0, 1.0),
                x2_cell: 0,
                x2_range: (0.0, 0.0),
                x3_cell: 0,
                x3_range: (0.0, 0.0),
    };

    // Print Error Checking
    //let test_prims = sim_func::init_prims(&init_setup, &drive);
    //let test_value = sim_func::cons_from_prims(&test_prims, &init_setup, &drive);
    //println!("{:?}", test_value);

    // Simulation Prerequisites
    let before = Instant::now();

    let mut t: f64 = 0.0; 
    let mut t_checkpoint = drive.checkpoint;
    let mut time_step_count: f64 = 0.0;
    let mut check_count: i64 = 0;

    let mut prims = sim_func::init_prims(&init_setup, &drive);
    sim_func::boundary_conds(&mut prims, &init_setup);
    let mut cons: sim_func::Conserved;
    let mut dt: f64;
    
    let _ = sim_func::write_checkpoint(&prims, &init_setup, &drive, t, check_count);
    check_count += 1;

    // Running the Simulation 
    while t < drive.tfinal {
        cons = sim_func::cons_from_prims(&prims, &init_setup, &drive);

        dt = sim_func::determine_t_step(&prims, &init_setup, &drive);

        cons = sim_func::rk_step(&prims, &cons, &init_setup, &drive, dt);
        prims = sim_func::prims_from_cons(&cons, &init_setup, &drive);
        sim_func::boundary_conds(&mut prims, &init_setup);

        t += dt;
        time_step_count += 1.0;

        if t >= t_checkpoint {
            let _ = sim_func::write_checkpoint(&prims, &init_setup, &drive, t, check_count);
            t_checkpoint += drive.checkpoint;
            check_count += 1;
        }
        
        println!("step={} t={} dt={}", time_step_count, t, dt);
    }

    // Simulation Performance
    let runtime = (before.elapsed().as_millis() as f64) / 1000.0;
    let performance = (time_step_count * (prims.p.len() as f64)) / runtime;

    println!("The number of timesteps: {:?}", time_step_count);
    println!("The runtime in seconds: {:?}", runtime);
    println!("The performance in zones per second: {:.2?}", performance);
}