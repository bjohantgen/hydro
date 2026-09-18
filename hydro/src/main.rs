// This is the main file of the hydro simulation.
//
// Author: Brayden JoHantgen
// Last Update: 9/17/2026

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

    // Initial Conditions
    let init_setup = sim_func::InitConds{
                        adiabatic_index: 1.4,
                        discontinuity: 0.0,
                        bh_mass: 1.0,
                        softness: 3.0,
                        p: (1.0, 0.0),
                        rho: (1.0, 0.0),
                        vx1: (0.6*(1.4_f64).sqrt(), 0.0), 
                        vx2: (0.0, 0.0),
                        vx3: (0.0, 0.0),
                        bx1: (0.0, 0.0),
                        bx2: (0.0, 0.0),
                        bx3: (0.0, 0.0),
    };

    let drive = sim_func::Driver{
                cfl: 0.3,
                tfinal: 1.001,
                checkpoint: 0.0125,
                x1_cell: 256,
                x1_range: (0.0, 1.0),
                x2_cell: 0,
                x2_range: (0.0, 0.0),
                x3_cell: 0,
                x3_range: (0.0, 0.0),
    };

    // Simulation Prerequisites
    let before = Instant::now();

    let mut t: f64 = 0.0; 
    let mut t_checkpoint = drive.checkpoint;
    let mut time_step_count: f64 = 0.0;
    let mut check_count: i64 = 0;

    let dr = (drive.x1_range.1 - drive.x1_range.0) / (drive.x1_cell as f64);
    let mut prims = sim_func::init_prims_1d(&init_setup, &drive);
    let mut cons: sim_func::Conserved;
    let mut dt: f64;

    println!("{:?}", prims.rho);
    
    let _ = sim_func::write_checkpoint(&prims, &drive, t, check_count);
    check_count += 1;

    // Running the Simulation 
    while t < drive.tfinal {
        cons = sim_func::cons_from_prim_1d(&prims, init_setup.adiabatic_index, dr, init_setup.bh_mass, init_setup.softness);

        dt = sim_func::determine_time_step_1d(&prims, &drive, init_setup.adiabatic_index, dr);

        cons = sim_func::rk_step_1d(&prims, &cons, init_setup.adiabatic_index, dr, dt, init_setup.bh_mass, init_setup.softness);
        prims = sim_func::prim_from_cons_1d(&cons, init_setup.adiabatic_index, dr, init_setup.bh_mass, init_setup.softness);

        if t >= t_checkpoint {
            let _ = sim_func::write_checkpoint(&prims, &drive, t, check_count);
            t_checkpoint += drive.checkpoint;
            check_count += 1;
        }
        t += dt;
        time_step_count += 1.0;
        
        println!("step={} t={} dt={}", time_step_count, t, dt);
    }

    // Simulation Performance
    let runtime = (before.elapsed().as_millis() as f64) / 1000.0;
    let performance = (time_step_count * (prims.p.len() as f64)) / runtime;

    println!("The number of timesteps: {:?}", time_step_count);
    println!("The runtime in seconds: {:?}", runtime);
    println!("The performance in zones per second: {:.2?}", performance);
    
}