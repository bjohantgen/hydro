// This file is full of the initial conditions of famous test problems
// to supplement the "hydro" code.
//
// Author: Brayden JoHantgen
// Last Update: 9/16/2026

//////////////
// Importing 
//////////////
use crate::sim_func;

/////////////////////////////////////
// Defining Functions of Test Cases
/////////////////////////////////////

// Sod Shock Tube
pub fn sod_shock() -> sim_func::InitConds {
    let sod = sim_func::InitConds{
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
    sod
}

// Harder Shock Tube
pub fn hard_shock() -> sim_func::InitConds {
    let hard = sim_func::InitConds{
        adiabatic_index: 1.4,
        discontinuity: 0.5,
        p: (100.0, 1.0),
        rho: (10.0, 1.0),
        vx1: (0.0, 0.0), 
        vx2: (0.0, 0.0),
        vx3: (0.0, 0.0),
        bx1: (0.0, 0.0),
        bx2: (0.0, 0.0),
        bx3: (0.0, 0.0),
    };
    hard
}

// Brio-Wu Shock Tube