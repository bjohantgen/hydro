// This file is full of the initial conditions of famous test problems
// to supplement the "hydro" code.
//
// Author: Brayden JoHantgen
// Last Update: 11/9/2026

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
        bh_mass: 0.0,
        softness: 0.0,
        length: 0.0,
        p: (1.0, 0.125),
        rho: (1.0, 0.1),
        v1: (0.0, 0.0), 
        v2: (0.0, 0.0),
        v3: (0.0, 0.0),
        b1: (0.0, 0.0),
        b2: (0.0, 0.0),
        b3: (0.0, 0.0),
    };
    sod
}

// Harder Shock Tube
pub fn hard_shock() -> sim_func::InitConds {
    let hard = sim_func::InitConds{
        adiabatic_index: 1.4,
        discontinuity: 0.5,
        bh_mass: 0.0,
        softness: 0.0,
        length: 0.0,
        p: (100.0, 1.0),
        rho: (10.0, 1.0),
        v1: (0.0, 0.0), 
        v2: (0.0, 0.0),
        v3: (0.0, 0.0),
        b1: (0.0, 0.0),
        b2: (0.0, 0.0),
        b3: (0.0, 0.0),
    };
    hard
}

// Brio-Wu Shock Tube
pub fn brio_wu_shock() -> sim_func::InitConds {
    let bw = sim_func::InitConds{
        adiabatic_index: 2.0,
        discontinuity: 0.5,
        bh_mass: 0.0,
        softness: 0.0,
        length: 0.0,
        p: (1.0, 0.1),
        rho: (1.0, 0.125),
        v1: (0.0, 0.0), 
        v2: (0.0, 0.0),
        v3: (0.0, 0.0),
        b1: (0.75, 0.75),
        b2: (1.0, -1.0),
        b3: (0.0, 0.0),
    };
    bw
}

// Bondi-Accretion
pub fn bondi_accretion_state_5() -> sim_func::InitConds {
    let bondi = sim_func::InitConds{
        adiabatic_index: 5.0 / 3.0,
        discontinuity: 0.0,
        bh_mass: 0.5,
        softness: 5.0,
        length: 64.0,
        p: (3.0 / 5.0, 0.0),
        rho: (1.0, 0.0),
        v1: (0.0, 0.0), 
        v2: (0.0, 0.0),
        v3: (0.0, 0.0),
        b1: (0.0, 0.0),
        b2: (0.0, 0.0),
        b3: (0.0, 0.0),
    };
    bondi
}