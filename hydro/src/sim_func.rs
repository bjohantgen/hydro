// This file is full of functions to supplement the "hydro" code.
//
// Author: Brayden JoHantgen
// Last Update: 11/9/2026

///////////////////////////////////
// Importing and Defining Structs
///////////////////////////////////
#![allow(unused_imports)]
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::f64::consts::PI;

use crate::math_func;

pub struct InitConds {
    pub adiabatic_index: f64,
    pub discontinuity: f64,
    pub bh_mass: f64,
    pub softness: f64,
    pub length: f64,
    pub p: (f64, f64),
    pub rho: (f64, f64),
    pub v1: (f64, f64),
    pub v2: (f64, f64),
    pub v3: (f64, f64),
    pub b1: (f64, f64),
    pub b2: (f64, f64),
    pub b3: (f64, f64),
}

pub struct Driver {
    pub cfl: f64,
    pub tfinal: f64,
    pub checkpoint: f64,
    pub x1_cell: usize,
    pub x1_range: (f64, f64),
    pub x2_cell: usize,
    pub x2_range: (f64, f64),
    pub x3_cell: usize,
    pub x3_range: (f64, f64),
}

#[derive(Debug)]
pub struct Primitives {
    pub p: Vec<f64>, 
    pub rho: Vec<f64>,  
    pub v1: Vec<f64>,   
    pub v2: Vec<f64>,   
    pub v3: Vec<f64>,  
}

#[derive(Debug)]
pub struct Conserved {
    pub rho: Vec<f64>,
    pub m1: Vec<f64>,
    pub m2: Vec<f64>,
    pub m3: Vec<f64>,
    pub e: Vec<f64>,
}

pub type Cell8 = (f64, f64, f64, f64, f64, f64, f64, f64);
pub type Cell5 = (f64, f64, f64, f64, f64);
pub type Cell3 = (f64, f64, f64);

//////////////////////////////
// One Dimensional Functions
//////////////////////////////

/// Input:
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     prims (Primitives): the primitive variables for each cell in the grid
/// Description:
///     This funciton uses drive to initialize the grid and then uses init_conds
///     to populate the grid with the initial primitive varibles. The grid is 
///     uniformly spaced.
pub fn init_prims(init_conds: &InitConds, drive: &Driver) -> Primitives {
    let n_cells = drive.x1_cell + 4;

    let mut init_p = vec![0.0; n_cells];
    let mut init_rho = vec![0.0; n_cells];
    let mut init_v1 = vec![0.0; n_cells];
    let mut init_v2 = vec![0.0; n_cells];
    let mut init_v3 = vec![0.0; n_cells];

    for i in 0..n_cells {
        init_p[i] = init_conds.p.0;
        init_rho[i] = init_conds.rho.0;
        init_v1[i] = init_conds.v1.0;
        init_v2[i] = init_conds.v2.0;
        init_v3[i] = init_conds.v3.0;
    }

    let prims = Primitives{
        p: init_p,
        rho: init_rho,
        v1: init_v1,
        v2: init_v2,
        v3: init_v3,
    };

    prims
}

/// Input:
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     delta (f64): the width of the cell in units of R_b
/// Description:
///     This function calculates the width of a cell using the length
///     and cell number given by the Driver construct.
pub fn cell_width(drive: &Driver) -> f64 {
    let delta = (drive.x1_range.1 - drive.x1_range.0) / drive.x1_cell as f64;
    delta
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     i (usize): the index of the desired cell
/// Output:
///     cell (Cell5): the primitive variables of the desired cell
/// Description:
///     This function takes the index of a cell to return the primitive variables of that cell on the grid
pub fn cell_loc(prims: &Primitives, i: usize) -> Cell5 {
    let cell = (prims.p[i], prims.rho[i], prims.v1[i], prims.v2[i], prims.v3[i]);
    cell
}

/// Input:
///     prims (&mut Primitives): the grid containing all the primitive variables for each cell
///     i_org (usize): the index of the original cell
///     i_new (usize): the index of the changed cell
/// Output:
///     the updated prims 
/// Description:
///     This function takes the primitive variables on the grid and replaces a
///     desired index with a different chosen index.
pub fn copy_cell(prims: &mut Primitives, i_org: usize, i_new: usize) {
    prims.p[i_new] = prims.p[i_org];
    prims.rho[i_new] = prims.rho[i_org];
    prims.v1[i_new] = prims.v1[i_org];
    prims.v2[i_new] = prims.v2[i_org];
    prims.v3[i_new] = prims.v3[i_org];
}

/// Input:
///     prims (&mut Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
/// Output:
///     the updated prims with adjusted boundary conditions
/// Description:
///     This function takes the prims grid and adds the boundary conditiond to the ghost cells on either
///     side of the grid. The left bounday is that of a black hole and requires the gas to fall through.
///     The right boundary pulls from an infinite gas reservoir with the same parameters as the initial conditions.
pub fn boundary_conds(prims: &mut Primitives, init_conds: &InitConds) {
    let n_cells = prims.p.len();
    let inner_active = 2;

    // left boundary
    copy_cell(prims, inner_active, 1);
    copy_cell(prims, inner_active, 0);

    // right boundary
    for &i in &[n_cells - 2, n_cells - 1] {
        prims.p[i] = init_conds.p.0;
        prims.rho[i] = init_conds.rho.0;
        prims.v1[i] = init_conds.v1.0;
        prims.v2[i] = init_conds.v2.0;
        prims.v3[i] = init_conds.v3.0;
    }
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
///     m_bh (f64): the mass of the black hole
///     i_num (usize): the index of the cell
///     s (f64): the softening parameter
///     delta (f64): the width of a cell
/// Output:
///     con (Cell5): the conserved variables of the cell
/// Description:
///     This function uses the total energy function to construct the conserved variables
///     using the primitive variables.
pub fn prim_to_con(prim: Cell5, a_ind: f64, m_bh: f64, i_num: usize, s: f64, delta: f64) -> Cell5 {
    let energy = math_func::total_energy(prim, a_ind, m_bh, i_num, s, delta);

    let u0 = prim.1;
    let u1 = prim.1 * prim.2;
    let u2 = prim.1 * prim.3;
    let u3 = prim.1 * prim.4;
    let u4 = energy;

    let con = (u0, u1, u2, u3, u4);
    con
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     cons (Conserved): the conserved variables on the grid.
/// Description:
///     This function uses the cell_width and prim_to_con functions to construct the conserved variables
///     on the simulation grid.
pub fn cons_from_prims(prims: &Primitives, init_conds: &InitConds, drive: &Driver) -> Conserved {
    let n = prims.p.len() as usize;

    let delta = cell_width(drive);

    let mut cons_rho = vec![0.0; n];
    let mut cons_m1 = vec![0.0; n];
    let mut cons_m2 = vec![0.0; n];
    let mut cons_m3 = vec![0.0; n];
    let mut cons_e = vec![0.0; n];

    for i in 0..n {
        let prim_fill = (prims.p[i], prims.rho[i], prims.v1[i], prims.v2[i], prims.v3[i]);
        let con_result = prim_to_con(prim_fill, init_conds.adiabatic_index, init_conds.bh_mass, i, init_conds.softness, delta);

        cons_rho[i] = con_result.0;
        cons_m1[i] = con_result.1;
        cons_m2[i] = con_result.2;
        cons_m3[i] = con_result.3;
        cons_e[i] = con_result.4;
    }

    let cons = Conserved{
        rho: cons_rho,
        m1: cons_m1,
        m2: cons_m2,
        m3: cons_m3,
        e: cons_e,
    };

    cons
}

/// Input:
///     con (Cell5): the conserved variables of a cell
///     a_ind (f64): the adiabatic index of the gas
///     m_bh (f64): the mass of the black hole
///     i_num (usize): the index of the cell
///     s (f64): the softening parameter
///     delta (f64): the width of a cell
/// Output:
///     prim (Cell5): the primitive variables of a cell
/// Description:
///     This function uses the distance and grav_potential functions to transform the conserved variables
///     to the primitive variables. There is a pressure and density floor implemented to prevent these 
///     values from being negative.
pub fn con_to_prim(con: Cell5, a_ind: f64, m_bh: f64, i_num: usize, s: f64, delta: f64) -> Cell5 {
    let rho = con.0.max(1e-8);

    let v1 = con.1 / rho;
    let v2 = con.2 / rho;
    let v3 = con.3 / rho;

    let r = math_func::distance(i_num, delta);
    let phi = math_func::grav_potential(m_bh, r, s, delta);

    let v_square = v1 * v1 + v2 * v2 + v3 * v3;

    let p_check = (a_ind - 1.0) * (con.4 - rho * phi - 0.5 * rho * v_square);
    let p = p_check.max(1e-8);

    let prim = (p, rho, v1, v2, v3);
    prim
}

/// Input:
///     cons (&Conserved): the conserved variables on the grid.
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     prims (Primitives): the grid containing all the primitive variables for each cell
/// Description:
///     This function uses the cell width function and the con_to_prim function to transform the
///     conserved grid to the primitive grid.
pub fn prims_from_cons(cons: &Conserved, init_conds: &InitConds, drive: &Driver) -> Primitives {
    let n = cons.rho.len() as usize;

    let delta = cell_width(drive);

    let mut prims_p = vec![0.0; n];
    let mut prims_rho = vec![0.0; n];
    let mut prims_v1 = vec![0.0; n];
    let mut prims_v2 = vec![0.0; n];
    let mut prims_v3 = vec![0.0; n];

    for i in 0..n {
        let con_fill = (cons.rho[i], cons.m1[i], cons.m2[i], cons.m3[i], cons.e[i]);
        let prim_result = con_to_prim(con_fill, init_conds.adiabatic_index, init_conds.bh_mass, i, init_conds.softness, delta);

        prims_p[i] = prim_result.0;
        prims_rho[i] = prim_result.1;
        prims_v1[i] = prim_result.2;
        prims_v2[i] = prim_result.3;
        prims_v3[i] = prim_result.4;
    }

    let prims = Primitives{
        p: prims_p,
        rho: prims_rho,
        v1: prims_v1,
        v2: prims_v2,
        v3: prims_v3,
    };
    
    prims
}

/// Input:
///     prim_1 (Cell5): the primitive variables in position 1
///     prim_2 (Cell5): the primitive variables in position 2
///     prim_3 (Cell5): the primitive variables in position 3
///     prim_4 (Cell5): the primitive variables in position 4
///     a_ind (f64): the adiabatic index of the gas
///     m_bh (f64): the mass of the black hole
///     i_num1 (usize): the index of cell 2
///     i_num2 (usize): the index of cell 3
///     s (f64): the softening parameter
///     delta (f64): the width of a cell
/// Output:
///     hll (Cell5): the hll plm reconstructed hll flux between cells 2 and 3
/// Description:
///     This function uses the plm_reconstruction, sound_speed, prim_to_con, flux, and tuple max functions
///     to calculate the hll flux.
pub fn hll(prim_1: Cell5, prim_2: Cell5, prim_3: Cell5, prim_4: Cell5, a_ind: f64, m_bh: f64, i_num1: usize, i_num2: usize, s: f64, delta: f64) -> Cell5 {
    let p_l = math_func::plm_reconstruction(prim_2.0, prim_3.0, prim_4.0, prim_1.0, true);
    let rho_l = math_func::plm_reconstruction(prim_2.1, prim_3.1, prim_4.1, prim_1.1, true);
    let v1_l = math_func::plm_reconstruction(prim_2.2, prim_3.2, prim_4.2, prim_1.2, true);
    let v2_l = math_func::plm_reconstruction(prim_2.3, prim_3.3, prim_4.3, prim_1.3, true);
    let v3_l = math_func::plm_reconstruction(prim_2.4, prim_3.4, prim_4.4, prim_1.4, true);

    let p_r = math_func::plm_reconstruction(prim_2.0, prim_3.0, prim_4.0, prim_1.0, false);
    let rho_r = math_func::plm_reconstruction(prim_2.1, prim_3.1, prim_4.1, prim_1.1, false);
    let v1_r = math_func::plm_reconstruction(prim_2.2, prim_3.2, prim_4.2, prim_1.2, false);
    let v2_r = math_func::plm_reconstruction(prim_2.3, prim_3.3, prim_4.3, prim_1.3, false);
    let v3_r = math_func::plm_reconstruction(prim_2.4, prim_3.4, prim_4.4, prim_1.4, false);

    let prim_l = (p_l, rho_l, v1_l, v2_l, v3_l);
    let prim_r = (p_r, rho_r, v1_r, v2_r, v3_r);

    let cf_l = math_func::sound_speed(prim_l, a_ind);
    let plus_l = v1_l + cf_l;
    let minus_l = v1_l - cf_l;
    let u_l = prim_to_con(prim_l, a_ind, m_bh, i_num1, s, delta);
    let f_l = math_func::flux(prim_l, a_ind, m_bh, i_num1, s, delta);

    let cf_r = math_func::sound_speed(prim_r, a_ind);
    let plus_r = v1_r + cf_r;
    let minus_r = v1_r - cf_r;
    let u_r = prim_to_con(prim_r, a_ind, m_bh, i_num2, s, delta);
    let f_r = math_func::flux(prim_r, a_ind, m_bh, i_num2, s, delta);

    let a_plus = math_func::tuple_max((0.0, plus_l, plus_r));
    let a_minus = math_func::tuple_max((0.0, -minus_l, -minus_r));

    let hll_0 = ((a_plus * f_l.0) + (a_minus * f_r.0) - (a_plus * a_minus * (u_r.0 - u_l.0))) / (a_minus + a_plus);
    let hll_1 = ((a_plus * f_l.1) + (a_minus * f_r.1) - (a_plus * a_minus * (u_r.1 - u_l.1))) / (a_minus + a_plus);
    let hll_2 = ((a_plus * f_l.2) + (a_minus * f_r.2) - (a_plus * a_minus * (u_r.2 - u_l.2))) / (a_minus + a_plus);
    let hll_3 = ((a_plus * f_l.3) + (a_minus * f_r.3) - (a_plus * a_minus * (u_r.3 - u_l.3))) / (a_minus + a_plus);
    let hll_4 = ((a_plus * f_l.4) + (a_minus * f_r.4) - (a_plus * a_minus * (u_r.4 - u_l.4))) / (a_minus + a_plus);

    let hll = (hll_0, hll_1, hll_2, hll_3, hll_4);
    hll
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     go (Vec<Cell5>): the vector of godonov fluxes
/// Description:
///     This function calculates all of the hll fluxes on the grid.
pub fn godonov(prims: &Primitives, init_conds: &InitConds, drive: &Driver) -> Vec<Cell5> {
    let n = prims.p.len() as usize;
    let mut go = vec![(0.0, 0.0, 0.0, 0.0, 0.0); n-1];

    let delta = cell_width(drive);

    for i in 1..(n-2) {
        go[i] = hll(cell_loc(prims, i-1), cell_loc(prims, i), cell_loc(prims, i+1), cell_loc(prims, i+2), init_conds.adiabatic_index, init_conds.bh_mass, i, i+1, init_conds.softness, delta);
    }
    go
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     l_vec (Conserved): this is the l function used to update the conserved variables.
/// Description:
///     This function uses the cell_width, godonov, and source functions to calculate the l_functions
///     to update the conserved variables.
pub fn l_function(prims: &Primitives, init_conds: &InitConds, drive: &Driver) -> Conserved {
    let n = prims.p.len() as usize;

    let delta = cell_width(drive);

    let mut l_rho = vec![0.0; n];
    let mut l_m1 = vec![0.0; n];
    let mut l_m2 = vec![0.0; n];
    let mut l_m3 = vec![0.0; n];
    let mut l_e = vec![0.0; n];

    let f_vec = godonov(prims, init_conds, drive);

    for &i in &[0usize, 1, n - 2, n - 1] {
        l_rho[i] = 0.0;
        l_m1[i] = 0.0;
        l_m2[i] = 0.0;
        l_m3[i] = 0.0;
        l_e[i] = 0.0;
    }

    let x_min = drive.x1_range.0;
    for i in 2..(n-2) {
        // active cell j = i-2 spans [x_min + j*delta, x_min + (j+1)*delta]
        let r_l = x_min + (i as f64 - 2.0) * delta;
        let r_r = r_l + delta;
        let r_c = r_l + 0.5 * delta;
        let a_l = r_l * r_l;
        let a_r = r_r * r_r;
        let vol = (r_r.powi(3) - r_l.powi(3)) / 3.0;

        let p_fill = (prims.p[i], prims.rho[i], prims.v1[i], prims.v2[i], prims.v3[i]);
        let s = math_func::source(p_fill, init_conds.bh_mass, r_c, init_conds.softness, delta);

        l_rho[i] = s.0 - (a_r * f_vec[i].0 - a_l * f_vec[i-1].0) / vol;
        l_m1[i]  = s.1 - (a_r * f_vec[i].1 - a_l * f_vec[i-1].1) / vol;
        l_m2[i]  = s.2 - (a_r * f_vec[i].2 - a_l * f_vec[i-1].2) / vol;
        l_m3[i]  = s.3 - (a_r * f_vec[i].3 - a_l * f_vec[i-1].3) / vol;
        l_e[i]   = s.4 - (a_r * f_vec[i].4 - a_l * f_vec[i-1].4) / vol;
    }

    let l_vec = Conserved{
        rho: l_rho,
        m1: l_m1,
        m2: l_m2,
        m3: l_m3,
        e: l_e,
    };

    l_vec
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
///     dt (f64): the size of the timestep
/// Output:
///     cons2 (Conserved): the time stepped conserved variable.
/// Description:
///     This function uses the l_function, prims_from_cons, and boundary_conds functions to advance the
///     conserved variables through and RK3 timestep.
pub fn rk_step(prims: &Primitives, cons: &Conserved, init_conds: &InitConds, drive: &Driver, dt: f64) -> Conserved {
    let n = prims.p.len() as usize;

    let l0 = l_function(prims, init_conds, drive);
    
    let mut con0_rho = vec![0.0; n];
    let mut con0_m1 = vec![0.0; n];
    let mut con0_m2 = vec![0.0; n];
    let mut con0_m3 = vec![0.0; n];
    let mut con0_e = vec![0.0; n];

    for i in 0..n {
        con0_rho[i] = cons.rho[i] + dt * l0.rho[i];
        con0_m1[i] = cons.m1[i] + dt * l0.m1[i];
        con0_m2[i] = cons.m2[i] + dt * l0.m2[i];
        con0_m3[i] = cons.m3[i] + dt * l0.m3[i];
        con0_e[i] = cons.e[i] + dt * l0.e[i];
    }

    let cons0 = Conserved{
        rho: con0_rho,
        m1: con0_m1, 
        m2: con0_m2,
        m3: con0_m3,
        e: con0_e,
    };

    let mut prims_1 = prims_from_cons(&cons0, init_conds, drive);
    boundary_conds(&mut prims_1, init_conds);

    let l1 = l_function(&prims_1, init_conds, drive);

    let mut con1_rho = vec![0.0; n];
    let mut con1_m1 = vec![0.0; n];
    let mut con1_m2 = vec![0.0; n];
    let mut con1_m3 = vec![0.0; n];
    let mut con1_e = vec![0.0; n];

    let dt25 = 0.25 * dt;

    for i in 0..n{
        con1_rho[i] = 0.75 * cons.rho[i] + 0.25 * cons0.rho[i] + dt25 * l1.rho[i];
        con1_m1[i] = 0.75 * cons.m1[i] + 0.25 * cons0.m1[i] + dt25 * l1.m1[i];
        con1_m2[i] = 0.75 * cons.m2[i] + 0.25 * cons0.m2[i] + dt25 * l1.m2[i];
        con1_m3[i] = 0.75 * cons.m3[i] + 0.25 * cons0.m3[i] + dt25 * l1.m3[i];
        con1_e[i] = 0.75 * cons.e[i] + 0.25 * cons0.e[i] + dt25 * l1.e[i];
    }

    let cons1 = Conserved{
        rho: con1_rho,
        m1: con1_m1, 
        m2: con1_m2,
        m3: con1_m3,
        e: con1_e,
    };

    let mut prims_2 = prims_from_cons(&cons1, init_conds, drive);
    boundary_conds(&mut prims_2, init_conds);

    let l2 = l_function(&prims_2, init_conds, drive);

    let mut con2_rho = vec![0.0; n];
    let mut con2_m1 = vec![0.0; n];
    let mut con2_m2 = vec![0.0; n];
    let mut con2_m3 = vec![0.0; n];
    let mut con2_e = vec![0.0; n];

    let dt66 = (2.0 / 3.0) * dt;
    let c33 = 1.0 / 3.0;
    let c66 = 2.0 / 3.0;

    for i in 0..n{
        con2_rho[i] = c33 * cons.rho[i] + c66 * cons1.rho[i] + dt66 * l2.rho[i];
        con2_m1[i] = c33 * cons.m1[i] + c66 * cons1.m1[i] + dt66 * l2.m1[i];
        con2_m2[i] = c33 * cons.m2[i] + c66 * cons1.m2[i] + dt66 * l2.m2[i];
        con2_m3[i] = c33 * cons.m3[i] + c66 * cons1.m3[i] + dt66 * l2.m3[i];
        con2_e[i] = c33 * cons.e[i] + c66 * cons1.e[i] + dt66 * l2.e[i];
    }

    let cons2 = Conserved{
        rho: con2_rho,
        m1: con2_m1, 
        m2: con2_m2,
        m3: con2_m3,
        e: con2_e,
    };

    cons2
}

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     init_conds (&InitConds): the initial conditions of the simulation
///     drive (&Driver): the simulation parameters
/// Output:
///     dt (f64): the time step
/// Description:
///     This function uses the cell_width, compute_time_step, distance, and gravity_time_step functions to 
///     calculate the time step satisfying the Courant-Friedrich-Levy condition. This function checks for 
///     non-finite values from non physical primitive variables.
pub fn determine_t_step(prims: &Primitives, init_conds: &InitConds, drive: &Driver) -> f64 {
    let mut dt = f64::INFINITY;
    let n = prims.p.len();

    let delta = cell_width(drive);

        for i in 1..(n - 2) {
            let prim_fill_1 = (prims.p[i], prims.rho[i], prims.v1[i], prims.v2[i], prims.v3[i]);
            let prim_fill_2 = (prims.p[i+1], prims.rho[i+1], prims.v1[i+1], prims.v2[i+1], prims.v3[i+1]);

            let dt_check = math_func::compute_time_step(prim_fill_1, prim_fill_2, init_conds.adiabatic_index, delta);

            assert!(
                dt_check.is_finite(),
                "non-finite dt at interface {}: p=({}, {}) rho=({}, {})",
                i, prims.p[i], prims.p[i + 1], prims.rho[i], prims.rho[i + 1]
            );

            dt = dt.min(dt_check);
        }

        for i in 2..(n - 2) {
            let r = math_func::distance(i, delta);
            let dt_grav = math_func::gravity_time_step(init_conds.bh_mass, r, init_conds.softness, delta);
            dt = dt.min(dt_grav);
        }

    dt = drive.cfl * dt;
    dt
}

//////////////////////
// General Functions
//////////////////////

/// Input:
///     prims (&Primitives): the grid containing all the primitive variables for each cell
///     drive (&Driver): the simulation parameters
///     t (f64): the time of the checkpoint
///     check_count (i64): the chcekpoint number
/// Output:
///     A file containing the recorded values
/// Description:
///     This function writes the number of cells, grid range, and all primitive variables into txt files for plotting.
pub fn write_checkpoint(prims: &Primitives, drive: &Driver, t: f64, check_count: i64) -> Result<(), Box<dyn std::error::Error>> {
    let file_num = check_count.to_string();
    let file_type = ".txt".to_string();
    let file_name = format!("{}{}", file_num, file_type);
    let file_path = "time_step_files/".to_string() + &file_name;
    let output_file_path = Path::new(&file_path);
    
    let t_fill = format!("{} {}", "t:".to_string(), &(t.to_string()+&" ".to_string()));
    
    let size = prims.p.len();

    let num_fill_x1 = format!("{} {}", "cell_num_x1:".to_string(), &(drive.x1_cell.to_string()+&" ".to_string()));
    let left_x1_fill = format!("{} {}", "left_x1_bound:".to_string(), &(drive.x1_range.0.to_string()+&" ".to_string()));
    let right_x1_fill = format!("{} {}", "right_x1_bound:".to_string(), &(drive.x1_range.1.to_string()+&" ".to_string()));
    
    let mut p_string = "p: ".to_string();
    let mut rho_string = "rho: ".to_string();
    let mut v1_string = "v1: ".to_string();
    let mut v2_string = "v2: ".to_string();
    let mut v3_string = "v3: ".to_string();

    //let mut mdot_string = "M_dot: ".to_string();
        
    for i in 0..size {
        let p_fill = prims.p[i].to_string() + &" ".to_string();
        p_string.push_str(&p_fill);
    
        let rho_fill = prims.rho[i].to_string() + &" ".to_string();
        rho_string.push_str(&rho_fill);
    
        let v1_fill = prims.v1[i].to_string() + &" ".to_string();
        v1_string.push_str(&v1_fill);
    
        let v2_fill = prims.v2[i].to_string() + &" ".to_string();
        v2_string.push_str(&v2_fill);
    
        let v3_fill = prims.v3[i].to_string() + &" ".to_string();
        v3_string.push_str(&v3_fill);
    }

    // let dr = (drive.x1_range.1 - drive.x1_range.0) / (drive.x1_cell as f64);
    // let r = (math_func::radial_distance((drive.x1_cell as f64 /2.0) as usize, dr) + math_func::radial_distance(((drive.x1_cell as f64 /2.0)+1.0) as usize, dr)) / 2.0;
    //let mdot_fill = math_func::accretion_rate(r, prims.rho[(drive.x1_cell as f64 /2.0) as usize], prims.vx1[(drive.x1_cell as f64 /2.0) as usize]).to_string() 
    //let index = (drive.x1_cell as f64 /2.0) as usize;
    //let mdot_fill = (4.0 * PI * r * r * hll_1d(cell_loc(prims, index-1), cell_loc(prims, index), cell_loc(prims, index+1), cell_loc(prims, index+2), init_setup.adiabatic_index, index, index+1, dr, init_setup.bh_mass, init_setup.softness).0).to_string() + &" ".to_string();
    //mdot_string.push_str(&mdot_fill);
     
    let file = fs::OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output_file_path)?;
    
    let mut writer = BufWriter::new(file);
    
    writeln!(writer, "{}", num_fill_x1)?;
    writeln!(writer, "{}", left_x1_fill)?;
    writeln!(writer, "{}", right_x1_fill)?;

    writeln!(writer, "{}", t_fill)?;

    writeln!(writer, "{}", p_string)?;
    writeln!(writer, "{}", rho_string)?;
    writeln!(writer, "{}", v1_string)?;
    writeln!(writer, "{}", v2_string)?;
    writeln!(writer, "{}", v3_string)?;

    //writeln!(writer, "{}", mdot_string)?;

    writer.flush()?;
    
    Ok(())
}