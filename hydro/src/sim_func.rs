// This file is full of functions to supplement the "hydro" code.
//
// Author: Brayden JoHantgen
// Last Update: 9/8/2026

///////////////////////////////////
// Importing and Defining Structs
///////////////////////////////////
#![allow(unused_imports)]
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::f64::consts::PI;

use crate::math_func;

pub struct Init_Conds {
    pub adiabatic_index: f64,
    pub discontinuity: f64,
    pub p: (f64, f64),
    pub rho: (f64, f64),
    pub vx: (f64, f64),
    pub vy: (f64, f64),
    pub vz: (f64, f64),
    pub bx: (f64, f64),
    pub by: (f64, f64),
    pub bz: (f64, f64),
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
    pub x3_range: (f64, f64)
}

#[derive(Debug)]
pub struct Primitives {
    pub p: Vec<f64>, 
    pub rho: Vec<f64>,  
    pub vx1: Vec<f64>,   
    pub vx2: Vec<f64>,   
    pub vx3: Vec<f64>,   
    pub bx1: Vec<f64>,    
    pub bx2: Vec<f64>,     
    pub bx3: Vec<f64>,   
}

#[derive(Debug)]
pub struct Conserved {
    pub rho: Vec<f64>,
    pub mx1: Vec<f64>,
    pub mx2: Vec<f64>,
    pub mx3: Vec<f64>,
    pub bx1: Vec<f64>,
    pub bx2: Vec<f64>,
    pub bx3: Vec<f64>,
    pub e: Vec<f64>,
}

pub type Cell8 = (f64, f64, f64, f64, f64, f64, f64, f64);
pub type Cell5 = (f64, f64, f64, f64, f64);
pub type Cell3 = (f64, f64, f64);

//////////////////////////////
// One Dimensional Functions
//////////////////////////////

/// Input:
/// Output:
/// Description:
pub fn init_prims_1d(init_conds: &Init_Conds, drive: &Driver) -> Primitives {
    let ncells = drive.x1_cell + 2;

    let mut init_p = vec![0.0; ncells];
    let mut init_rho = vec![0.0; ncells];
    let mut init_vx1 = vec![0.0; ncells];
    let mut init_vx2 = vec![0.0; ncells];
    let mut init_vx3 = vec![0.0; ncells];
    let mut init_bx1 = vec![0.0; ncells];
    let mut init_bx2 = vec![0.0; ncells];
    let mut init_bx3 = vec![0.0; ncells];

    for i in 0..ncells {
        if i <=(((drive.x1_cell+1) as f64) * init_conds.discontinuity) as usize {
            init_p[i] = init_conds.p.0;
            init_rho[i] = init_conds.rho.0;
            init_vx1[i] = init_conds.vx1.0;
            init_vx2[i] = init_conds.vx2.0;
            init_vx3[i] = init_conds.vx3.0;
            init_bx1[i] = init_conds.bx1.0;
            init_bx2[i] = init_conds.bx2.0;
            init_bx3[i] = init_conds.bx3.0;
        } else {
            init_p[i] = init_conds.p.1;
            init_rho[i] = init_conds.rho.1;
            init_vx1[i] = init_conds.vx1.1;
            init_vx2[i] = init_conds.vx2.1;
            init_vx3[i] = init_conds.vx3.1;
            init_bx1[i] = init_conds.bx1.1;
            init_bx2[i] = init_conds.bx2.1;
            init_bx3[i] = init_conds.bx3.1;
        }
    }

    let prims = Primitives{
        p: init_p,
        rho: init_rho,
        vx1: init_vx1,
        vx2: init_vx2,
        vx3: init_vx3,
        bx1: init_bx1,
        bx2: init_bx2,
        bx3: init_bx3,
    };

    prims
}

/// Input:
/// Output:
/// Description:
pub fn prim_to_cons_1d(prim: Cell8, a_index: f64) -> Cell8 {
    let e = math_func::total_energy(prim, a_index);
    let con = (prim.1, prim.1 * prim.2, prim.1 * prim.3, prim.1 * prim.4, prim.5, prim.6, prim.7, e);
    con
}

/// Input:
/// Output:
/// Description:
pub fn cons_from_prim_1d(prims: &Primitives, a_index: f64) -> Conserved {
    let size = prims.p.len() as usize;

    let mut cons_rho = vec![0.0; size];
    let mut cons_mx1 = vec![0.0; size];
    let mut cons_mx2 = vec![0.0; size];
    let mut cons_mx3 = vec![0.0; size];
    let mut cons_bx1 = vec![0.0; size];
    let mut cons_bx2 = vec![0.0; size];
    let mut cons_bx3 = vec![0.0; size];
    let mut cons_e = vec![0.0; size];

    for i in 0..size {
        let prim_fill = (prims.p[i], prims.rho[i], prims.vx1[i], prims.vx2[i], prims.vx3[i], prims.bx1[i], prims.bx2[i], prims.bx3[i]);
        let con_result = prim_to_cons_1d(prim_fill, a_index);

        cons_rho[i] = con_result.0;
        cons_mx1[i] = con_result.1;
        cons_mx2[i] = con_result.2;
        cons_mx3[i] = con_result.3;
        cons_bx1[i] = con_result.4;
        cons_bx2[i] = con_result.5;
        cons_bx3[i] = con_result.6;
        cons_e[i] = con_result.7;
    }

    let cons = Conserved{
        rho: cons_rho,
        mx1: cons_mx1,
        mx2: cons_mx2,
        mx3: cons_mx3,
        bx1: cons_bx1,
        bx2: cons_bx2,
        bx3: cons_bx3,
        e: cons_e,
    };

    cons
}

/// Input:
///     con: the eight conserved variables of a single cell in the 1d grid (Cell8).
///     a_index: the adiabatic index of the gas in the simulation.
/// Output:
///     prim: the physical components of a single cell in the 1d grid (Cell8).
/// Description:
///     This function takes the conserved variables of a single cell and reconstructs the primitive variables for that cell.
pub fn cons_to_prim_1d(con: Cell8, a_index: f64) -> Cell8 {
    let vx1 = con.1 / con.0;
    let vx2 = con.2 / con.0;
    let vx3 = con.3 / con.0;
    let p = (a_index - 1.0) * (con.7 - 0.5 * (con.0 * (vx1 * vx1 + vx2 * vx2 + vx3 * vx3) + (con.4 * con.4 + con.5 * con.5 + con.6 * con.6)));
    let prim = (p, con.0, vx1, vx2, vx3, con.4, con.5, con.6);
    prim
}

/// Input:
/// Output:
/// Description:
pub fn prim_from_cons_1d(cons: &Conserved, a_index: f64) -> Primitives {
    let size = cons.rho.len() as usize;

    let mut prims_p = vec![0.0; size];
    let mut prims_rho = vec![0.0; size];
    let mut prims_vx1 = vec![0.0; size];
    let mut prims_vx2 = vec![0.0; size];
    let mut prims_vx3 = vec![0.0; size];
    let mut prims_bx1 = vec![0.0; size];
    let mut prims_bx2 = vec![0.0; size];
    let mut prims_bx3 = vec![0.0; size];

    for i in 0..size {
        let con_fill = (cons.rho[i], cons.mx1[i], cons.mx2[i], cons.mx3[i], cons.bx1[i], cons.bx2[i], cons.bx3[i], cons.e[i]);
        let prim_result = cons_to_prim_1d(con_fill, a_index);

        prims_p[i] = prim_result.0;
        prims_rho[i] = prim_result.1;
        prims_vx1[i] = prim_result.2;
        prims_vx2[i] = prim_result.3;
        prims_vx3[i] = prim_result.4;
        prims_bx1[i] = prim_result.5;
        prims_bx2[i] = prim_result.6;
        prims_bx3[i] = prim_result.7;
    }

    let prims = Primitives{
        p: prims_p,
        rho: prims_rho,
        vx1: prims_vx1,
        vx2: prims_vx2,
        vx3: prims_vx3,
        bx1: prims_bx1,
        bx2: prims_bx2,
        bx3: prims_bx3,
    };
    
    prims
}

/// Input:
/// Output:
/// Description:
pub fn hll_flux_1d(prim_1: Cell8, prim_2: Cell8, prim_3: Cell8, prim_4: Cell8, a_index: f64) -> Cell8 {
    let p_l = math_func::plm_reconstruction(prim_2.0, prim_3.0, prim_4.0, prim_1.0, true);
    let rho_l = math_func::plm_reconstruction(prim_2.1, prim_3.1, prim_4.1, prim_1.1, true);
    let vx1_l = math_func::plm_reconstruction(prim_2.2, prim_3.2, prim_4.2, prim_1.2, true);
    let vx2_l = math_func::plm_reconstruction(prim_2.3, prim_3.3, prim_4.3, prim_1.3, true);
    let vx3_l = math_func::plm_reconstruction(prim_2.4, prim_3.4, prim_4.4, prim_1.4, true);
    let bx1_l = math_func::plm_reconstruction(prim_2.5, prim_3.5, prim_4.5, prim_1.5, true);
    let bx2_l = math_func::plm_reconstruction(prim_2.6, prim_3.6, prim_4.6, prim_1.6, true);
    let bx3_l = math_func::plm_reconstruction(prim_2.7, prim_3.7, prim_4.7, prim_1.7, true);

    let p_r = math_func::plm_reconstruction(prim_2.0, prim_3.0, prim_4.0, prim_1.0, false);
    let rho_r = math_func::plm_reconstruction(prim_2.1, prim_3.1, prim_4.1, prim_1.1, false);
    let vx1_r = math_func::plm_reconstruction(prim_2.2, prim_3.2, prim_4.2, prim_1.2, false);
    let vx2_r = math_func::plm_reconstruction(prim_2.3, prim_3.3, prim_4.3, prim_1.3, false);
    let vx3_r = math_func::plm_reconstruction(prim_2.4, prim_3.4, prim_4.4, prim_1.4, false);
    let bx1_r = math_func::plm_reconstruction(prim_2.5, prim_3.5, prim_4.5, prim_1.5, false);
    let bx2_r = math_func::plm_reconstruction(prim_2.6, prim_3.6, prim_4.6, prim_1.6, false);
    let bx3_r = math_func::plm_reconstruction(prim_2.7, prim_3.7, prim_4.7, prim_1.7, false);

    let prim_l = (p_l, rho_l, vx1_l, vx2_l, vx3_l, bx1_l, bx2_l, bx3_l);
    let prim_r = (p_r, rho_r, vx1_r, vx2_r, vx3_r, bx1_r, bx2_r, bx3_r);

    let cf_l = math_func::fast_magsonic_speed_x(prim_l, a_index);
    let plus_l = vx1_l + cf_l;
    let minus_l = vx1_l - cf_l;
    let u_l = prim_to_cons_1d(prim_l, a_index);
    let f_l = math_func::flux_x(prim_l, a_index);

    let cf_r = math_func::fast_magsonic_speed_x(prim_r, a_index);
    let plus_r = vx1_r + cf_r;
    let minus_r = vx1_r - cf_r;
    let u_r = prim_to_cons_1d(prim_r, a_index);
    let f_r = math_func::flux_x(prim_r, a_index);

    let a_plus = math_func::tuple_max((0.0, plus_l, plus_r));
    let a_minus = math_func::tuple_max((0.0, -minus_l, -minus_r));

    let hll_0 = ((a_plus * f_l.0) + (a_minus * f_r.0) - (a_plus * a_minus * (u_r.0 - u_l.0))) / (a_minus + a_plus);
    let hll_1 = ((a_plus * f_l.1) + (a_minus * f_r.1) - (a_plus * a_minus * (u_r.1 - u_l.1))) / (a_minus + a_plus);
    let hll_2 = ((a_plus * f_l.2) + (a_minus * f_r.2) - (a_plus * a_minus * (u_r.2 - u_l.2))) / (a_minus + a_plus);
    let hll_3 = ((a_plus * f_l.3) + (a_minus * f_r.3) - (a_plus * a_minus * (u_r.3 - u_l.3))) / (a_minus + a_plus);
    let hll_4 = ((a_plus * f_l.4) + (a_minus * f_r.4) - (a_plus * a_minus * (u_r.4 - u_l.4))) / (a_minus + a_plus);
    let hll_5 = ((a_plus * f_l.5) + (a_minus * f_r.5) - (a_plus * a_minus * (u_r.5 - u_l.5))) / (a_minus + a_plus);
    let hll_6 = ((a_plus * f_l.6) + (a_minus * f_r.6) - (a_plus * a_minus * (u_r.6 - u_l.6))) / (a_minus + a_plus);
    let hll_7 = ((a_plus * f_l.7) + (a_minus * f_r.7) - (a_plus * a_minus * (u_r.7 - u_l.7))) / (a_minus + a_plus);
    let hll = (hll_0, hll_1, hll_2, hll_3, hll_4, hll_5, hll_6, hll_7);
    hll
}

/// Input:
/// Output:
/// Description:
pub fn godonov_1d(prims: &Primitives, a_index: f64) ->  Vec<Cell8> {
    let size = prims.p.len() as usize;
    let mut go = vec![(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); size-1];
    for i in 0..(size-1) {
        if i == 0 {
            let prim_0 = (prims.p[0], prims.rho[0], prims.vx1[0], prims.vx2[0], prims.vx3[0], prims.bx1[0], prims.bx2[0], prims.bx3[0]);
            let prim_1 = (prims.p[1], prims.rho[1], prims.vx1[1], prims.vx2[1], prims.vx3[1], prims.bx1[1], prims.bx2[1], prims.bx3[1]);
            let prim_2 = (prims.p[2], prims.rho[2], prims.vx1[2], prims.vx2[2], prims.vx3[2], prims.bx1[2], prims.bx2[2], prims.bx3[2]);
            let prim_3 = (prims.p[3], prims.rho[3], prims.vx1[3], prims.vx2[3], prims.vx3[3], prims.bx1[3], prims.bx2[3], prims.bx3[3]);
            go[i] = hll_flux_1d(prim_0, prim_1, prim_2, prim_3, a_index);
        } else if i < (size - 2) {
            let prim_fst = (prims.p[i-1], prims.rho[i-1], prims.vx1[i-1], prims.vx2[i-1], prims.vx3[i-1], prims.bx1[i-1], prims.bx2[i-1], prims.bx3[i-1]);
            let prim_s = (prims.p[i], prims.rho[i], prims.vx1[i], prims.vx2[i], prims.vx3[i], prims.bx1[i], prims.bx2[i], prims.bx3[i]);
            let prim_t = (prims.p[i+1], prims.rho[i+1], prims.vx1[i+1], prims.vx2[i+1], prims.vx3[i+1], prims.bx1[i+1], prims.bx2[i+1], prims.bx3[i+1]);
            let prim_frt = (prims.p[i+2], prims.rho[i+2], prims.vx1[i+2], prims.vx2[i+2], prims.vx3[i+2], prims.bx1[i+2], prims.bx2[i+2], prims.bx3[i+2]);
            go[i] = hll_flux_1d(prim_fst, prim_s, prim_t, prim_frt, a_index);
        } else {
            let prim_0 = (prims.p[size-4], prims.rho[size-4], prims.vx1[size-4], prims.vx2[size-4], prims.vx3[size-4], prims.bx1[size-4], prims.bx2[size-4], prims.bx3[size-4]);
            let prim_1 = (prims.p[size-3], prims.rho[size-3], prims.vx1[size-3], prims.vx2[size-3], prims.vx3[size-3], prims.bx1[size-3], prims.bx2[size-3], prims.bx3[size-3]);
            let prim_2 = (prims.p[size-2], prims.rho[size-2], prims.vx1[size-2], prims.vx2[size-2], prims.vx3[size-2], prims.bx1[size-2], prims.bx2[size-2], prims.bx3[size-2]);
            let prim_3 = (prims.p[size-1], prims.rho[size-1], prims.vx1[size-1], prims.vx2[size-1], prims.vx3[size-1], prims.bx1[size-1], prims.bx2[size-1], prims.bx3[size-1]);
            go[i] = hll_flux_1d(prim_0, prim_1, prim_2, prim_3, a_index);
        }
    }
    go
}

/// Input:
/// Output:
/// Description:
pub fn l_function_1d(prims: &Primitives, a_index: f64, dx: f64) -> Conserved {
    let size = prims.p.len() as usize;

    let mut l_rho = vec![0.0; size];
    let mut l_mx1 = vec![0.0; size];
    let mut l_mx2 = vec![0.0; size];
    let mut l_mx3 = vec![0.0; size];
    let mut l_bx1 = vec![0.0; size];
    let mut l_bx2 = vec![0.0; size];
    let mut l_bx3 = vec![0.0; size];
    let mut l_e = vec![0.0; size];

    let f_vec = godonov_1d(prims, a_index);

    l_rho[0] = 0.0;
    l_mx1[0] = 0.0;
    l_mx2[0] = 0.0;
    l_mx3[0] = 0.0;
    l_bx1[0] = 0.0;
    l_bx2[0] = 0.0;
    l_bx3[0] = 0.0;
    l_e[0] = 0.0;

    for i in 1..(size-1) {
        l_rho[i] = - (f_vec[i].0 - f_vec[i-1].0) / dx;
        l_mx1[i] = - (f_vec[i].1 - f_vec[i-1].1) / dx;
        l_mx2[i] = - (f_vec[i].2 - f_vec[i-1].2) / dx;
        l_mx3[i] = - (f_vec[i].3 - f_vec[i-1].3) / dx;
        l_bx1[i] = - (f_vec[i].4 - f_vec[i-1].4) / dx;
        l_bx2[i] = - (f_vec[i].5 - f_vec[i-1].5) / dx;
        l_bx3[i] = - (f_vec[i].6 - f_vec[i-1].6) / dx;
        l_e[i] = - (f_vec[i].7 - f_vec[i-1].7) / dx;
    }

    l_rho[size-1] = 0.0;
    l_mx1[size-1] = 0.0;
    l_mx2[size-1] = 0.0;
    l_mx3[size-1] = 0.0;
    l_bx1[size-1] = 0.0;
    l_bx2[size-1] = 0.0;
    l_bx3[size-1] = 0.0;
    l_e[size-1] = 0.0;

    let l_vec = Conserved{
        rho: l_rho,
        mx1: l_mx2,
        mx2: l_mx2,
        mx3: l_mx3,
        bx1: l_bx1,
        bx2: l_bx2,
        bx3: l_bx3,
        e: l_e,
    };

    l_vec
}

/// Input:
/// Output:
/// Description:
pub fn rk_step_1d(prims: &Primitives, cons: &Conserved, a_index: f64, dx: f64, dt: f64) -> Conserved {
    let size = prims.p.len() as usize;
    let l0 = l_function_1d(prims, a_index, dx);
    
    let mut con0_rho = vec![0.0; size];
    let mut con0_mx1 = vec![0.0; size];
    let mut con0_mx2 = vec![0.0; size];
    let mut con0_mx3 = vec![0.0; size];
    let mut con0_bx1 = vec![0.0; size];
    let mut con0_bx2 = vec![0.0; size];
    let mut con0_bx3 = vec![0.0; size];
    let mut con0_e = vec![0.0; size];

    for i in 0..size {
        con0_rho[i] = cons.rho[i] + dt * l0.rho[i];
        con0_mx1[i] = cons.mx1[i] + dt * l0.mx1[i];
        con0_mx2[i] = cons.mx2[i] + dt * l0.mx2[i];
        con0_mx3[i] = cons.mx3[i] + dt * l0.mx3[i];
        con0_bx1[i] = cons.bx1[i] + dt * l0.bx1[i];
        con0_bx2[i] = cons.bx2[i] + dt * l0.bx2[i];
        con0_bx3[i] = cons.bx3[i] + dt * l0.bx3[i];
        con0_e[i] = cons.e[i] + dt * l0.e[i];
    }

    let cons0 = Conserved{
        rho: con0_rho,
        mx1: con0_mx1, 
        mx2: con0_mx2,
        mx3: con0_mx3,
        bx1: con0_bx1,
        bx2: con0_bx2,
        bx3: con0_bx3,
        e: con0_e,
    };

    let prims_1 = prim_from_cons_1d(&cons0, a_index);
    let l1 = l_function_1d(&prims_1, a_index, dx);

    let mut con1_rho = vec![0.0; size];
    let mut con1_mx1 = vec![0.0; size];
    let mut con1_mx2 = vec![0.0; size];
    let mut con1_mx3 = vec![0.0; size];
    let mut con1_bx1 = vec![0.0; size];
    let mut con1_bx2 = vec![0.0; size];
    let mut con1_bx3 = vec![0.0; size];
    let mut con1_e = vec![0.0; size];

    let dt25 = 0.25 * dt;

    for i in 0..size{
        con1_rho[i] = 0.75 * cons.rho[i] + 0.25 * cons0.rho[i] + dt25 * l1.rho[i];
        con1_mx1[i] = 0.75 * cons.mx1[i] + 0.25 * cons0.mx1[i] + dt25 * l1.mx1[i];
        con1_mx2[i] = 0.75 * cons.mx2[i] + 0.25 * cons0.mx2[i] + dt25 * l1.mx2[i];
        con1_mx3[i] = 0.75 * cons.mx3[i] + 0.25 * cons0.mx3[i] + dt25 * l1.mx3[i];
        con1_bx1[i] = 0.75 * cons.bx1[i] + 0.25 * cons0.bx1[i] + dt25 * l1.bx1[i];
        con1_bx2[i] = 0.75 * cons.bx2[i] + 0.25 * cons0.bx2[i] + dt25 * l1.bx2[i];
        con1_bx3[i] = 0.75 * cons.bx3[i] + 0.25 * cons0.bx3[i] + dt25 * l1.bx3[i];
        con1_e[i] = 0.75 * cons.e[i] + 0.25 * cons0.e[i] + dt25 * l1.e[i];
    }

    let cons1 = Conserved{
        rho: con1_rho,
        mx1: con1_mx1, 
        mx2: con1_mx2,
        mx3: con1_mx3,
        bx1: con1_bx1,
        bx2: con1_bx2,
        bx3: con1_bx3,
        e: con1_e,
    };

    let prims_2 = prim_from_cons_1d(&cons1, a_index);
    let l2 = l_function_1d(&prims_2, a_index, dx);

    let mut con2_rho = vec![0.0; size];
    let mut con2_mx1 = vec![0.0; size];
    let mut con2_mx2 = vec![0.0; size];
    let mut con2_mx3 = vec![0.0; size];
    let mut con2_bx1 = vec![0.0; size];
    let mut con2_bx2 = vec![0.0; size];
    let mut con2_bx3 = vec![0.0; size];
    let mut con2_e = vec![0.0; size];

    let dt66 = (2.0 / 3.0) * dt;
    let c33 = 1.0 / 3.0;
    let c66 = 2.0 / 3.0;

    for i in 0..size{
        con2_rho[i] = c33 * cons.rho[i] + c66 * cons1.rho[i] + dt66 * l2.rho[i];
        con2_mx1[i] = c33 * cons.mx1[i] + c66 * cons1.mx1[i] + dt66 * l2.mx1[i];
        con2_mx2[i] = c33 * cons.mx2[i] + c66 * cons1.mx2[i] + dt66 * l2.mx2[i];
        con2_mx3[i] = c33 * cons.mx3[i] + c66 * cons1.mx3[i] + dt66 * l2.mx3[i];
        con2_bx1[i] = c33 * cons.bx1[i] + c66 * cons1.bx1[i] + dt66 * l2.bx1[i];
        con2_bx2[i] = c33 * cons.bx2[i] + c66 * cons1.bx2[i] + dt66 * l2.bx2[i];
        con2_bx3[i] = c33 * cons.bx3[i] + c66 * cons1.bx3[i] + dt66 * l2.bx3[i];
        con2_e[i] = c33 * cons.e[i] + c66 * cons1.e[i] + dt66 * l2.e[i];
    }

    let cons2 = Conserved{
        rho: con2_rho,
        mx1: con2_mx1, 
        mx2: con2_mx2,
        mx3: con2_mx3,
        bx1: con2_bx1,
        bx2: con2_bx2,
        bx3: con2_bx3,
        e: con2_e,
    };

    cons2
}

//////////////////////////////
// Two Dimensional Functions
//////////////////////////////