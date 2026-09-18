// This file is full of functions to supplement the "hydro" code.
//
// Author: Brayden JoHantgen
// Last Update: 9/17/2026

////////////
// Imports
////////////
use std::f64::consts::PI;

/////////////////////
// Define Datatypes
/////////////////////
pub type Cell8 = (f64, f64, f64, f64, f64, f64, f64, f64);
pub type Cell5 = (f64, f64, f64, f64, f64);
pub type Cell3 = (f64, f64, f64);

////////////////////////////////////////////////
// Defining Coordinate Free Physical Functions
////////////////////////////////////////////////

/// Input:
/// Output:
/// Description:
pub fn sound_speed(prim: Cell8, a_index: f64) -> f64 {
    let a = (a_index * prim.0 / prim.1).sqrt();
    a
}

/// Input:
/// Output:
/// Description:
pub fn bondi_accretion_rate(m: f64, a_index: f64, p_infinity: f64, rho_infinity: f64) -> f64 {
    let l_c = f64::powf(0.5, (a_index + 1.0) / (2.0 * (a_index - 1.0))) * f64::powf((5.0 - 3.0 * a_index) / 4.0, -(5.0 - 3.0 * a_index) / (2.0 * (a_index - 1.0)));
    let c_s = (a_index * p_infinity / rho_infinity).sqrt();
    let m_dot = 4.0 * PI * l_c * f64::powf(m, 2.0) * rho_infinity / (f64::powf(c_s, 3.0));
    m_dot
}

/// Input:
/// Output:
/// Description:
pub fn gravity(m:f64, r: f64, e: f64, dx1: f64) -> f64 {
    let numerator = -1.0 * m;
    let denominator = f64::powf(r, 2.0) + (f64::powf(e, 2.0)) * (f64::powf(dx1, 2.0)) * f64::exp(-f64::powf(r, 2.0) / (f64::powf(e, 2.0) * f64::powf(dx1, 2.0)));
    let p = numerator / (denominator.sqrt());
    p
}

//////////////////////////////////////////
// Defining Cartesean Physical Functions
//////////////////////////////////////////

/// Input:
/// Output:
/// Description:
pub fn total_pressure_cart(prim: Cell8) -> f64 {
    let p = prim.0 + 0.5 * (prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7);
    p
}

/// Input:
/// Output:
/// Description:
pub fn total_energy_cart(prim: Cell8, a_index: f64, i_num: usize, dx1: f64, m: f64, e: f64) -> f64 {
    let e = 0.5 * prim.1 * (prim.2 * prim.2 + prim.3 * prim.3 + prim.4 * prim.4) + prim.0 / (a_index - 1.0) + 0.5 * (prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7);
    let position = radial_distance(i_num, dx1);
    let grav = gravity(m, position, e, dx1);
    let tot_energy = e + grav;
    tot_energy
}

/// Input:
/// Output:
/// Description:
pub fn flux_x(prim: Cell8, a_index: f64, i_num: usize, dx1: f64, m: f64, e: f64) -> Cell8 {
    let p = total_pressure_cart(prim);
    let energy = total_energy_cart(prim, a_index, i_num, dx1, m, e);
    let f0 = prim.1 * prim.2;
    let f1 = prim.1 * prim.2 * prim.2 + p - prim.5 * prim.5;
    let f2 = prim.1 * prim.2 * prim.3 - prim.5 * prim.6;
    let f3 = prim.1 * prim.2 * prim.4 - prim.5 * prim.7;
    let f4 = 0.0;
    let f5 = prim.6 * prim.2 - prim.5 * prim.3;
    let f6 = prim.7 * prim.2 - prim.5 * prim.4;  
    let f7 = (energy + p) * prim.2 - prim.5 * (prim.2 * prim.5 + prim.3 * prim.6 + prim.4 * prim.7);
    let f = (f0, f1, f2, f3, f4, f5, f6, f7);
    f
}

/// Input:
/// Output:
/// Description:
pub fn flux_y(prim:Cell8, a_index: f64, i_num: usize, dx1: f64, m: f64, e: f64) -> Cell8 {
    let p = total_pressure_cart(prim);
    let energy = total_energy_cart(prim, a_index, i_num, dx1, m, e);
    let f0 = prim.1 * prim.3;
    let f1 = prim.1 * prim.2 * prim.3 - prim.5 * prim.6;
    let f2 = prim.1 * prim.3 * prim.3 + p - prim.6 * prim.6;
    let f3 = prim.1 * prim.3 * prim.4 - prim.6 * prim.7;
    let f4 = prim.5 * prim.3 - prim.6 * prim.2;
    let f5 = 0.0;
    let f6 = prim.7 * prim.3 - prim.6 * prim.4;
    let f7 = (energy + p) * prim.3 - prim.6 * (prim.2 * prim.5 + prim.3 * prim.6 + prim.4 * prim.7);
    let f_y = (f0, f1, f2, f3, f4, f5, f6, f7);
    f_y
}

//////////////////////////////////////////
// Defining Spherical Physical Functions
//////////////////////////////////////////
pub fn radial_distance(i_num: usize, dr: f64) -> f64 {
    let r = (i_num as f64) * dr;
    r 
}

pub fn source(prim:Cell8, a_index: f64, i_num: usize, dx1: f64, m: f64, e: f64) -> Cell8 {
    let r = radial_distance(i_num, dx1);
    let energy = total_energy_cart(prim, a_index, i_num, dx1, m, e);
    let s0 = (-2.0 / r) * prim.1 * prim.2;
    let s1 = (-2.0 / r) * prim.1 * f64::powf(prim.2, 2.0);
    let s2 = (-2.0 / r) * prim.2 * (energy + prim.0);
    let s = (s0, s1, s2, 0.0, 0.0, 0.0, 0.0, 0.0);
    s
}

//////////////////////////////////
// Defining Simulation Functions
//////////////////////////////////

/// Input:
/// Output:
/// Description:
pub fn tuple_max(tup: Cell3) -> f64 {
    let arr = [tup.0, tup.1, tup.2];
    let mut max_check: f64 = f64::NEG_INFINITY;
    for i in 0..3 {
        debug_assert!(!arr[i].is_nan(), "tuple_max received NaN input");
        if arr[i] > max_check {
            max_check = arr[i];
        }
    }
    max_check
}

/// Input:
/// Output:
/// Description:
pub fn tuple_min(tup: Cell3) -> f64 {
    let arr = [tup.0, tup.1, tup.2];
    let mut min_check: f64 = f64::INFINITY;
    for i in 0..3 {
        debug_assert!(!arr[i].is_nan(), "tuple_min received NaN input");
        if arr[i] < min_check {
            min_check = arr[i]
        }
    }
    min_check
}

/// Input:
/// Output:
/// Description:
pub fn sgn(num: f64) -> f64 {
    let sign: f64;
    if num > 0.0 {
        sign = 1.0;
    } else if num < 0.0 {
        sign = -1.0;
    }
    else {
        sign = 0.0;
    }
    sign
}

/// Input:
/// Output:
/// Description:
pub fn minmod(x: f64, y: f64, z: f64) -> f64 {
    let mm_1 = (sgn(x) + sgn(y)).abs();
    let mm_2 = sgn(x) + sgn(z);
    let mm_3 = tuple_min((x.abs(), y.abs(), z.abs()));
    let mm = 0.25 * mm_1 * mm_2 * mm_3;
    mm
}

/// Input:
/// Output:
/// Description:
pub fn plm_reconstruction(ci: f64, cip1: f64, cip2: f64, cim1: f64, left: bool) -> f64 {
    let a: f64;
    if left == true {
        let num_1 = 1.5 * (ci - cim1);
        let num_2 = 0.5 * (cip1 - cim1);
        let num_3 = 1.5 * (cip1 - ci);
        a = ci + 0.5 * minmod(num_1, num_2, num_3);
    } else {
        let num_1 = 1.5 * (cip1 - ci);
        let num_2 = 0.5 * (cip2 - ci);
        let num_3 = 1.5 * (cip2 - cip1);
        a = cip1 - 0.5 * minmod(num_1, num_2, num_3);
    }
    a
}

/// Input:
/// Output:
/// Description:
pub fn compute_time_step_1d(prim_l: Cell8, prim_r: Cell8, a_index: f64, dx1: f64) -> f64 {
    let cf_l = sound_speed(prim_l, a_index);
    let plus_l = cf_l + prim_l.2;
    let minus_l = prim_l.2 - cf_l;

    let cf_r = sound_speed(prim_r, a_index); 
    let plus_r = cf_r + prim_r.2; 
    let minus_r = prim_r.2 - cf_r;

    let a_plus = tuple_max((0.0, plus_l, plus_r));
    let a_minus = tuple_max((0.0, -minus_l, -minus_r));

    let mut dt: f64 = 0.0;
    
    if a_minus > a_plus {
        dt += dx1 / a_minus;
    } else {
        dt += dx1 / a_plus;
    }    
    dt
}