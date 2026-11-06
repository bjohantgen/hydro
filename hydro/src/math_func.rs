// This file is full of functions to supplement the "hydro" code.
//
// Author: Brayden JoHantgen
// Last Update: 11/6/2026

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

////////////////////////////////
// Defining Physical Functions
////////////////////////////////

/// Input:
///     prim (Cell8): the primitive variables of a single cell
///     a_index (f64): the adiabatic index of the gas
/// Output:
///     a (f64): the sound speed of the gas
/// Description:
///     Equation: c_s = sqrt(gamma * P / rho)
///     This function calculates the sound speed using the equation above.
///     This function does not check for negative arguments of the square root.
pub fn sound_speed(prim: Cell8, a_index: f64) -> f64 {
    let a = (a_index * prim.0 / prim.1).sqrt();
    a
}

/// Input:
///     m (f64): the black hole mass
///     a_index (f64): the adiabatic index of the gas
///     p_infinity (f64): the pressure of the gas at infinity
///     rho_infinity (f64): the density of the gas at infinity
/// Output:
///     rb (f64): the Bondi radius of the black hole with mass m
/// Description:
///     Equation: R_b = G * M / (c_s)^2
///     This function calculates the Bondi radius of a black hole and assumes G = 1.
///     This function does not check for negative arguments of the square root.
pub fn bondi_radius(m: f64, a_index: f64, p_infinity: f64, rho_infinity: f64) -> f64 {
    let c_s = (a_index * p_infinity / rho_infinity).sqrt();
    let rb = m / (f64::powf(c_s, 2.0));
    rb
}

/// Input:
///     m (f64): the black hole mass
///     a_index (f64): the adiabatic index of the gas
///     p_infinity (f64): the pressure of the gas at infinity
///     rho_infinity (f64): the density of the gas at infinity
/// Output:
///     m_dot (f64): the Bondi accretion rate
/// Description:
///     Equations: M_dot = 4 * pi * lamda_c * (G * M)^2 * c_s^-3 * rho_infinity
///                c_s = sqrt(gamma * p_infinity / rho_infinity)
///                lamda_c = (1/2)^((gamma + 1) / (2 * (gamma - 1))) * ((5 - 3 * gamma) / 4)^(-(5 - 3 * gamma) / (2 * (gamma - 1)))
///     This function calculates the Bondi accretion rate and assumes G = 1.
///     This function does not check for negative arguments of the square root.
pub fn bondi_accretion_rate(m: f64, a_index: f64, p_infinity: f64, rho_infinity: f64) -> f64 {
    let l_c = f64::powf(0.5, (a_index + 1.0) / (2.0 * (a_index - 1.0))) * f64::powf((5.0 - 3.0 * a_index) / 4.0, -(5.0 - 3.0 * a_index) / (2.0 * (a_index - 1.0)));
    let c_s = (a_index * p_infinity / rho_infinity).sqrt();
    let m_dot = 4.0 * PI * l_c * f64::powf(m, 2.0) * rho_infinity / (f64::powf(c_s, 3.0));
    m_dot
}

/// Input:
///     m (f64): the black hole mass
///     a_index (f64): the adiabatic index of the gas
///     p_infinity (f64): the pressure of the gas at infinity
///     rho_infinity (f64): the density of the gas at infinity
/// Output:
///     t (f64): the characteristic time of the simulation
/// Description:
///     Equation: t = r_b / c_infinity
///     This function uses the Bondi radius function to calculate the characteristic
///     time of the simulation.
pub fn time_unit(m: f64, a_index: f64, p_infinity: f64, rho_infinity: f64) -> f64 {
    let r_b = bondi_radius(m, a_index, p_infinity, rho_infinity);
    let c_s = (a_index * p_infinity / rho_infinity).sqrt();
    let t = r_b / c_s;
    t
}

/// Input:
/// Output:
/// Description:
pub fn radial_distance(i_num: usize, dr: f64) -> f64 {
    let r = (i_num as f64) * dr;
    r 
}

/// Input:
/// Output:
/// Description:
pub fn gravity(m:f64, r: f64, e: f64, dx1: f64) -> f64 {
    let numerator = -1.0 * m;
    let denominator = f64::powf(r, 2.0) + (f64::powf(e, 2.0))* (f64::powf(dx1, 2.0)) * f64::exp(-f64::powf(r, 2.0) / (f64::powf(e, 2.0) * f64::powf(dx1, 2.0)));
    let p = numerator / (denominator.sqrt());
    p
}

/// Input:
///     prim (Cell8): the primitive variables of the cell
///     a_index (f64): the adiabatic index
///     i_num (usize): the index of the cell
///     dx1 (f64): the width of a cell
///     m (f64): the mass of the black hole
///     e (f64): the softening parameter
/// Output:
///     tot_energy (f64): the total energy density
/// Description:
///     Equation: E = rho * (u^2/2 + epsilon + phi)
///     This function uses the radial distance and gravity functions to calculate
///     the toal energy density of the desired cell. 
pub fn total_energy_cart(prim: Cell8, a_index: f64, i_num: usize, dx1: f64, m: f64, e: f64) -> f64 {
    let energy = 0.5 * prim.1 * (prim.2 * prim.2 + prim.3 * prim.3 + prim.4 * prim.4) + prim.0 / (a_index - 1.0) + 0.5 * (prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7);
    let position = radial_distance(i_num, dx1);
    let grav = gravity(m, position, e, dx1);
    let tot_energy = energy + prim.1 * grav;
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
pub fn source(prim:Cell8, i_num: usize, dx1: f64, m: f64, e: f64) -> Cell8 {
    let r = radial_distance(i_num, dx1);
    let r_soft = (r*r + e*e*dx1*dx1).sqrt();
    let grav = gravity(m, r, e, dx1);
    let s0 = 0.0;
    let s1 = (-2.0 / r_soft) * prim.0 + prim.1 * grav;
    let s2 = (-2.0 / r_soft) * prim.0 + prim.2 * grav;
    let s3 = (-2.0 / r_soft) * prim.0 + prim.3 * grav;
    let s7 = prim.1 * prim.2 * grav;
    let s = (s0, s1, s2, s3, 0.0, 0.0, 0.0, s7);
    s
}

//////////////////////////////////
// Defining Simulation Functions
//////////////////////////////////

/// Input:
///     tup (Cell3): a tuple with three elements
/// Output:
///     max_check (f64): the maximum element of the tuple
/// Description:
///     This function compares the elements of a tuple with negative infinity
///     and stores the variable if it is larger than the previous. The function
///     checks for NaN elements of the tuple.
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
///     tup (Cell3): a tuple with three elements
/// Output:
///     min_check (f64): the minimum element of the tuple
/// Description:
///     This function compares the elements of a tuple with infinity
///     and stores the variable if it is smaller than the previous. The 
///     function checks for NaN elements of the tuple.
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
///     num (f64): any number
/// Output:
///     sign (f64): either 1.0, -1.0, or 0.0 depending on the 
///                 sign of the number
/// Description:
///     If the input number equals 0 then the out put is 0. If the 
///     input number is positive then the number equals 1. If the
///     input number is negative then the output is -1.
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
///     x (f64): argument 1
///     y (f64): argument 2
///     z (f64): argument 3
/// Output:
///     mm (f64): the value calculated by the minmod funciton
/// Description:
///     Equation: minmod(x,y,z) = 0.25 * |sgn(x) + sgn(y)| * (sgn(x) + sgn(z)) * min(|x|, |y|, |z|)
///     This equation uses the sgn function and the tuple_min function to calculate the minmod function
///     of a minmod slope limiter.
pub fn minmod(x: f64, y: f64, z: f64) -> f64 {
    let mm_1 = (sgn(x) + sgn(y)).abs();
    let mm_2 = sgn(x) + sgn(z);
    let mm_3 = tuple_min((x.abs(), y.abs(), z.abs()));
    let mm = 0.25 * mm_1 * mm_2 * mm_3;
    mm
}

/// Input:
///     ci (f64): desired prim of the ith cell on the grid
///     cip1 (f64): desired prim of the ith+1 cell on the grid
///     cip2 (f64): desired prim of the ith+2 cell on the grid
///     cim1 (f64): desired prim of the ith-1 cell on the grid
///     left (bool): if this is true the left-biased interface value of the left state will be reconstructed.
/// Output:
///     a (f64): the reconstructed prim value
/// Description:
///     This function uses the minmod function to reconstruct the primitive variables
///     of multiple cells using the piecewise linear method in order to get higher order
///     spacial calculations.
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
///     prim_l (Cell8): the primitive variables for the left cell
///     prim_r (Cell8): the primitive variables for the right cell
///     a_index (f64): the adiabatic index
///     dx1 (f64): the width of a cell
/// Output:
///     dt (f64): the timestep allowed for those cells
/// Description:
///     This function uses the sound_speed and tuple_max functions to determine
///     the timestep satisfied by the Courant-Friedrich-Levy condition.
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

/// Input:
/// Output:
/// Description:
pub fn gravity_time_step(m: f64, r: f64, e: f64, dx1: f64) -> f64 {
    let g = gravity(m, r, e, dx1).abs();
    if g > 0.0 {
        (dx1 / g).sqrt()
    } else {
        f64::INFINITY
    }
}