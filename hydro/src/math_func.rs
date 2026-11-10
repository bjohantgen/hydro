// This file is full of functions to supplement the "hydro" code.
//
// Author: Brayden JoHantgen
// Last Update: 11/9/2026

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
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
/// Output:
///     c (f64): the sound speed of the gas
/// Description:
///     Equation: c^2 = gamma * P / rho
///     This function calculates the sound speed of the gas.
///     It does no error checking.
pub fn sound_speed(prim: Cell5, a_ind: f64) -> f64 {
    let c = (a_ind * prim.0 / prim.1).sqrt();
    c
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
///     m_bh (f64): the mass of the black hole
/// Output:
///     rb (f64): the Bondi radius
/// Description:
///     Equation: Rb = G * M / c^2
///     This function uses the sound_speed function to calculate the 
///     Bondi radius. G = 1 and this function does no error checking.
pub fn bondi_radius(prim: Cell5, a_ind: f64, m_bh: f64) -> f64 {
    let c = sound_speed(prim, a_ind);

    let rb = m_bh / (f64::powf(c, 2.0));
    rb
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
///     m_bh (f64): the mass of the black hole
/// Output:
///     m_dot (f64): the Bondi accretion rate
/// Description:
///     Equation: M_dot = 4 * pi * lambda * Rb^2 * rho * c
///     This function uses the sound_speed and the bondi_radius functions
///     to calculate the Bondi accretion rate. This function does no error checking.
pub fn bondi_accretion(prim: Cell5, a_ind: f64, m_bh: f64) -> f64 {
    let c = sound_speed(prim, a_ind);
    let rb = bondi_radius(prim, a_ind, m_bh);

    let exponent_1 = (a_ind + 1.0) / (2.0 * (a_ind - 1.0));
    let exponent_2 = (3.0 * a_ind - 5.0) / (2.0 * (a_ind - 1.0));
    let lambda = (f64::powf(0.5, exponent_1)) * (f64::powf((5.0 - 3.0 * a_ind) / 4.0, exponent_2));

    let m_dot = 4.0 * PI * lambda * (f64::powf(rb, 2.0)) * prim.1 * c;
    m_dot
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
///     m_bh (f64): the mass of the black hole
/// Output:
///     bt (f64): the characteristic Bondi timescale
/// Description:
///     Equation: T_b = R_b / c
///     This function calculates the Bondi timescale. It does no
///     error checking.
pub fn bondi_time(prim: Cell5, a_ind: f64, m_bh: f64) -> f64 {
    let c = sound_speed(prim, a_ind);
    let rb = bondi_radius(prim, a_ind, m_bh);

    let bt = rb / c;
    bt
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
/// Output:
///     epsilon (f64): the internal energy density
/// Description:
///     Equation: epsilon = P / (rho * (gamma - 1))
///     This equation calcualtes the internal energy density and
///     does not check for errors.
pub fn internal_energy(prim: Cell5, a_ind: f64) -> f64 {
    let epsilon = prim.0 / (prim.1 * (a_ind - 1.0));
    epsilon
}

/// Input:
///     i_num (usize): the index number of the cell
///     delta (f64): the width of a cell
/// Output:
///     r (f64): the radial distance 
/// Description:
///     This function solves for the distance from the leftmost edge of
///     a one-dimensional grid. There is no error checking.
pub fn distance(i_num: usize, delta: f64) -> f64 {
    let r = (i_num as f64) * delta;
    r
}

/// Input:
///     m_bh (f64): the mass of the black hole
///     r (f64): the distance from the black hole
///     s (f64): the softening parameter
///     delta (f64): the width of a cell
/// Output:
///     phi (f64): the gravitational potential 
/// Description:
///     Equation: phi = - G * M / sqrt(r^2 + r_soft)
///               r_soft = s^2 * delta^2 * exp(- r^2 * s^-2 * delta^-2)
///     This function calculates the gravitational potential a distance r
///     away from the black hole. There is no error checking.
pub fn grav_potential(m_bh: f64, r: f64, s: f64, delta: f64) -> f64 {
    let r_soft = (f64::powf(s, 2.0)) * (f64::powf(delta, 2.0)) * f64::exp(-f64::powf(r, 2.0) / (f64::powf(s, 2.0) * f64::powf(delta, 2.0)));
    let phi = - m_bh / ((f64::powf(r, 2.0)) + r_soft).sqrt();
    phi
}

/// Input:
///     m_bh (f64): the mass of the black hole
///     r (f64): the distance from the black hole
///     s (f64): the softening parameter
///     delta (f64): the width of a cell
/// Output:
///     g (f64): the gravitational acceleration
/// Description:
///     Equation: g = - G * M * r * (1 - e) / (r^2 + s^2 * delta^2 * e)^(3/2)
///               e = exp(- r^2 * s^-2 * delta^-2)
///     This function calculates the gravitational acceleration. It does no 
///     error checking.
pub fn grav_acceleration(m_bh: f64, r: f64, s: f64, delta: f64) -> f64 {
    let soft_length = f64::powf(s * delta, 2.0);
    let exponent = f64::exp(- (f64::powf(r, 2.0)) / soft_length);
    let denominator = f64::powf(r, 2.0) + soft_length * exponent;
    let g = - m_bh * r * (1.0 - exponent) / (f64::powf(denominator, 1.5));
    g
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
///     m_bh (f64): the mass of the black hole
///     i_num (usize): the index of the cell
///     s (f64): the softening parameter
///     delta (f64): the width of a cell
/// Output:
///     energy (f64): the total energy of the cell
/// Description:
///     Equation: E = rho * (v^2/2 + epsilon + phi)
///     This functin uses the distance, grav_potential, and internal_energy functions
///     to calculate the total energy of the cell. This function does no error checking.
pub fn total_energy(prim: Cell5, a_ind: f64, m_bh: f64, i_num: usize, s: f64, delta: f64) -> f64 {
    let r = distance(i_num, delta);
    let phi = grav_potential(m_bh, r, s, delta);
    let epsilon = internal_energy(prim, a_ind);

    let v2 = prim.2 * prim.2 + prim.3 * prim.3 + prim.4 * prim.4;
    let energy = prim.1 * (0.5 * v2 + epsilon + phi);
    energy
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
/// Output:
///     m (f64): the mach number of the gas
/// Description:
///     Equation: M = v / c
///     This function uses the sound speed function to calculate
///     the mach number. There is no error checking.
pub fn mach_number(prim: Cell5, a_ind:f64) -> f64 {
    let c = sound_speed(prim, a_ind);
    
    let v = (prim.2 * prim.2 + prim.3 * prim.3 + prim.4 * prim.4).sqrt();
    let m = v / c;
    m
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     a_ind (f64): the adiabatic index of the gas
///     m_bh (f64): the mass of the black hole
///     i_num (f64): the index of the cell
///     s (f64): the softening parameter
///     delta (f64): the width of a cell
/// Output:
///     flux (Cell5): the flux of a cell
/// Description:
///     This function uses the total_energy function to calculate the flux.
///     There is no error checking.
pub fn flux(prim: Cell5, a_ind: f64, m_bh: f64, i_num: usize, s: f64, delta: f64) -> Cell5 {
    let e = total_energy(prim, a_ind, m_bh, i_num, s, delta);

    let f0 = prim.1 * prim.2;
    let f1 = prim.1 * prim.2 * prim.2 + prim.0;
    let f2 = prim.1 * prim.2 * prim.3;
    let f3 = prim.1 * prim.2 * prim.4;
    let f4 = prim.2 * (e + prim.0);

    let flux = (f0, f1, f2, f3, f4);
    flux
}

/// Input:
///     prim (Cell5): the primitive variables of a cell
///     m_bh (f64): the mass of the black hole
///     i_num (f64): the index of the cell
///     s (f64): the softening parameter
///     delta (f64): the width of a cell
/// Output:
///     source (f64): the source of the cell
/// Description:
///     This function uses the distance and grav_acceleration functions to calculate the source
///     terms. There is no error checking.
pub fn source(prim: Cell5, m_bh: f64, i_num: usize, s: f64, delta: f64) -> Cell5 {
    let r = distance(i_num, delta);
    let g = grav_acceleration(m_bh, r, s, delta);

    let s0 = 0.0;
    let s1 = (- 2.0 * prim.0 / r) + prim.1 * g;
    let s2 = 0.0;
    let s3 = 0.0;
    let s4 = prim.1 * prim.2 * g;

    let source = (s0, s1, s2, s3, s4);
    source
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
///     prim_l (Cell5): the primitive variables for the left cell
///     prim_r (Cell5): the primitive variables for the right cell
///     a_ind (f64): the adiabatic index
///     dx (f64): the width of a cell
/// Output:
///     dt (f64): the timestep allowed for those cells
/// Description:
///     This function uses the sound_speed and tuple_max functions to determine
///     the timestep satisfied by the Courant-Friedrich-Levy condition.
pub fn compute_time_step(prim_l: Cell5, prim_r: Cell5, a_ind: f64, dx: f64) -> f64 {
    let cf_l = sound_speed(prim_l, a_ind);

    let plus_l = cf_l + prim_l.2;
    let minus_l = prim_l.2 - cf_l;

    let cf_r = sound_speed(prim_r, a_ind); 

    let plus_r = cf_r + prim_r.2; 
    let minus_r = prim_r.2 - cf_r;

    let a_plus = tuple_max((0.0, plus_l, plus_r));
    let a_minus = tuple_max((0.0, -minus_l, -minus_r));

    let mut dt: f64 = 0.0;
    
    if a_minus > a_plus {
        dt += dx / a_minus;
    } else {
        dt += dx / a_plus;
    }    
    dt
}

/// Input:
///     m (f64): mass of the black hole
///     r (f64): the distance from the black hole
///     s (f64): the softness parameter
///     dr (f64): the size of a cell
/// Output:
///     the gravity timestep
/// Description:
///     This function uses the grav_acceleration to calculate the gravity
///     timestep. There is no error checking.
pub fn gravity_time_step(m: f64, r: f64, s: f64, dr: f64) -> f64 {
    let g = grav_acceleration(m, r, s, dr).abs();

    if g > 0.0 {
        (dr / g).sqrt()
    } else {
        f64::INFINITY
    }
}