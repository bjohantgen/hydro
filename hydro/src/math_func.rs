// This file is full of functions to supplement the "hydro" code.
//
// Author: Brayden JoHantgen
// Last Update: 8/4/2026

//////////////////////////
// Initialized Datatypes
//////////////////////////
pub type Cell8 = (f64, f64, f64, f64, f64, f64, f64, f64);
pub type Cell5 = (f64, f64, f64, f64, f64);
pub type Cell3 = (f64, f64, f64);

//////////////////////////////////////////
// Defining Cartesean Physical Functions
//////////////////////////////////////////
#[inline(always)]
pub fn total_pressure(prim: Cell8) -> f64 {
    let p = prim.0 + 0.5 * (prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7);
    p
}

#[inline(always)]
pub fn total_energy(prim: Cell8, a_index: f64) -> f64 {
    let e = 0.5 * prim.1 * (prim.2 * prim.2 + prim.3 * prim.3 + prim.4 * prim.4) + prim.0 / (a_index - 1.0) + 0.5 * (prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7);
    e
}

pub fn flux_x(prim: Cell8, a_index: f64) -> Cell8 {
    let p = total_pressure(prim);
    let e = total_energy(prim, a_index);
    let f0 = prim.1 * prim.2;
    let f1 = prim.1 * prim.2 * prim.2 + p - prim.5 * prim.5;
    let f2 = prim.1 * prim.2 * prim.3 - prim.5 * prim.6;
    let f3 = prim.1 * prim.2 * prim.4 - prim.5 * prim.7;
    let f4 = 0.0;
    let f5 = prim.6 * prim.2 - prim.5 * prim.3;
    let f6 = prim.7 * prim.2 - prim.5 * prim.4;  
    let f7 = (e + p) * prim.2 - prim.5 * (prim.2 * prim.5 + prim.3 * prim.6 + prim.4 * prim.7);
    let f = (f0, f1, f2, f3, f4, f5, f6, f7);
    f
}

pub fn flux_y(prim:Cell8, a_index: f64) -> Cell8 {
    let p = total_pressure(prim);
    let e = total_energy(prim, a_index);
    let f0 = prim.1 * prim.3;
    let f1 = prim.1 * prim.2 * prim.3 - prim.5 * prim.6;
    let f2 = prim.1 * prim.3 * prim.3 + p - prim.6 * prim.6;
    let f3 = prim.1 * prim.3 * prim.4 - prim.6 * prim.7;
    let f4 = prim.5 * prim.3 - prim.6 * prim.2;
    let f5 = 0.0;
    let f6 = prim.7 * prim.3 - prim.6 * prim.4;
    let f7 = (e + p) * prim.3 - prim.6 * (prim.2 * prim.5 + prim.3 * prim.6 + prim.4 * prim.7);
    let f_y = (f0, f1, f2, f3, f4, f5, f6, f7);
    f_y
}

#[inline(always)] 
pub fn sound_speed(prim: Cell8, a_index: f64) -> f64 {
    let a = (a_index * prim.0 / prim.1).sqrt();
    a
}

#[inline(always)]
pub fn fast_magsonic_speed_x(prim: Cell8, a_index: f64) -> f64 {
    let b_squared = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7;
    let a = a_index * prim.0;
    let disc = (a + b_squared) * (a + b_squared) - 4.0 * a * prim.5 * prim.5;
    let numerator = a + b_squared + disc.sqrt();
    let cfx = (numerator / (2.0 * prim.1)).sqrt();
    cfx
}

#[inline(always)]
pub fn fast_magsonic_speed_y(prim: Cell8, a_index: f64) -> f64 {
    let b_squared = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7;
    let a = a_index * prim.0;
    let disc = (a + b_squared) * (a + b_squared) - 4.0 * a * prim.6 * prim.6;
    let numerator = a + b_squared + disc.sqrt();
    let cfy = (numerator / (2.0 * prim.1)).sqrt();
    cfy
}

#[inline(always)]
pub fn fast_magsonic_speed_z(prim: Cell8, a_index: f64) -> f64 {
    let b_squared = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7;
    let a = a_index * prim.0;
    let disc = (a + b_squared) * (a + b_squared) - 4.0 * a * prim.7 * prim.7;
    let numerator = a + b_squared + disc.sqrt();
    let cfz = (numerator / (2.0 * prim.1)).sqrt();
    cfz
}

pub fn slow_magsonic_speed_x(prim: Cell8, a_index: f64) -> f64 {
    let b_squared = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7;
    let a = a_index * prim.0;
    let disc = (a + b_squared) * (a + b_squared) - 4.0 * a * prim.5 * prim.5;
    let numerator = a + b_squared - disc.sqrt();
    let csx = (numerator / (2.0 * prim.1)).sqrt();
    csx
}

pub fn slow_magsonic_speed_y(prim: Cell8, a_index: f64) -> f64 {
    let b_squared = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7;
    let a = a_index * prim.0;
    let disc = (a + b_squared) * (a + b_squared) - 4.0 * a * prim.6 * prim.6;
    let numerator = a + b_squared - disc.sqrt();
    let csy = (numerator / (2.0 * prim.1)).sqrt();
    csy
}

pub fn slow_magsonic_speed_z(prim: Cell8, a_index: f64) -> f64 {
    let b_squared = prim.5 * prim.5 + prim.6 * prim.6 + prim.7 * prim.7;
    let a = a_index * prim.0;
    let disc = (a + b_squared) * (a + b_squared) - 4.0 * a * prim.7 * prim.7;
    let numerator = a + b_squared - disc.sqrt();
    let csz = (numerator / (2.0 * prim.1)).sqrt();
    csz
}

pub fn alfven_speed_x(prim: Cell8) -> f64 {
    let cax = (prim.5.abs()) / (prim.1.sqrt());
    cax
}

pub fn alfven_speed_y(prim: Cell8) -> f64 {
    let cay = (prim.6.abs()) / (prim.1.sqrt());
    cay
}

pub fn alfven_speed_z(prim: Cell8) -> f64 {
    let caz = (prim.7.abs()) / (prim.1.sqrt());
    caz
}

pub fn compute_time_step_x(prim_l: Cell8, prim_r: Cell8, a_index: f64, dx: f64) -> f64 {
    let cf_l = fast_magsonic_speed_x(prim_l, a_index);
    let plus_l = cf_l + prim_l.2;
    let minus_l = prim_l.2 - cf_l;

    let cf_r = fast_magsonic_speed_x(prim_r, a_index); 
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

pub fn compute_time_step_y(prim_t: Cell8, prim_b: Cell8, a_index: f64, dy: f64) -> f64 {
    let cf_t = fast_magsonic_speed_y(prim_t, a_index);
    let plus_t = cf_t + prim_t.3;
    let minus_t = prim_t.3 - cf_t;

    let cf_b = fast_magsonic_speed_y(prim_b, a_index);
    let plus_b = cf_b + prim_b.3;
    let minus_b = prim_b.3 - cf_b;

    let a_plus = tuple_max((0.0, plus_t, plus_b));
    let a_minus = tuple_max((0.0, -minus_t, -minus_b));

    let mut dt: f64 = 0.0;
    
    if a_minus > a_plus {
        dt += dy / a_minus;
    } else {
        dt += dy / a_plus;
    }    
    dt
}

#[inline(always)]
pub fn emf(vx: f64, vy: f64, bx: f64, by: f64) -> f64 {
    let ez = vy * bx - vx * by;
    ez
}

#[inline(always)]
pub fn bx_center(bx_face: &[f64], i: usize, j: usize, ynum: usize) -> f64 {
    0.5 * (bx_face[i*ynum + j] + bx_face[(i+1)*ynum + j])
}

#[inline(always)]
pub fn by_center(by_face: &[f64], i: usize, j: usize, ynum: usize) -> f64 {
    0.5 * (by_face[i*(ynum+1) + j] + by_face[i*(ynum+1) + j + 1])
}

//////////////////////////////////
// Defining Simulation Functions
//////////////////////////////////
#[inline(always)]
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

#[inline(always)]
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

#[inline(always)]
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

#[inline(always)]
pub fn minmod(x: f64, y: f64, z: f64) -> f64 {
    let mm_1 = (sgn(x) + sgn(y)).abs();
    let mm_2 = sgn(x) + sgn(z);
    let mm_3 = tuple_min((x.abs(), y.abs(), z.abs()));
    let mm = 0.25 * mm_1 * mm_2 * mm_3;
    mm
}

#[inline(always)]
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

#[inline(always)]
pub fn ppm_reconstruction(cj: f64, cjp1: f64, cjp2: f64, cjm1: f64, cjm2: f64, left: bool) -> f64 {
    let a: f64;
    if left == true {
        a = (7.0 / 12.0) * (cj + cjp1) - (1.0 / 12.0) * (cjp2  + cjm1);
    } else {
        a = (7.0 / 12.0) * (cjm1 + cj) - (1.0 / 12.0) * (cjp1 + cjm2);
    }
    a
}

/// Input:
///     i: the index in the x direction 
///     j: the index in the y direction
///     ny: the length of the grid in the y direction
/// Output:
///     The flattened index of the 2D grid
/// Description:
///     This function has been inlined due to its frequent use. It flattens the 
///     two dimensional index to one dimension.
#[inline(always)]
pub fn index_2d(i: usize, j: usize, ynum: usize) -> usize {
    i * ynum + j
}