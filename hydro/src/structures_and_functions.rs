/*
This file contains custom data types for this hydrocode. It also contains some functions used to run the 
simulations.

Author: Brayden JoHantgen
Last Update: 9/7/2026
*/

////////////
// Imports
////////////
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::f64::consts::PI;

//use crate::

/////////////////////
// Defining Structs
/////////////////////
pub struct Init_Conds {
    pub adiabatic_index: f64,
    pub discontinuity: f64,
    pub p: (f64, f64),
    pub rho: (f64, f64),
    pub vx1: (f64, f64),
    pub vx2: (f64, f64),
    pub vx3: (f64, f64),
    pub bx1: (f64, f64),
    pub bx2: (f64, f64),
    pub bx3: (f64, f64),
}

pub struct Driver {
    pub cfl: f64,
    pub t_final: f64,
    pub checkpoint: f64,
    pub x1_cell: usize,
    pub x1_range: (f64, f64),
    pub x2_cell: usize,
    pub x2_range: (f64, f64),
    pub x3_cell: usize,
    pub x3_range: (f64, f64),
}

pub struct Primatives {
    pub p: Vec<f64>,
    pub rho: Vec<f64>,
    pub vx1: Vec<f64>,
    pub vx2: Vec<f64>,
    pub vx3: Vec<f64>,
    pub bx1: Vec<f64>,
    pub bx2: Vec<f64>,
    pub bx3: Vec<f64>,
}

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

///////////////////////
// Defining Functions
///////////////////////
/// Input:
/// Output:
/// Description:
pub fn write_checkpoint(prims: &Primitives, drive: &Driver, t: f64, check_count: i64) -> Result<(), Box<dyn std::error::Error>> {
    let file_num = check_count.to_string();
    let file_type = ".txt".to_string();
    let file_name = format!("{}{}", file_num, file_type);
    let file_path = "time_step_files/".to_string() + &file_name;
    let output_file_path = Path::new(&file_path);
    
    let t_fill = format!("{} {}", "t:".to_string(), &(t.to_string()+&" ".to_string()));
    
    let size = prims.p.len();

    let num_fill_x = format!("{} {}", "cell_num_x:".to_string(), &(drive.x_cell.to_string()+&" ".to_string()));
    let left_x_fill = format!("{} {}", "left_x_bound:".to_string(), &(drive.x_range.0.to_string()+&" ".to_string()));
    let right_x_fill = format!("{} {}", "right_x_bound:".to_string(), &(drive.x_range.1.to_string()+&" ".to_string()));

    let num_fill_y = format!("{} {}", "cell_num_y:".to_string(), &(drive.y_cell.to_string()+&" ".to_string()));
    let left_y_fill = format!("{} {}", "left_y_bound:".to_string(), &(drive.y_range.0.to_string()+&" ".to_string()));
    let right_y_fill = format!("{} {}", "right_y_bound:".to_string(), &(drive.y_range.1.to_string()+&" ".to_string()));
    
    let mut p_string = "p: ".to_string();
    let mut rho_string = "rho: ".to_string();
    let mut vx_string = "vx: ".to_string();
    let mut vy_string = "vy: ".to_string();
    let mut vz_string = "vz: ".to_string();
    let mut bx_string = "Bx: ".to_string();
    let mut by_string = "By: ".to_string();
    let mut bz_string = "Bz: ".to_string();
        
    for i in 0..size {
        let p_fill = prims.p[i].to_string() + &" ".to_string();
        p_string.push_str(&p_fill);
    
        let rho_fill = prims.rho[i].to_string() + &" ".to_string();
        rho_string.push_str(&rho_fill);
    
        let vx_fill = prims.vx[i].to_string() + &" ".to_string();
        vx_string.push_str(&vx_fill);
    
        let vy_fill = prims.vy[i].to_string() + &" ".to_string();
        vy_string.push_str(&vy_fill);
    
        let vz_fill = prims.vz[i].to_string() + &" ".to_string();
        vz_string.push_str(&vz_fill);
    
        let bx_fill = prims.bx[i].to_string() + &" ".to_string();
        bx_string.push_str(&bx_fill);
    
        let by_fill = prims.by[i].to_string() + &" ".to_string();
        by_string.push_str(&by_fill);
    
        let bz_fill = prims.bz[i].to_string() + &" ".to_string();
        bz_string.push_str(&bz_fill);
    }
     
    let file = fs::OpenOptions::new()
        .write(true)
        .create(true)
        .open(output_file_path)?;
    
    let mut writer = BufWriter::new(file);
    
    writeln!(writer, "{}", num_fill_x)?;
    writeln!(writer, "{}", num_fill_y)?;
    writeln!(writer, "{}", left_x_fill)?;
    writeln!(writer, "{}", right_x_fill)?;
    writeln!(writer, "{}", left_y_fill)?;
    writeln!(writer, "{}", right_y_fill)?;

    writeln!(writer, "{}", t_fill)?;

    writeln!(writer, "{}", p_string)?;
    writeln!(writer, "{}", rho_string)?;
    writeln!(writer, "{}", vx_string)?;
    writeln!(writer, "{}", vy_string)?;
    writeln!(writer, "{}", vz_string)?;
    writeln!(writer, "{}", bx_string)?;
    writeln!(writer, "{}", by_string)?;
    writeln!(writer, "{}", bz_string)?;

    writer.flush()?;
    
    Ok(())
}