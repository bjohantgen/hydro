// This is the main file of the hydro simulation.
//
// Author: Brayden JoHantgen
// Last Update: 9/8/2026

////////////
// Imports
////////////
#![allow(unused_imports)]
use std::time::Instant;
use std::f64::consts::PI;

pub mod math_func;
pub mod sim_func;




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

///////////////
// Simulation
///////////////
fn main() {

// Initial Conditions
    
    let init_conds = sim_func::Conditions{
                        adiabatic_index: 1.4, 
                        p: (2.5, 2.5), 
                        rho: (2.0, 1.0), 
                        vx: (0.5, -0.5), 
                        vy: (0.0, 0.0), 
                        vz: (0.0, 0.0), 
                        bx: ((4.0*PI).sqrt()*(0.5), (4.0*PI).sqrt()*(0.5)), 
                        by: (0.0, 0.0), 
                        bz: (0.0, 0.0)};

    let drive = sim_func::Driver{
                        cfl: 0.3, 
                        tfinal: 1.001, 
                        checkpoint: 0.0125, 
                        x_cell: 5,
                        x_range: (-0.5, 0.5),
                        y_cell: 5, 
                        y_range: (-0.5, 0.5), 
                        discontinuity: 0.25};

// Simulation
    let onedimension: bool;
        if drive.y_cell == 0 {
            onedimension = true;
        } else {
            onedimension = false;
        }

    let before = Instant::now();

    let mut t: f64 = 0.0; 
    let mut t_checkpoint = drive.checkpoint;
    let mut time_step_count: f64 = 0.0;
    let mut check_count: i64 = 0;

    let onedimension: bool;
    if drive.y_cell == 0 {
        onedimension = true;
    } else {
        onedimension = false;
    }

    let before = Instant::now();

    let mut t: f64 = 0.0; 
    let mut t_checkpoint = drive.checkpoint;
    let mut time_step_count: f64 = 0.0;
    let mut check_count: i64 = 0;

    if onedimension == true {
        let dx = (drive.x_range.1 - drive.x_range.0) / (drive.x_cell as f64);

        let mut prims = sim_func::init_prims_1d(&init_conds, &drive);
        //let test_1 = sim_func::cons_from_prim_1d(&prims, init_conds.adiabatic_index);
        //let test_2 = sim_func::prim_from_cons_1d(&test_1, init_conds.adiabatic_index);
        let mut cons: sim_func::Conserved;
        //println!("{:?}", test_2);
        
        let _ = sim_func::write_checkpoint(&prims, &drive, t, check_count);
        check_count +=1;

        while t < drive.tfinal {
            cons = sim_func::cons_from_prim_1d(&prims, init_conds.adiabatic_index);

            let mut dt = f64::INFINITY;
            for i in 0..(drive.x_cell-1) {
                let prim_fill_1 = (prims.p[i], prims.rho[i], prims.vx[i], prims.vy[i], prims.vz[i], prims.bx[i], prims.by[i], prims.bz[i]);
                let prim_fill_2 = (prims.p[i+1], prims.rho[i+1], prims.vx[i+1], prims.vy[i+1], prims.vz[i+1], prims.bx[i+1], prims.by[i+1], prims.bz[i+1]);
                let dt_check = math_func::compute_time_step_x(prim_fill_1, prim_fill_2, init_conds.adiabatic_index, dx);
                dt = dt.min(dt_check);
            }
            dt = drive.cfl * dt;

            cons = sim_func::rk_step_1d(&prims, &cons, init_conds.adiabatic_index, dx, dt);
            prims = sim_func::prim_from_cons_1d(&cons, init_conds.adiabatic_index);

            if t >= t_checkpoint {
                let _ = sim_func::write_checkpoint(&prims, &drive, t, check_count);
                t_checkpoint += drive.checkpoint;
                check_count += 1;
                }
            t += dt;
            time_step_count += 1.0;
            
            println!(
                "step={} t={} dt={}",
                time_step_count,
                t,
                dt
                );
        }
        let runtime = (before.elapsed().as_millis() as f64) / 1000.0;
        let performance = (time_step_count * (prims.p.len() as f64)) / runtime;
        println!("The number of timesteps: {:?}", time_step_count);
        println!("The runtime in seconds: {:?}", runtime);
        println!("The performance in zones per second: {:.2?}", performance);
    } else {
        let dx = (drive.x_range.1 - drive.x_range.0) / (drive.x_cell as f64);
        let dy = (drive.y_range.1 - drive.y_range.0) / (drive.y_cell as f64);

        let mut prims = sim_func::init_prims_2d(&init_conds, &drive);
        println!("{:?}", prims);
    }
}
/*
    } else {
        let dx = (drive.x_range.1 - drive.x_range.0) / (drive.x_cell as f64);
        let dy = (drive.y_range.1 - drive.y_range.0) / (drive.y_cell as f64);

        let mut prims = sim_func::init_prims_2d(&init_conds, &drive);
        let mut cons: sim_func::Conserved;

        let _ = sim_func::write_checkpoint(&prims, &drive, t, check_count);
        check_count +=1;

        while t < drive.tfinal {
            cons = sim_func::cons_from_prim_2d(&prims, drive.x_cell, drive.y_cell, init_conds.adiabatic_index);
        
            let mut dt = f64::INFINITY;
            for i in 0..drive.x_cell {
                for j in 0..drive.y_cell {
                    let ip1 = (i + 1) % drive.x_cell;
                    let jp1 = (j + 1) % drive.y_cell;

                    let ind_a = sim_func::index(i,j,drive.y_cell);
                    let ind_b = sim_func::index(ip1,j,drive.y_cell);
                    let ind_c = sim_func::index(i,j,drive.y_cell);
                    let ind_d = sim_func::index(i,jp1,drive.y_cell);

                    let prim_fill_1 = (prims.p[ind_a], prims.rho[ind_a], prims.vx[ind_a], prims.vy[ind_a], prims.vz[ind_a], prims.bx[ind_a], prims.by[ind_a], prims.bz[ind_a]);
                    let prim_fill_2 = (prims.p[ind_b], prims.rho[ind_b], prims.vx[ind_b], prims.vy[ind_b], prims.vz[ind_b], prims.bx[ind_b], prims.by[ind_b], prims.bz[ind_b]);
                    
                    let prim_fill_3 = (prims.p[ind_c], prims.rho[ind_c], prims.vx[ind_c], prims.vy[ind_c], prims.vz[ind_c], prims.bx[ind_c], prims.by[ind_c], prims.bz[ind_c]);
                    let prim_fill_4 = (prims.p[ind_d], prims.rho[ind_d], prims.vx[ind_d], prims.vy[ind_d], prims.vz[ind_d], prims.bx[ind_d], prims.by[ind_d], prims.bz[ind_d]);

                    let dt_check1 = math_func::compute_time_step_x(prim_fill_1, prim_fill_2, init_conds.adiabatic_index, dx);
                    let dt_check2 = math_func::compute_time_step_y(prim_fill_3, prim_fill_4, init_conds.adiabatic_index, dy);
                    dt = dt.min(dt_check1);
                    dt = dt.min(dt_check2);
                }
            }
            //if dt < 1e-8 {
            //    dt = 1e-8
            //}
            dt = drive.cfl * dt;  
        
            cons = sim_func::rk_step_2d(&prims, &cons, drive.x_cell, drive.y_cell, init_conds.adiabatic_index, dx, dy, dt);

            let div_b = sim_func::max_div_b(&cons.bx, &cons.by, drive.x_cell, drive.y_cell, dx, dy);

            println!("step={} t={} dt={} max_div_b={}", time_step_count, t, dt, div_b);

            prims = sim_func::prim_from_cons_2d(&cons, drive.x_cell, drive.y_cell, init_conds.adiabatic_index);

            if t >= t_checkpoint {
                let _ = sim_func::write_checkpoint(&prims, &drive, t, check_count);
                t_checkpoint += drive.checkpoint;
                check_count += 1;
            }
            t += dt;
            time_step_count += 1.0;

            if time_step_count >= 100.0 {
                break;
            }

        }
        let runtime = (before.elapsed().as_millis() as f64) / 1000.0;
        let performance = (time_step_count * (prims.p.len() as f64)) / runtime;
        println!("The number of timesteps: {:?}", time_step_count);
        println!("The runtime in seconds: {:?}", runtime);
        println!("The performance in zones per second: {:.2?}", performance);
    }
}
*/