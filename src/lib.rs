use core::f32;

use matrix::format::{Compressed, Conventional};

pub fn get_enter_var(obj_coef: &Compressed<f32>) -> (i32, f32) {
	let mut min_value = f32::INFINITY;
	let mut min_index = -1;

	for (i, &value) in obj_coef.values.iter().enumerate() {
		if value < min_value && value < 0.0 {
			min_value = value;
			min_index = i as i32;
		}
	}

	(min_index, min_value)
}

pub fn get_leaving_var(
	con_coef: &Compressed<f32>,
	rhs_coef: &Compressed<f32>,
	enter_idx: i32,
) -> (i32, f32) {
	let mut min_ratio = f32::INFINITY;
	let mut min_index = -1;

	for (i, &rhs) in rhs_coef.values.iter().enumerate() {
		if rhs > 0.0 {
			let row_coef = con_coef.get((i, enter_idx as usize));
			let ratio = rhs / row_coef;
			if ratio < min_ratio {
				if ratio == 0.0 && row_coef != 1.0 {
					continue;
				}

				min_ratio = ratio;
				min_index = i as i32;
			}
		}
	}

	(min_index, min_ratio)
}

pub fn pivot_coef(
	con_coef: &mut Compressed<f32>,
	con_rhs_coef: &mut Compressed<f32>,
	obj_coef: &mut Compressed<f32>,
	rhs_coef: &mut Compressed<f32>,
	enter_idx: i32,
	leaving_idx: i32,
) {
	let rhs_value = con_rhs_coef.get((leaving_idx as usize, 0));
	let filled_con_coef = Conventional::from(con_coef.clone());

	// get row values for the leaving variable
	let pivot_row = (0..filled_con_coef.columns)
		.map(|j| filled_con_coef[(leaving_idx as usize, j)])
		.collect::<Vec<_>>();
	let pivot_value = pivot_row[enter_idx as usize];

	for j in 0..con_coef.columns {
		if j != enter_idx as usize {
			let new_value = pivot_row[j] / pivot_value;
			con_coef.set((leaving_idx as usize, j), new_value);
		}
	}

	con_coef.set((leaving_idx as usize, enter_idx as usize), 1.0);
	// set corresponding rhs value in pivot row
	con_rhs_coef.set((leaving_idx as usize, 0), rhs_value / pivot_value);

	for i in 0..con_coef.rows {
		if i == leaving_idx as usize {
			continue;
		}

		let factor = con_coef.get((i, enter_idx as usize));
		for j in 0..con_coef.columns {
			let pivot_coef_val = con_coef.get((leaving_idx as usize, j));
			con_coef.set((i, j), con_coef.get((i, j)) - factor * pivot_coef_val);
		}

		let pivot_rhs_val = con_rhs_coef.get((leaving_idx as usize, 0));
		con_rhs_coef.set((i, 0), con_rhs_coef.get((i, 0)) - factor * pivot_rhs_val);
	}

	let obj_factor = obj_coef.get((0, enter_idx as usize));
	for j in 0..con_coef.columns {
		let pivot_coef_val = con_coef.get((leaving_idx as usize, j));
		obj_coef.set((0, j), obj_coef.get((0, j)) - obj_factor * pivot_coef_val);
	}

	let pivot_rhs_val = con_rhs_coef.get((leaving_idx as usize, 0));
	rhs_coef.set((0, 0), rhs_coef.get((0, 0)) - obj_factor * pivot_rhs_val)
}

pub fn get_optimal_primal(
	con_coef: &mut Compressed<f32>,
	con_rhs_coef: &mut Compressed<f32>,
	obj_coef: &mut Compressed<f32>,
	rhs_coef: &mut Compressed<f32>,
) {
	loop {
		// maximisation checks if all coefficients in the objective function are non-negative
		if obj_coef.values.iter().all(|&x| x >= 0.0) {
			println!("Optimal solution found, breaking.");
			break;
		}

		let (enter_idx, _) = get_enter_var(&obj_coef);
		let (leaving_idx, _) = get_leaving_var(&con_coef, &con_rhs_coef, enter_idx);

		// pivot the tableau
		pivot_coef(
			con_coef,
			con_rhs_coef,
			obj_coef,
			rhs_coef,
			enter_idx,
			leaving_idx,
		);
	}
}

pub fn print_comp_matrix(matrix: &Compressed<f32>) {
	for i in 0..matrix.rows {
		for j in 0..matrix.columns {
			print!("{:8.2} ", matrix.get((i, j)));
		}
		println!();
	}
}

pub fn print_conv_matrix(matrix: &Conventional<f32>) {
	for i in 0..matrix.rows {
		for j in 0..matrix.columns {
			print!("{:8.2} ", matrix[(i, j)]);
		}
		println!();
	}
}
